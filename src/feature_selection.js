import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as qc_module from "./quality_control.js";
import * as norm_module from "./normalization.js";

/**
 * Feature selection is performed by modelling the per-gene variance and finding highly variable genes.
 * This wraps the `modelGeneVar` function from [**scran.js**](https://github.com/jkanche/scran.js).
 *
 * The parameters in {@linkcode runAnalysis} should be an object containing:
 *
 * - `span`: a number between 0 and 1 specifying the span for the LOWESS smoother.
 *
 * Calling the **`results()`** method for the relevant state instance will return an object containing:
 *
 * - `means`: a `Float64Array` containing the mean log-expression for each gene.
 * - `vars`: a `Float64Array` containing the variance in log-expression for each gene.
 * - `fitted`: a `Float64Array` containing the fitted value of the mean-variance trend for each gene.
 * - `resids`: a `Float64Array` containing the residuals from the mean-variance trend for each gene.
 * 
 * Methods not documented here are not part of the stable API and should not be used by applications.
 *
 * @namespace feature_selection
 */

export class State {
    #qc;
    #norm;
    #cache;
    #parameters;

    constructor(qc, norm, parameters = null, cache = null) {
        if (!(qc instanceof qc_module.State)) {
            throw new Error("'qc' should be a State object from './quality_control.js'");
        }
        this.#qc = qc;

        if (!(norm instanceof norm_module.State)) {
            throw new Error("'norm' should be a State object from './normalization.js'");
        }
        this.#norm = norm;

        this.#parameters = (parameters === null ? {} : parameters);
        this.#cache = (cache === null ? {} : cache);
        this.changed = false;
    }

    free() {
        utils.freeCache(this.#cache.matrix);
    }

    /***************************
     ******** Getters **********
     ***************************/

    fetchSortedResiduals() {
        return this.#cache.sorted_residuals;
    }

    fetchResiduals({ unsafe = false } = {}) {
        return this.#cache.results.residuals({ copy: !unsafe });
    }

    /***************************
     ******** Compute **********
     ***************************/

    compute(span) {
        this.changed = false;
        
        if (this.#norm.changed || span != this.#parameters.span) {
            utils.freeCache(this.#cache.results);

            let mat = this.#norm.fetchNormalizedMatrix();
            let block = this.#qc.fetchFilteredBlock();
            this.#cache.results = scran.modelGeneVar(mat, { span: span, block: block });

            this.#cache.sorted_residuals = this.#cache.results.residuals().slice(); // a separate copy.
            this.#cache.sorted_residuals.sort();

            this.#parameters.span = span;
            this.changed = true;
        }

        return;
    }

    /***************************
     ******** Results **********
     ***************************/

    #format_results({ copy = true } = {}) {
        copy = (copy ? true : "view");
        return {
            "means": this.#cache.results.means({ copy: copy }),
            "vars": this.#cache.results.variances({ copy: copy }),
            "fitted": this.#cache.results.fitted({ copy: copy }),
            "resids": this.#cache.results.residuals({copy: copy })
        };
    }

    results() {
        return this.#format_results();
    }

    /*************************
     ******** Saving *********
     *************************/

    serialize(handle) {
        let ghandle = handle.createGroup("feature_selection");

        {
            let phandle = ghandle.createGroup("parameters"); 
            phandle.writeDataSet("span", "Float64", [], this.#parameters.span);
        }

        {
            let res = this.#format_results({ copy: false });
            let rhandle = ghandle.createGroup("results"); 
            for (const [k, v] of Object.entries(res)) {
                rhandle.writeDataSet(k, "Float64", null, v);
            }
        }
    }
}

/**************************
 ******** Loading *********
 **************************/

class ModelGeneVarMimic {
    constructor(means, vars, fitted, resids) {
        this.means_ = means;
        this.vars_ = vars;
        this.fitted_ = fitted;
        this.resids_ = resids;
    }

    means({copy}) {
        return utils.mimicGetter(this.means_, copy);
    }

    variances({copy}) {
        return utils.mimicGetter(this.vars_, copy);
    }

    fitted({copy}) {
        return utils.mimicGetter(this.fitted_, copy);
    }

    residuals({copy}) {
        return utils.mimicGetter(this.resids_, copy);
    }

    free() {}
}

export function unserialize(handle, permuter, qc, norm) {
    let ghandle = handle.open("feature_selection");

    let parameters;
    {
        let phandle = ghandle.open("parameters");
        parameters = {
            span: phandle.open("span", { load: true }).values[0]
        };
    }

    let cache = {};
    {
        let rhandle = ghandle.open("results");
        let reloaded = {};

        // Possibly permuting it to match the new permutation order;
        // see 'unserialize' in 'inputs.js'.
        for (const key of [ "means", "vars", "fitted", "resids" ]) {
            let value = rhandle.open(key, { load: true }).values;
            permuter(value);
            reloaded[key] = value;
        }

        cache.results = new ModelGeneVarMimic(reloaded.means, reloaded.vars, reloaded.fitted, reloaded.resids);
    }

    cache.sorted_residuals = cache.results.residuals({ copy: true });
    cache.sorted_residuals.sort();

    return { 
        state: new State(qc, norm, parameters, cache),
        parameters: { ...parameters }
    };
}

