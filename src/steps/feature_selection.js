import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as filter_module from "./cell_filtering.js";
import * as norm_module from "./normalization.js";

/**
 * Feature selection is performed by modelling the per-gene variance and finding highly variable genes.
 * This wraps the `modelGeneVar` function from [**scran.js**](https://github.com/jkanche/scran.js).
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class FeatureSelectionState {
    #filter;
    #norm;
    #cache;
    #parameters;

    constructor(filter, norm, parameters = null, cache = null) {
        if (!(filter instanceof filter_module.CellFilteringState)) {
            throw new Error("'filter' should be a State object from './cell_filtering.js'");
        }
        this.#filter = filter;

        if (!(norm instanceof norm_module.NormalizationState)) {
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

    fetchParameters() {
        return { ...this.#parameters }; // avoid pass-by-reference activity.
    }

    /***************************
     ******** Compute **********
     ***************************/
 
    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     * Each argument is taken from the property of the same name in the `feature_selection` property of the `parameters` of {@linkcode runAnalysis}.
     *
     * @param {number} span - Value between 0 and 1 specifying the span for the LOWESS smoother.
     *
     * @return The object is updated with the new results.
     */
    compute(span) {
        this.changed = false;
        
        if (this.#norm.changed || span != this.#parameters.span) {
            utils.freeCache(this.#cache.results);

            let mat = this.#norm.fetchNormalizedMatrix();
            let block = this.#filter.fetchFilteredBlock();
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
        return {
            "means": this.#cache.results.means({ copy: copy }),
            "vars": this.#cache.results.variances({ copy: copy }),
            "fitted": this.#cache.results.fitted({ copy: copy }),
            "resids": this.#cache.results.residuals({copy: copy })
        };
    }

    /**
     * Obtain a summary of the state, typically for display on a UI like **kana**.
     *
     * @return An object containing:
     *
     * - `means`: a Float64Array containing the mean log-expression for each gene.
     * - `vars`: a Float64Array containing the variance in log-expression for each gene.
     * - `fitted`: a Float64Array containing the fitted value of the mean-variance trend for each gene.
     * - `resids`: a Float64Array containing the residuals from the mean-variance trend for each gene.
     */
    summary() {
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
            let res = this.#format_results({ copy: "hdf5" }); 
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

export function unserialize(handle, permuter, filter, norm) {
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

    return new FeatureSelectionState(filter, norm, parameters, cache);
}

