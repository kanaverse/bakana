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

    /**
     * @return {ModelGeneVarResults} Variance modelling results,
     * available after running {@linkcode FeatureSelectionState#compute compute}.
     */
    fetchResults() {
        return this.#cache.results;
    }

    /**
     * @return {Float64Array} Array of length equal to the number of genes,
     * containing the sorted residuals after fitting a mean-dependent trend to the variances.
     * Available after running {@linkcode FeatureSelectionState#compute compute}.
     */
    fetchSortedResiduals() {
        return this.#cache.sorted_residuals;
    }

    /**
     * @return {object} Object containing the parameters.
     */
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
            let res = {
                "means": this.#cache.results.means({ copy: "hdf5" }),
                "vars": this.#cache.results.variances({ copy: "hdf5" }),
                "fitted": this.#cache.results.fitted({ copy: "hdf5" }),
                "resids": this.#cache.results.residuals({copy: "hdf5" })
            };

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

        cache.results = scran.emptyModelGeneVarResults(reloaded.means.length, 1);
        cache.results.means({ copy : false }).set(reloaded.means);
        cache.results.variances({ copy : false }).set(reloaded.variances);
        cache.results.fitted({ copy : false }).set(reloaded.fitted);
        cache.results.residuals({ copy : false }).set(reloaded.resids);
    }

    cache.sorted_residuals = cache.results.residuals({ copy: true });
    cache.sorted_residuals.sort();

    return new FeatureSelectionState(filter, norm, parameters, cache);
}

