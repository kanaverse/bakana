import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as filter_module from "./cell_filtering.js";
import * as norm_module from "./rna_normalization.js";

/**
 * Results of per-gene variance modelling,
 * see [here](https://kanaverse.github.io/scran.js/ModelGeneVarResults.html) for details.
 *
 * @external ModelGeneVarResults
 */

/**
 * Feature selection is performed by modelling the per-gene variance and finding highly variable genes.
 * This wraps the [`modelGeneVariances`](https://kanaverse.github.io/scran.js/global.html#modelGeneVariances) function 
 * from [**scran.js**](https://github.com/kanaverse/scran.js).
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
            throw new Error("'filter' should be a CellFilteringState object");
        }
        this.#filter = filter;

        if (!(norm instanceof norm_module.RnaNormalizationState)) {
            throw new Error("'norm' should be an RnaNormalizationState object");
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

    valid() {
        return this.#norm.valid();
    }

    /**
     * @return {external:ModelGeneVarResults} Variance modelling results,
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
     *
     * @param {object} parameters - Parameter object, equivalent to the `feature_selection` property of the `parameters` of {@linkcode runAnalysis}.
     * @param {number} parameters.span - Value between 0 and 1 specifying the span for the LOWESS smoother.
     *
     * @return The object is updated with the new results.
     */
    compute(parameters) {
        let { span } = parameters;
        this.changed = false;
        
        if (this.#norm.changed || span != this.#parameters.span) {
            utils.freeCache(this.#cache.results);

            if (this.valid()) {
                let mat = this.#norm.fetchNormalizedMatrix();
                let block = this.#filter.fetchFilteredBlock();
                this.#cache.results = scran.modelGeneVariances(mat, { span: span, block: block });

                this.#cache.sorted_residuals = this.#cache.results.residuals().slice(); // a separate copy.
                this.#cache.sorted_residuals.sort();

                this.changed = true;
            }

            this.#parameters.span = span;
        }

        return;
    }
}
