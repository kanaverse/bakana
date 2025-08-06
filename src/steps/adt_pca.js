import * as scran from "scran.js";
import * as utils from "./utils/general.js";
import * as filter_module from "./cell_filtering.js";
import * as norm_module from "./adt_normalization.js";

export const step_name = "adt_pca";

/**
 * This step performs a principal components analysis (PCA) to compact and denoise ADT data.
 * The resulting PCs can be used as input to various per-cell analyses like clustering and dimensionality reduction.
 * It wraps the [`runPca`](https://kanaverse.github.io/scran.js/global.html#runPca) function
 * from [**scran.js**](https://github.com/kanaverse/scran.js).
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class AdtPcaState {
    #filter;
    #norm;
    #cache;
    #parameters;

    constructor(filter, norm, parameters = null, cache = null) {
        if (!(filter instanceof filter_module.CellFilteringState)) {
            throw new Error("'filter' should be a CellFilteringState object");
        }
        this.#filter = filter;

        if (!(norm instanceof norm_module.AdtNormalizationState)) {
            throw new Error("'norm' should be a AdtNormalizationState object");
        }
        this.#norm = norm;

        this.#parameters = (parameters === null ? {} : parameters);
        this.#cache = (cache === null ? {} : cache);
        this.changed = false;
    }

    free() {
        utils.freeCache(this.#cache.pcs);
    }

    /***************************
     ******** Getters **********
     ***************************/

    valid() {
        return this.#norm.valid();
    }

    /**
     * @return {external:RunPCAResults} Results of the PCA on the normalized ADT matrix,
     * available after running {@linkcode AdtPcaState#compute compute}.
     */
    fetchPCs() {
        return this.#cache.pcs;
    }

    /**
     * @return {object} Object containing the parameters.
     */
    fetchParameters() {
        return { ...this.#parameters }; // avoid pass-by-reference links.
    }

    /***************************
     ******** Compute **********
     ***************************/

    /**
     * @return {object} Object containing default parameters,
     * see the `parameters` argument in {@linkcode AdtPcaState#compute compute} for details.
     */
    static defaults() {
        return {
            num_pcs: 20,
            block_method: "none"
        };
    }

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     *
     * @param {object} parameters - Parameter object, equivalent to the `adt_pca` property of the `parameters` of {@linkcode runAnalysis}.
     * @param {number} [parameters.num_pcs] - Number of PCs to return.
     * @param {string} [parameters.block_method] - Blocking method to use when dealing with multiple samples.
     * This can be `"none"`, `"regress"` or `"project"`, see comments in {@linkplain RnaPcaState}.
     *
     * @return The object is updated with the new results.
     */
    compute(parameters) {
        parameters = utils.defaultizeParameters(parameters, AdtPcaState.defaults());
        this.changed = false;

        // For back-compatibility:
        if (parameters.block_method == "weight") {
            parameters.block_method = "project";
        }

        if (this.#norm.changed || parameters.num_pcs !== this.#parameters.num_pcs || parameters.block_method !== this.#parameters.block_method) { 
            if (this.valid()) {
                let block = this.#filter.fetchFilteredBlock();
                var mat = this.#norm.fetchNormalizedMatrix();
                utils.freeCache(this.#cache.pcs);
                this.#cache.pcs = scran.runPca(mat, { numberOfPCs: parameters.num_pcs, block: block, blockMethod: parameters.block_method });

                this.changed = true;
            }

        }

        this.#parameters = parameters;
        return;
    }
}
