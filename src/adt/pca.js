import * as scran from "scran.js";
import * as utils from "./utils/general.js";
import * as putils from "./utils/pca.js";
import * as filter_module from "./cell_filtering.js";
import * as norm_module from "./normalization.js";
import * as feat_module from "./feature_selection.js";

/**
 * This step performs a principal components analysis (PCA) to compact and denoise ADT data.
 * The resulting PCs can be used as input to various per-cell analyses like clustering and dimensionality reduction.
 * It wraps the `runPCA` function from [**scran.js**](https://github.com/jkanche/scran.js).
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class AdtPcaState extends putils.PcaStateBase {
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

    fetchPCs() {
        return putils.formatPCs(this.#cache.pcs);
    }

    /***************************
     ******** Compute **********
     ***************************/

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     * Each argument is taken from the property of the same name in the `pca` property of the `parameters` of {@linkcode runAnalysis}.
     *
     * @param {number} num_pcs - Number of PCs to return.
     * @param {string} block_method - Blocking method to use when dealing with multiple samples.
     * This can be `"none"`, `"block"` or `"weight"`.
     *
     * @return The object is updated with the new results.
     */
    compute(num_pcs, block_method) {
        this.changed = false;

        if (this.#norm.changed || num_pcs !== this.#parameters.num_pcs || block_method !== this.#parameters.block_method) { 
            let block = this.#filter.fetchFilteredBlock();
            if (block_method == "none") {
                block = null;
            }

            var mat = this.#norm.fetchNormalizedMatrix();
            utils.freeCache(this.#cache.pcs);
            this.#cache.pcs = scran.runPCA(mat, { numberOfPCs: num_pcs, block: block, blockMethod: block_method });

            this.#parameters.num_pcs = num_pcs;
            this.#parameters.block_method = block_method;
            this.changed = true;
        }

        return;
    }

    /***************************
     ******** Results **********
     ***************************/

    /**
     * Obtain a summary of the state, typically for display on a UI like **kana**.
     *
     * @return An object containing:
     *
     * - `var_exp`: a `Float64Array` of length equal to `num_pcs`, containing the proportion of variance explained for each successive PC.
     */
    summary() {
        return putils.formatSummary(this.#cache.pcs);
    }
}
