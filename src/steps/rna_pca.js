import * as scran from "scran.js";
import * as utils from "./utils/general.js";
import * as filter_module from "./cell_filtering.js";
import * as norm_module from "./rna_normalization.js";
import * as feat_module from "./feature_selection.js";

export const step_name = "rna_pca";

/**
 * Results of running PCA on some input matrix,
 * see [here](https://kanaverse.github.io/scran.js/RunPCAResults.html) for details.
 *
 * @external RunPCAResults
 */

/**
 * This step performs a principal components analysis (PCA) to compact and denoise the data.
 * The resulting PCs can be used as input to various per-cell analyses like clustering and dimensionality reduction.
 * It wraps the [`runPca`](https://kanaverse.github.io/scran.js/global.html#runPca) function
 * from [**scran.js**](https://github.com/kanaverse/scran.js).
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class RnaPcaState { 
    #filter;
    #norm;
    #feat;
    #cache;
    #parameters;

    constructor(filter, norm, feat, parameters = null, cache = null) {
        if (!(filter instanceof filter_module.CellFilteringState)) {
            throw new Error("'filter' should be a CellFilteringState object");
        }
        this.#filter = filter;

        if (!(norm instanceof norm_module.RnaNormalizationState)) {
            throw new Error("'norm' should be an RnaNormalizationState object");
        }
        this.#norm = norm;

        if (!(feat instanceof feat_module.FeatureSelectionState)) {
            throw new Error("'feat' should be a FeatureSelectionState object");
        }
        this.#feat = feat;

        this.#parameters = (parameters === null ? {} : parameters);
        this.#cache = (cache === null ? {} : cache);
        this.changed = false;
    }

    free() {
        utils.freeCache(this.#cache.hvg_buffer);
        utils.freeCache(this.#cache.pcs);
    }

    /***************************
     ******** Getters **********
     ***************************/

    valid() {
        return this.#norm.valid();
    }

    /**
     * @return {external:RunPCAResults} Results of the PCA on the normalized gene expression values.
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
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     *
     * @param {object} parameters - Parameter object, equivalent to the `rna_pca` property of the `parameters` of {@linkcode runAnalysis}.
     * @param {number} [parameters.num_pcs] - Number of PCs to return.
     * @param {number} [parameters.num_hvgs] - Number of highly variable genes (see {@linkplain FeatureSelectionState}) to use in the PCA.
     * @param {string} [parameters.block_method] - Blocking method to use when dealing with multiple samples.
     * This can be one of:
     *
     * - `"none"`, in which case nothing is done using the sample information. 
     * - `"regress"`, where linear regression is applied to remove mean differences between samples.
     * - `"project"`, where samples are weighted so that they contribute equally regardless of the number of cells.
     *
     * @return The object is updated with the new results.
     */
    compute(parameters) {
        parameters = utils.defaultizeParameters(parameters, RnaPcaState.defaults());
        this.changed = false;

        // For back-compatibility.
        if (parameters.block_method == "weight") {
            parameters.block_method = "project";
        }

        if (this.#feat.changed || parameters.num_hvgs !== this.#parameters.num_hvgs) {
            if (this.valid()) {
                choose_hvgs(parameters.num_hvgs, this.#feat, this.#cache);
                this.changed = true;
            }
        }

        if (this.changed || this.#norm.changed || parameters.num_pcs !== this.#parameters.num_pcs || parameters.block_method !== this.#parameters.block_method) { 
            utils.freeCache(this.#cache.pcs);

            if (this.valid()) {
                let sub = this.#cache.hvg_buffer;
                let block = this.#filter.fetchFilteredBlock();
                var mat = this.#norm.fetchNormalizedMatrix();
                this.#cache.pcs = scran.runPca(mat, { features: sub, numberOfPCs: parameters.num_pcs, block: block, blockMethod: parameters.block_method });
                this.changed = true;
            }
        }

        this.#parameters = parameters;
        return;
    }

    /**
     * @return {object} Object containing default parameters,
     * see the `parameters` argument in {@linkcode RnaPcaState#compute compute} for details.
     */
    static defaults() {
        return {
            num_hvgs: 2000,
            num_pcs: 20,
            block_method: "none"
        };
    }
}

/**************************
 ******* Internals ********
 **************************/

function choose_hvgs(num_hvgs, feat, cache) {
    var sorted_resids = feat.fetchSortedResiduals();
    var sub = utils.allocateCachedArray(sorted_resids.length, "Uint8Array", cache, "hvg_buffer");

    if (num_hvgs < sorted_resids.length) {
        var threshold_at = sorted_resids[sorted_resids.length - num_hvgs];
        var unsorted_resids = feat.fetchResults().residuals({ copy: false });
        sub.array().forEach((element, index, array) => {
            array[index] = unsorted_resids[index] >= threshold_at;
        });
    } else {
        sub.fill(1);
    }

    return sub;
}
