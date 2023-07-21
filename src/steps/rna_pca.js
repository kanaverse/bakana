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
     * @param {number} parameters.num_pcs - Number of PCs to return.
     * @param {number} parameters.num_hvgs - Number of highly variable genes (see {@linkplain FeatureSelectionState}) to use in the PCA.
     * @param {string} parameters.block_method - Blocking method to use when dealing with multiple samples.
     * This can be one of:
     *
     * - `"none"`, in which case nothing is done using the sample information. 
     * - `"regress"`, where linear regression is applied to remove mean differences between samples.
     * - `"project"`, where samples are weighted so that they contribute equally regardless of the number of cells.
     *
     * @return The object is updated with the new results.
     */
    compute(parameters) {
        let { num_hvgs, num_pcs, block_method } = parameters;
        if (block_method == "weight") {
            block_method = "project";
        }
        this.changed = false;

        if (this.#feat.changed || num_hvgs !== this.#parameters.num_hvgs) {
            if (this.valid()) {
                choose_hvgs(num_hvgs, this.#feat, this.#cache);
                this.changed = true;
            }

            this.#parameters.num_hvgs = num_hvgs;
        }

        if (this.changed || this.#norm.changed || num_pcs !== this.#parameters.num_pcs || block_method !== this.#parameters.block_method) { 
            utils.freeCache(this.#cache.pcs);

            if (this.valid()) {
                let sub = this.#cache.hvg_buffer;
                let block = this.#filter.fetchFilteredBlock();
                var mat = this.#norm.fetchNormalizedMatrix();
                this.#cache.pcs = scran.runPca(mat, { features: sub, numberOfPCs: num_pcs, block: block, blockMethod: block_method });
                this.changed = true;
            }

            this.#parameters.num_pcs = num_pcs;
            this.#parameters.block_method = block_method;
        }

        return;
    }

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

/**************************
 ******** Loading *********
 **************************/

export function unserialize(handle, filter, norm, feat) {
    let ghandle = handle.open("rna_pca" in handle.children ? "rna_pca" : "pca");

    let parameters = {};
    {
        let phandle = ghandle.open("parameters"); 
        parameters = { 
            num_hvgs: phandle.open("num_hvgs", { load: true }).values[0],
            num_pcs: phandle.open("num_pcs", { load: true }).values[0]
        };

        // For back-compatibility.
        if ("block_method" in phandle.children) {
            parameters.block_method = phandle.open("block_method", { load: true }).values[0];
            if (parameters.block_method == "mnn") {
                parameters.block_method = "project";
            }
        } else {
            parameters.block_method = "none";
        }
    }

    let output;
    let cache = {};
    try {
        if (feat.valid()) {
            choose_hvgs(parameters.num_hvgs, feat, cache);

            let rhandle = ghandle.open("results");
            if ("pcs" in rhandle.children) {
                let pcs_handle = rhandle.open("pcs", { load: true });
                let pcs = pcs_handle.values;
                let var_exp = rhandle.open("var_exp", { load: true }).values;

                cache.pcs = scran.emptyRunPCAResults(pcs_handle.shape[0], pcs_handle.shape[1]);
                cache.pcs.principalComponents({ fillable: true }).set(pcs);
                cache.pcs.varianceExplained({ fillable: true }).set(var_exp);
                cache.pcs.setTotalVariance(1); // because the file only stores proportions.
            }
        }

        output = new RnaPcaState(filter, norm, feat, parameters, cache);
    } catch (e) {
        utils.freeCache(cache.hvg_buffer);
        utils.freeCache(cache.pcs);
        utils.freeCache(output);
        throw e;
    }

    return output;
}
