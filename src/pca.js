import * as scran from "scran.js";
import * as utils from "./utils/general.js";
import * as putils from "./utils/pca.js";
import * as filter_module from "./cell_filtering.js";
import * as norm_module from "./normalization.js";
import * as feat_module from "./feature_selection.js";

export const step_name = "pca";

/**
 * This step performs a principal components analysis (PCA) to compact and denoise the data.
 * The resulting PCs can be used as input to various per-cell analyses like clustering and dimensionality reduction.
 * It wraps the `runPCA` function from [**scran.js**](https://github.com/jkanche/scran.js).
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class PcaState extends putils.PcaStateBase {
    #filter;
    #norm;
    #feat;
    #cache;
    #parameters;

    constructor(filter, norm, feat, parameters = null, cache = null) {
        super();

        if (!(filter instanceof filter_module.CellFilteringState)) {
            throw new Error("'filter' should be a State object from './cell_filtering.js'");
        }
        this.#filter = filter;

        if (!(norm instanceof norm_module.NormalizationState)) {
            throw new Error("'norm' should be a State object from './normalization.js'");
        }
        this.#norm = norm;

        if (!(feat instanceof feat_module.FeatureSelectionState)) {
            throw new Error("'feat' should be a State object from './feature_selection.js'");
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
        return true;
    }

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
     * @param {number} num_hvgs - Number of highly variable genes (see {@linkplain FeatureSelectionState}) to use in the PCA.
     * @param {string} block_method - Blocking method to use when dealing with multiple samples.
     * This can be `"none"`, `"block"` or `"weight"`.
     *
     * @return The object is updated with the new results.
     */
    compute(num_hvgs, num_pcs, block_method) {
        this.changed = false;

        if (this.#feat.changed || num_hvgs !== this.#parameters.num_hvgs) {
            choose_hvgs(num_hvgs, this.#feat, this.#cache);
            this.#parameters.num_hvgs = num_hvgs;
            this.changed = true;
        }

        if (this.changed || this.#norm.changed || num_pcs !== this.#parameters.num_pcs || block_method !== this.#parameters.block_method) { 
            let sub = this.#cache.hvg_buffer;
            let block = this.#filter.fetchFilteredBlock();
            if (block_method == "none") {
                block = null;
            }

            var mat = this.#norm.fetchNormalizedMatrix();
            utils.freeCache(this.#cache.pcs);
            this.#cache.pcs = scran.runPCA(mat, { features: sub, numberOfPCs: num_pcs, block: block, blockMethod: block_method });

            this.#parameters.num_pcs = num_pcs;
            this.#parameters.block_method = block_method;
            this.changed = true;
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

    /*************************
     ******** Saving *********
     *************************/

    serialize(handle) {
        let ghandle = handle.createGroup(step_name);

        {
            let phandle = ghandle.createGroup("parameters"); 
            phandle.writeDataSet("num_hvgs", "Int32", [], this.#parameters.num_hvgs);
            phandle.writeDataSet("num_pcs", "Int32", [], this.#parameters.num_pcs);
            phandle.writeDataSet("block_method", "String", [], this.#parameters.block_method);
        }

        {
            let rhandle = ghandle.createGroup("results");

            let ve = this.summary().var_exp;
            rhandle.writeDataSet("var_exp", "Float64", null, ve);

            let pcs = this.fetchPCs();
            rhandle.writeDataSet("pcs", "Float64", [pcs.num_obs, pcs.num_pcs], pcs.pcs); // remember, it's transposed.
        }
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
        var unsorted_resids = feat.fetchResiduals({ unsafe: true });
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
    let ghandle = handle.open(step_name);

    let parameters = {};
    {
        let phandle = ghandle.open("parameters"); 
        parameters = { 
            num_hvgs: phandle.open("num_hvgs", { load: true }).values[0],
            num_pcs: phandle.open("num_pcs", { load: true }).values[0]
        };

        // For back-compatibility.
        if ("block_method" in phandle.children) {
            parameters.block_method = phandle.open("block_method", { load: true }).values[0]
        } else {
            parameters.block_method = "none";
        }
    }

    let output;
    let cache = {};
    try {
        choose_hvgs(parameters.num_hvgs, feat, cache);

        let rhandle = ghandle.open("results");
        let var_exp = rhandle.open("var_exp", { load: true }).values;
        let pcs = rhandle.open("pcs", { load: true }).values;
        cache.pcs = new putils.PcaMimic(pcs, var_exp);

        output = new PcaState(filter, norm, feat, parameters, cache);
    } catch (e) {
        utils.freeCache(cache.hvg_buffer);
        utils.freeCache(cache.pcs);
        utils.freeCache(output);
        throw e;
    }

    return {
        state: output,
        parameters: parameters
    };
}
