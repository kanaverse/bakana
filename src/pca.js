import * as scran from "scran.js";
import * as utils from "./utils/general.js";
import * as qc_module from "./quality_control.js";
import * as norm_module from "./normalization.js";
import * as feat_module from "./feature_selection.js";

/**
 * This step performs a principal components analysis (PCA) to compact and denoise the data.
 * The resulting PCs can be used as input to various per-cell analyses like clustering and dimensionality reduction.
 * It wraps the `runPCA` function from [**scran.js**](https://github.com/jkanche/scran.js).
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class PcaState {
    #qc;
    #norm;
    #feat;
    #cache;
    #parameters;

    constructor(qc, norm, feat, parameters = null, cache = null) {
        if (!(qc instanceof qc_module.QualityControlState)) {
            throw new Error("'qc' should be a State object from './quality_control.js'");
        }
        this.#qc = qc;

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
        utils.freeCache(this.#cache.corrected);
    }

    /***************************
     ******** Getters **********
     ***************************/

    fetchPCs({ original = false } = {}) {
        let pcs;
        if (!original && this.#parameters.block_method == "mnn") {
            pcs = this.#cache.corrected;
        } else {
            pcs = this.#cache.pcs.principalComponents({ copy: "view" });
        }
        return {
            "pcs": pcs,
            "num_pcs": this.#parameters.num_pcs,
            "num_obs": pcs.length / this.#parameters.num_pcs
        };
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
     * This can be `"none"`, `"block"` or `"mnn"`.
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

            let block = this.#qc.fetchFilteredBlock();
            let block_type = "block";
            if (block_method == "none") {
                block = null;
            } else if (block_method == "mnn") {
                block_type = "weight";
            }

            var mat = this.#norm.fetchNormalizedMatrix();

            utils.freeCache(this.#cache.pcs);
            this.#cache.pcs = scran.runPCA(mat, { features: sub, numberOfPCs: num_pcs, block: block, blockMethod: block_type });

            if (block_method == "mnn") {
                let pcs = this.#cache.pcs.principalComponents({ copy:"view" });
                let corrected = utils.allocateCachedArray(pcs.length, "Float64Array", this.#cache, "corrected");
                scran.mnnCorrect(this.#cache.pcs, block, { buffer: corrected });
            }

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
        var pca_output = this.#cache.pcs;
        var var_exp = pca_output.varianceExplained();
        var total_var = pca_output.totalVariance();
        var_exp.forEach((x, i) => {
            var_exp[i] = x/total_var;
        });
        return { "var_exp": var_exp };
    }

    /*************************
     ******** Saving *********
     *************************/

    serialize(handle) {
        let ghandle = handle.createGroup("pca");

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

            let pcs = this.fetchPCs({ original: true });
            rhandle.writeDataSet("pcs", "Float64", [pcs.num_obs, pcs.num_pcs], pcs.pcs); // remember, it's transposed.

            if (this.#parameters.block_method == "mnn") {
                let corrected = this.#cache.corrected;
                rhandle.writeDataSet("corrected", "Float64", [pcs.num_obs, pcs.num_pcs], corrected); 
            }
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

class PCAMimic { 
    constructor(pcs, var_exp) {
        this.var_exp = var_exp;
        this.pcs = scran.createFloat64WasmArray(pcs.length);
        this.pcs.set(pcs);
    }

    principalComponents({ copy }) {
        return utils.mimicGetter(this.pcs, copy);
    }

    varianceExplained({ copy = true } = {}) {
        return utils.mimicGetter(this.var_exp, copy);
    }

    totalVariance () {
        return 1;
    }

    free() {
        this.pcs.free();
    }
}

export function unserialize(handle, qc, norm, feat) {
    let ghandle = handle.open("pca");

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

    let cache = {};
    choose_hvgs(parameters.num_hvgs, feat, cache);

    {
        let rhandle = ghandle.open("results");
        let var_exp = rhandle.open("var_exp", { load: true }).values;
        let pcs = rhandle.open("pcs", { load: true }).values;
        cache.pcs = new PCAMimic(pcs, var_exp);

        if (parameters.block_method == "mnn") {
            let corrected = rhandle.open("corrected", { load: true }).values;
            let corbuffer = utils.allocateCachedArray(corrected.length, "Float64Array", cache, "corrected");
            corbuffer.set(corrected);
        }
    }

    return {
        state: new PcaState(qc, norm, feat, parameters, cache),
        parameters: { ...parameters }
    };
}
