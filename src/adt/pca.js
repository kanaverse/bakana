import * as scran from "scran.js";
import * as utils from "../utils/general.js";
import * as putils from "../utils/pca.js";
import * as filter_module from "../cell_filtering.js";
import * as norm_module from "./normalization.js";

export const step_name = "adt_pca";

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
        super();

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

    fetchPCs() {
        if (this.valid()) {
            return putils.formatPCs(this.#cache.pcs);
        } else {
            return null;
        }
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
            if (this.valid()) {
                let block = this.#filter.fetchFilteredBlock();
                if (block_method == "none") {
                    block = null;
                }

                var mat = this.#norm.fetchNormalizedMatrix();
                utils.freeCache(this.#cache.pcs);
                this.#cache.pcs = scran.runPCA(mat, { numberOfPCs: num_pcs, block: block, blockMethod: block_method });

                this.changed = true;
            }

            this.#parameters.num_pcs = num_pcs;
            this.#parameters.block_method = block_method;
        }

        return;
    }

    static defaults() {
        return {
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
        if (this.valid()) {
            return putils.formatSummary(this.#cache.pcs);
        } else {
            return null;
        }
    }

    /*************************
     ******** Saving *********
     *************************/

    serialize(handle) {
        let ghandle = handle.createGroup(step_name);

        {
            let phandle = ghandle.createGroup("parameters"); 
            phandle.writeDataSet("num_pcs", "Int32", [], this.#parameters.num_pcs);
            phandle.writeDataSet("block_method", "String", [], this.#parameters.block_method);
        }

        {
            let rhandle = ghandle.createGroup("results");

            if (this.valid()) {
                let ve = this.summary().var_exp;
                rhandle.writeDataSet("var_exp", "Float64", null, ve);

                let pcs = this.fetchPCs();
                rhandle.writeDataSet("pcs", "Float64", [pcs.num_obs, pcs.num_pcs], pcs.pcs); // remember, it's transposed.
            }
        }
    }
}

/**************************
 ******** Loading *********
 **************************/

export function unserialize(handle, filter, norm) {
    let cache = {};
    let parameters = AdtPcaState.defaults();
    let output;

    if (step_name in handle.children) {
        let ghandle = handle.open(step_name);

        let phandle = ghandle.open("parameters"); 
        parameters.num_pcs = phandle.open("num_pcs", { load: true }).values[0];
        parameters.block_method = phandle.open("block_method", { load: true }).values[0];

        try {
            let rhandle = ghandle.open("results");

            if ("var_exp" in rhandle.children) {
                let var_exp = rhandle.open("var_exp", { load: true }).values;
                let pcs = rhandle.open("pcs", { load: true }).values;
                cache.pcs = new putils.PcaMimic(pcs, var_exp);
            }

            output = new AdtPcaState(filter, norm, parameters, cache);
        } catch (e) {
            utils.freeCache(cache.pcs);
            utils.freeCache(output);
            throw e;
        }
    } else {
        output = new AdtPcaState(filter, norm);
    }

    return {
        state: output,
        parameters: { ...parameters } // make a copy to avoid pass-by-reference links with state's internal parameters
    };
}
