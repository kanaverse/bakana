import * as scran from "scran.js";
import * as utils from "./utils/general.js";
import * as putils from "./utils/pca.js";

export const step_name = "combine_embeddings";

/**
 * This step combines multiple embeddings from different modalities into a single matrix for downstream analysis.
 * It wraps the `scaleByNeighbors` function from [**scran.js**](https://jkanche.com/scran.js).
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class CombineEmbeddingsState {
    #pca_states;
    #parameters;
    #cache;

    constructor(pca_states, parameters = null, cache = null) {
        for (const x of Object.values(pca_states)) {
            if (!(x instanceof putils.PcaStateBase)) {
                throw new Error("each entry of 'pc_states' should be a PcaStateBase object");
            }
        }
        this.#pca_states = pca_states;
        this.#parameters = (parameters === null ? {} : parameters);
        this.#cache = (cache === null ? {} : cache);
        this.changed = false;
    }

    free() {
        utils.freeCache(this.#cache.combined_buffer);
    }

    /***************************
     ******** Getters **********
     ***************************/

    fetchPCs() {
        return {
            "pcs": this.#cache.combined_buffer,
            "num_obs": this.#cache.num_cells,
            "num_pcs": this.#cache.total_dims
        };
    }

    /***************************
     ******** Compute **********
     ***************************/

    static defaults() {
        return { 
            weights: null,
            approximate: true
        };
    }

    static createPcsView(cache, upstream) {
        utils.freeCache(cache.combined_buffer);
        cache.combined_buffer = upstream.pcs.view();
        cache.num_cells = upstream.num_obs;
        cache.total_dims = upstream.num_pcs;
    }

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     *
     * @param {?Object} weights - Object containing the weights for each modality.
     * Keys should be the name of the modality (e.g., `"RNA"`, `"ADT"`) while values should be non-negative numbers.
     * Alternatively, if `null`, all modalities are given equal weight.
     * @param {boolean} approximate - Whether an approximate nearest neighbor search should be used by `scaleByNeighbors`.
     *
     * @return The object is updated with new results.
     */
    compute(weights, approximate) {
        this.changed = false;

        for (const [k, v] of Object.entries(this.#pca_states)) {
            if (v.changed) { // include possible changes from valid to invalid.
                this.changed = true;
            }
        }
        let to_use = utils.findValidUpstreamStates(this.#pca_states, "PCA");
        let needs_combining = to_use.length > 1;

        if (this.changed || approximate !== this.#parameters.approximate || utils.changedParameters(weights, this.#parameters.weights)) {
            if (needs_combining) {
                let collected = [];
                let total = 0;
                let ncells = null;

                for (const k of to_use) {
                    let curpcs = this.#pca_states[k].fetchPCs();
                    collected.push(curpcs.pcs);
                    if (ncells == null) {
                        ncells = curpcs.num_obs;
                    } else if (ncells !== curpcs.num_obs) {
                        throw new Error("number of cells should be consistent across all embeddings");
                    }
                    total += curpcs.num_pcs;
                }

                let weight_arr = null;
                if (weights !== null) {
                    weight_arr = [];
                    for (const x of to_use) {
                        if (!(x in weights)) {
                            throw new Error("no weight specified for '" + x + "'");
                        }
                        weight_arr.push(weights[x]);
                    }
                }

                let buffer = utils.allocateCachedArray(ncells * total, "Float64Array", this.#cache, "combined_buffer");
                scran.scaleByNeighbors(collected, ncells, { buffer: buffer, weights: weight_arr, approximate: approximate });
                this.#cache.num_cells = ncells;
                this.#cache.total_dims = total;
                this.changed = true;
            }
        }

        if (this.changed) {
            if (!needs_combining) {
                // If there's only one embedding, we shouldn't respond to changes
                // in parameters, because they won't have any effect.
                let pcs = this.#pca_states[to_use[0]].fetchPCs();
                this.constructor.createPcsView(this.#cache, pcs);
            }
        }

        // Updating all parameters anyway.
        this.#parameters.weights = weights;
        this.#parameters.approximate = approximate;
        return;
    }

    /***************************
     ******** Results **********
     ***************************/

    summary() {
        return {};
    }

    /*************************
     ******** Saving *********
     *************************/

    serialize(handle) {
        let ghandle = handle.createGroup(step_name);

        {
            let phandle = ghandle.createGroup("parameters"); 
            phandle.writeDataSet("approximate", "Uint8", [], this.#parameters.approximate);
            let whandle = phandle.createGroup("weights");
            if (this.#parameters.weights !== null) {
                for (const [k, v] of Object.entries(this.#parameters.weights)) {
                    whandle.writeDataSet(k, "Float64", [], v);
                }
            }
        }

        {
            let rhandle = ghandle.createGroup("results");
            let pcs = this.fetchPCs();
            if (pcs.pcs.owner === null) {
                // If it's not a view, we save it; otherwise we assume
                // that we can recover it from the upstream PCA states.
                rhandle.writeDataSet("combined", "Float64", [pcs.num_obs, pcs.num_pcs], pcs.pcs); // remember, it's transposed.
            }
        }
    }
}

export function unserialize(handle, pca_states) {
    let cache = {};
    let parameters = CombineEmbeddingsState.defaults();
    let output;

    if (step_name in handle.children) {
        let ghandle = handle.open(step_name);

        {
            let phandle = ghandle.open("parameters");
            parameters.approximate = phandle.open("approximate", { load: true }).values[0] > 0;

            let whandle = phandle.open("weights");
            let keys = Object.keys(whandle.children);

            // If it's empty, we just use the default of null.
            if (keys.length) {
                parameters.weights = {};
                for (const k of keys) {
                    parameters.weights[k] = whandle.open(k, { load: true }).values[0];
                }
            }
        }

        try {
            let rhandle = ghandle.open("results");

            if ("combined" in rhandle.children) {
                let phandle = rhandle.open("combined", { load: true });
                cache.num_cells = phandle.shape[0];
                cache.total_dims = phandle.shape[1];

                let vals = phandle.values;
                cache.combined_buffer = scran.createFloat64WasmArray(vals.length);
                cache.combined_buffer.set(vals);
            } else {
                // Creating a view from the valid upstream PCA state.
                let to_use = utils.findValidUpstreamStates(pca_states, "PCA");
                if (to_use.length != 1) {
                    throw new Error("only one upstream PCA state should be valid if 'discards' is not available");
                }
                let pcs = pca_states[to_use[0]].fetchPCs();
                CombineEmbeddingsState.createPcsView(cache, pcs);
            }

            output = new CombineEmbeddingsState(pca_states, parameters, cache);
        } catch (e) {
            utils.freeCache(cache.combined_buffer);
            utils.freeCache(output);
            throw e;
        }
    } else {
        // Fallback for v1.
        output = new CombineEmbeddingsState(pca_states);
    }

    return {
        state: output,

        // Make a copy to avoid pass-by-reference links with state's internal
        // parameters. Arguments include nested dicts so we have to do a deep
        // copy via stringify.
        parameters: JSON.parse(JSON.stringify(parameters)) 
    };
}
