import * as scran from "scran.js";
import * as utils from "./utils/general.js";
import * as rna_pca_module from "./rna_pca.js";
import * as adt_pca_module from "./adt_pca.js";

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
        if (!(pca_states.RNA instanceof rna_pca_module.RnaPcaState)) {
            throw new Error("'pca_states.RNA' should be an RnaPcaState object");
        }
        if (!(pca_states.ADT instanceof adt_pca_module.AdtPcaState)) {
            throw new Error("'pca_states.ADT' should be an AdtPcaState object");
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

    /**
     * @return {Float64WasmArray} Buffer containing the combined embeddings as a column-major dense matrix,
     * where the rows are the dimensions and the columns are the cells.
     * This is available after running {@linkcode CombineEmbeddingsState#compute compute}.
     */
    fetchCombined() {
        return this.#cache.combined_buffer;
    }

    /**
     * @return {number} Number of cells in {@linkcode CombineEmbeddingsState#fetchCombined fetchCombined},
     * available after running {@linkcode CombineEmbeddingsState#compute compute}.
     */
    fetchNumberOfCells() {
        return this.#cache.num_cells;
    }

    /**
     * @return {number} Number of dimensions in {@linkcode CombineEmbeddingsState#fetchCombined fetchCombined},
     * available after running {@linkcode CombineEmbeddingsState#compute compute}.
     */
    fetchNumberOfDimensions() {
        return this.#cache.total_dims;
    }

    /**
     * @return {object} Object containing the parameters.
     */
    fetchParameters() {
        // Avoid any pass-by-reference activity.
        let out = { ...this.#parameters };
        if (out.weights !== null) {
            out.weights = { ...out.weights };
        }
        return out;
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
        cache.combined_buffer = upstream.principalComponents({ copy: "view" }).view();
        cache.num_cells = upstream.numberOfCells();
        cache.total_dims = upstream.numberOfPCs();
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

        for (const v of Object.values(this.#pca_states)) {
            if (v.changed) {
                this.changed = true;
                break;
            }
        }

        let to_use = utils.findValidUpstreamStates(this.#pca_states);
        let needs_combining = to_use.length > 1;

        if (needs_combining) {
            if (this.changed || approximate !== this.#parameters.approximate || utils.changedParameters(weights, this.#parameters.weights)) {
                let weight_arr = null;

                if (weights !== null) {
                    weight_arr = [];
                    let has_pos_weight = [];

                    for (const x of to_use) {
                        if (!(x in weights)) {
                            throw new Error("no weight specified for '" + x + "'");
                        }
                        if (weights[x] > 0) {
                            weight_arr.push(weights[x]);
                            has_pos_weight.push(x);
                        }
                    }

                    to_use = has_pos_weight;
                }

                if (to_use.length == 1) {
                    // If only one modality has positive weight,
                    // we can skip the calculation of combined embeddings.
                    let pcs = this.#pca_states[to_use[0]].fetchPCs();
                    this.constructor.createPcsView(this.#cache, pcs);

                } else {
                    let collected = [];
                    let total = 0;
                    let ncells = null;

                    for (const k of to_use) {
                        let curpcs = this.#pca_states[k].fetchPCs();
                        collected.push(curpcs.principalComponents({ copy: "view" }));
                        if (ncells == null) {
                            ncells = curpcs.numberOfCells();
                        } else if (ncells !== curpcs.numberOfCells()) {
                            throw new Error("number of cells should be consistent across all embeddings");
                        }
                        total += curpcs.numberOfPCs();
                    }

                    let buffer = utils.allocateCachedArray(ncells * total, "Float64Array", this.#cache, "combined_buffer");
                    scran.scaleByNeighbors(collected, ncells, { buffer: buffer, weights: weight_arr, approximate: approximate });
                    this.#cache.num_cells = ncells;
                    this.#cache.total_dims = total;
                }

                this.changed = true;
            }
        } else {
            if (this.changed) {
                // If there's only one embedding, we shouldn't respond to changes
                // in parameters, because they won't have any effect.
                let pcs = this.#pca_states[to_use[0]].fetchPCs();
                this.constructor.createPcsView(this.#cache, pcs);
            }
        }

        // Updating all parameters anyway. This requires us to take ownership
        // of 'weights' to avoid pass-by-reference shenanigans.
        if (weights !== null) {
            weights = { ...weights };
        }
        this.#parameters.weights = weights;
        this.#parameters.approximate = approximate;
        return;
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
            let pcs = this.fetchCombined();
            if (pcs.owner === null) {
                // If it's not a view, we save it; otherwise we assume
                // that we can recover it from the upstream PCA states.
                rhandle.writeDataSet(
                    "combined", 
                    "Float64", 
                    [this.fetchNumberOfCells(), this.fetchNumberOfDimensions()], // remember, it's transposed.
                    pcs
                ); 
            }
        }
    }
}

export function unserialize(handle, pca_states) {
    let cache = {};
    let parameters = CombineEmbeddingsState.defaults();
    let output;

    try {
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

            let rhandle = ghandle.open("results");

            if ("combined" in rhandle.children) {
                let phandle = rhandle.open("combined", { load: true });
                cache.num_cells = phandle.shape[0];
                cache.total_dims = phandle.shape[1];

                let vals = phandle.values;
                cache.combined_buffer = scran.createFloat64WasmArray(vals.length);
                cache.combined_buffer.set(vals);
            }
        }

        if (!("combined_buffer" in cache)) {
            // This only happens if there was only one upstream PCA state; in which case, 
            // we figure out which upstream PCA state contains the PC vector
            // and create a view on it so that our fetchPCs() works properly.
            // (v1 and earlier also implicitly falls in this category.)

            let to_use = utils.findValidUpstreamStates(pca_states);

            if (to_use.length > 1 && parameters.weights !== null) {
                let has_nonzero_weight = [];
                for (const k of to_use) {
                    if (parameters.weights[k] > 0) {
                        has_nonzero_weight.push(k);
                    }
                }
                to_use = has_nonzero_weight;
            }

            if (to_use.length != 1) {
                throw new Error("only one upstream PCA state should be valid with non-zero weight if 'combined' is not available");
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

    return output;
}
