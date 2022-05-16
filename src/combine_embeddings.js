import * as scran from "scran.js";
import * as utils from "./utils/general.js";
import * as putils from "./utils/pca.js";

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
        let keys = Object.keys(this.#pca_states);
        if (keys.length > 1) {
            return {
                "pcs": this.#cache.combined_buffer,
                "num_obs": this.#cache.num_cells,
                "num_pcs": this.#cache.total_dims
            };
        } else {
            return this.#pca_states[keys[0]].fetchPCs();
        }
    }

    /***************************
     ******** Compute **********
     ***************************/

    compute(weights) {
        this.changed = false;

        for (const x of Object.values(this.#pca_states)) {
            if (x.changed) {
                this.changed = true;
                break;
            }
        }

        if (this.changed) {
            let keys = Object.keys(this.#pca_states);
            
            if (keys.length > 1) {
                let used = [];
                let collected = [];
                let total = 0;
                let ncells = null;

                for (const x of keys) {
                    let state = this.#pca_states[x];
                    if (!state.valid()) {
                        continue;
                    }
                    
                    let curpcs = state.fetchPCs();
                    collected.push(curpcs.pcs);
                    used.push(x);

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
                    for (const x of used) {
                        if (!(x in weights)) {
                            throw new Error("no weight specified for '" + x + "'");
                        }
                        weight_arr.push(weights[x]);
                    }
                }

                let buffer = utils.allocateCachedArray(ncells * total, "Float64Array", this.#cache, "combined_buffer");
                scran.scaleByNeighbors(collected, ncells, { buffer: buffer, weights: weight_arr });
                this.#cache.num_cells = ncells;
                this.#cache.total_dims = total;
            }

            this.#parameters.weights = weights;
            this.changed = true;
        }

        return;
    }

    /***************************
     ******** Results **********
     ***************************/

    summary() {
        return {};
    }
}
