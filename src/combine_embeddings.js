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
        return {
            "pcs": this.#cache.combined_buffer,
            "num_obs": this.#cache.num_cells,
            "num_pcs": this.#cache.total_dims
        };
    }

    /***************************
     ******** Compute **********
     ***************************/

    compute(weights) {
        this.changed = false;

        for (const x of Object.values(this.#pca_states)) {
            if (x.valid() && x.changed) {
                this.changed = true;
                break;
            }
        }

        if (this.changed) {
            let to_use = [];
            for (const [k, v] of Object.entries(this.#pca_states)) {
                if (v.valid()) {
                    to_use.push(k);
                }
            }

            // We merge the embeddings if there's more than one. If there's
            // just one, we create a view on it directly and avoid allocating
            // an unnecessary buffer.
            if (to_use.length > 1) {
                let collected = [];
                let total = 0;
                let ncells = null;

                for (const x of to_use) {
                    let curpcs = this.#pca_states[x].fetchPCs();
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
            } else {
                let pcs = this.#pca_states[to_use[0]].fetchPCs();
                utils.freeCache(this.#cache.combined_buffer);
                this.#cache.combined_buffer = pcs.pcs.view();
                this.#cache.num_cells = pcs.num_obs;
                this.#cache.total_dims = pcs.num_pcs;
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
