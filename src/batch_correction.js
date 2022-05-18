import * as scran from "scran.js";
import * as utils from "./utils/general.js";
import * as filter_module from "./cell_filtering.js";
import * as combine_module from "./combine_embeddings.js";

export const step_name = "batch_correction";

export class BatchCorrectionState {
    #filter;
    #combined;
    #parameters;
    #cache;

    constructor(filter, combined, parameters = null, cache = null) {
        if (!(filter instanceof filter_module.CellFilteringState)) {
            throw new Error("'filter' should be a CellFilteringState object");
        }
        this.#filter = filter;

        if (!(combined instanceof combine_module.CombineEmbeddingsState)) {
            throw new Error("'pca' should be a CombineEmbeddingsState object");
        }
        this.#combined = combined;

        this.#parameters = (parameters === null ? {} : parameters);
        this.#cache = (cache === null ? {} : cache);
        this.changed = false;
    }

    free() {
        utils.freeCache(this.#cache.corrected);
    }

    /***************************
     ******** Getters **********
     ***************************/

    fetchPCs() {
        let upstream = this.#combined.fetchPCs();
        upstream.pcs = this.#cache.corrected;
        return upstream;
    }

    /***************************
     ******** Compute **********
     ***************************/

    compute(method, num_neighbors, approximate) {
        this.changed = false;

        if (this.#filter.changed || this.#combined.changed) {
            this.changed = true;
        }
        let block = this.#filter.fetchFilteredBlock();
        let needs_correction = (method == "mnn" && block !== null);

        if (this.changed || method !== this.#parameters.method || num_neighbors !== this.#parameters.num_neighbors || approximate !== this.#parameters.approximate) { 
            if (needs_correction) {
                let corrected = utils.allocateCachedArray(pcs.length, "Float64Array", this.#cache, "corrected");
                let pcs = this.#combined.fetchPCs();
                scran.mnnCorrect(pcs.pcs, block, { k: num_neighbors, buffer: corrected, numberOfCells: pcs.num_obs, numberOfDims: pcs.num_pcs, approximate: approximate });
                this.changed = true;
            }
        }

        if (this.changed) {
            // If no correction is actually required, we shouldn't respond to
            // changes in parameters, because they won't have any effect.
            if (!needs_correction) {
                utils.freeCache(this.#cache.corrected);
                let upstream = this.#combined.fetchPCs();
                this.#cache.corrected = upstream.pcs.view();
            }
        }

        // Updating all parameters, even if they weren't used.
        this.#parameters.method = method;
        this.#parameters.num_neighbors = num_neighbors;
        this.#parameters.approximate = approximate;
        return;
    }

    static defaults() {
        return {
            method: "mnn",
            num_neighbors: 15,
            approximate: true
        };
    }

    /**************************
     ******** Results**********
     **************************/

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
            phandle.writeDataSet("method", "String", [], this.#parameters.method);
            phandle.writeDataSet("num_neighbors", "Int32", [], this.#parameters.num_neighbors);
            phandle.writeDataSet("approximate", "Uint8", [], Number(this.#parameters.approximate));
        }

        {
            let rhandle = ghandle.createGroup("results");
            let pcs = this.fetchPCs();
            if (pcs.pcs.owner !== null) {
                // If it's not a view, we save it; otherwise we assume
                // that we can recover it from the upstream state.
                rhandle.writeDataSet("corrected", "Float64", [pcs.num_obs, pcs.num_pcs], pcs.pcs); // remember, it's transposed.
            }
        }
    }
}

/**************************
 ******** Loading *********
 **************************/

export function unserialize(handle, filter, combined) {
    let ghandle = handle.open(step_name);

    let parameters = BatchCorrectionState.defaults();
    {
        let phandle = ghandle.open("parameters"); 
        parameters.method = phandle.open("method", { load: true }).values[0];
        parameters.num_neighbors = phandle.open("num_neighbors", { load: true }).values[0];
        parameters.approximate = phandle.open("approximate", { load: true }).values[0] > 0;
    }

    let output;
    let cache = {};
    try {
        let rhandle = ghandle.open("results");

        if ("corrected" in rhandle.children) {
            let corrected = rhandle.open("corrected", { load: true }).values;
            cache.corrected = scran.createFloat64WasmArray(corrected.length);
            cache.corrected.set(corrected);
        } else {
            // Creating a view from the upstream combined state.
            let pcs = combined.fetchPCs();
            cache.corrected = pcs.pcs.view();
        }

        output = new BatchCorrectionState(filter, combined, parameters, cache);
    } catch (e) {
        utils.freeCache(cache.corrected);
        utils.freeCache(output);
        throw e;
    }

    return {
        state: output,
        parameters: parameters
    };
}
