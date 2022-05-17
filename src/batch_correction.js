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

    /***************************
     ******** Getters **********
     ***************************/

    fetchPCs() {
        let upstream = this.#combined.fetchPCs();
        if (this.#parameters.method == "mnn") {
            upstream.pcs = this.#cache.corrected;
        } 
        return upstream;
    }

    /***************************
     ******** Compute **********
     ***************************/

    compute(method, num_neighbors, approximate) {
        this.changed = false;

        if (this.#filter.changed || this.#combined.changed || 
            method !== this.#parameters.method || 
            num_neighbors !== this.#parameters.num_neighbors || 
            approximate !== this.#parameters.approximate) 
        { 
            let block = this.#filter.fetchFilteredBlock();
            if (method == "mnn" && block !== null) {
                let corrected = utils.allocateCachedArray(pcs.length, "Float64Array", this.#cache, "corrected");
                let pcs = this.#combined.fetchPCs();
                scran.mnnCorrect(pcs.pcs, block, { k: num_neighbors, buffer: corrected, numberOfCells: pcs.num_obs, numberOfDims: pcs.num_pcs, approximate: approximate });
            } else {
                utils.freeCache(this.#cache.corrected);
                delete this.#cache.corrected;
            }

            this.#parameters.method = method;
            this.#parameters.num_neighbors = num_neighbors;
            this.#parameters.approximate = approximate;
            this.changed = true;
        }
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
            rhandle.writeDataSet("corrected", "Float64", [pcs.num_obs, pcs.num_pcs], pcs.pcs); // remember, it's transposed.
        }
    }
}

/**************************
 ******** Loading *********
 **************************/

export function unserialize(handle, filter, combined) {
    let ghandle = handle.open(step_name);

    let parameters = BatchCorrection.defaults();
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
        let corrected = rhandle.open("corrected", { load: true }).values;
        cache.corrected = scran.createFloat64WasmArray(corrected.length);
        cache.corrected.set(corrected);
        output = new BatchCorrectionState(filter, combined, parameters, cache);
    } catch (e) {
        utils.freeCache(cache.pcs);
        utils.freeCache(output);
        throw e;
    }

    return {
        state: output,
        parameters: parameters
    };
}
