import * as scran from "scran.js";
import * as utils from "./utils/general.js";
import * as pca_module from "./pca.js";

/**
 * This step assembles the neighbor search indices from the PCs (see {@linkcode pca}) in preparation for nearest neighbor searches in downstream steps.
 * It wraps the `buildNeighborSearchIndex` function from [**scran.js**](https://github.com/jkanche/scran.js).
 *
 * The parameters in {@linkcode runAnalysis} should be an object containing:
 *
 * - `approximate`: boolean indicating whether to create an approximate search index.
 *
 * Calling the **`results()`** method for the relevant state instance will return an empty object.
 * 
 * Methods not documented here are not part of the stable API and should not be used by applications.
 *
 * @namespace neighbor_index
 */

export class State {
    #pca;
    #parameters;
    #cache;

    constructor(pca, parameters = null, cache = null) {
        if (!(pca instanceof pca_module.State)) {
            throw new Error("'pca' should be a State object from './pca.js'");
        }
        this.#pca = pca;

        this.#parameters = (parameters === null ? {} : parameters);
        this.#cache = (cache === null ? {} : cache);
        this.changed = false;
    }

    free() {
        utils.freeCache(this.#cache.raw);
    }

    /***************************
     ******** Getters **********
     ***************************/

    fetchIndex() {
        if (!("raw" in this.#cache)) {
            this.#raw_compute(this.#parameters.approximate);
        }
        return this.#cache.raw;
    }

    /***************************
     ******** Compute **********
     ***************************/

    #raw_compute(approximate) {
        var pcs = this.#pca.fetchPCs();
        this.#cache.raw = scran.buildNeighborSearchIndex(pcs.pcs, { approximate: approximate, numberOfDims: pcs.num_pcs, numberOfCells: pcs.num_obs });
        return;
    }

    compute(approximate) {
        this.changed = false;

        if (this.#pca.changed || approximate != this.#parameters.approximate) {
            utils.freeCache(this.#cache.raw);
            this.#raw_compute(approximate);
            this.#parameters.approximate = approximate;
            this.changed = true;
        }

        return;
    }

    /***************************
     ******** Results **********
     ***************************/

    results() {
        return {};
    }

    /*************************
     ******** Saving *********
     *************************/

    serialize(handle) {
        let ghandle = handle.createGroup("neighbor_index");

        {
            let phandle = ghandle.createGroup("parameters");
            phandle.writeDataSet("approximate", "Uint8", [], Number(this.#parameters.approximate));
        }

        ghandle.createGroup("results");
        return;
    }
}

/**************************
 ******** Loading *********
 **************************/

export function unserialize(handle, pca) {
    let ghandle = handle.open("neighbor_index");

    let parameters = {};
    {
        let phandle = ghandle.open("parameters");
        parameters = {
            approximate: phandle.open("approximate", { load: true }).value > 0
        };
    }

    let cache = {};

    return { 
        state: new State(pca, parameters, cache),
        parameters: {...parameters }
    };
}