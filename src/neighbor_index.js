import * as scran from "scran.js";
import * as utils from "./utils/general.js";
import * as pca_module from "./pca.js";

/**
 * This step assembles the neighbor search indices from the PCs (see {@linkplain PcaState}) in preparation for nearest neighbor searches in downstream steps.
 * It wraps the `buildNeighborSearchIndex` function from [**scran.js**](https://github.com/jkanche/scran.js).
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class NeighborIndexState {
    #pca;
    #parameters;
    #cache;

    constructor(pca, parameters = null, cache = null) {
        if (!(pca instanceof pca_module.PcaState)) {
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

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     * Each argument is taken from the property of the same name in the `neighbor_index` property of the `parameters` of {@linkcode runAnalysis}.
     *
     * @param {boolean} approximate - Whether to create an approximate search index.
     * If `false`, an exact index is used.
     *
     * @return The object is updated with the new results.
     */
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

    /**
     * Obtain a summary of the state, typically for display on a UI like **kana**.
     *
     * @return An empty object.
     * This is just provided for consistency with the other classes.
     */
    summary() {
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
        state: new NeighborIndexState(pca, parameters, cache),
        parameters: {...parameters }
    };
}
