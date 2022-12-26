import * as scran from "scran.js";
import * as utils from "./utils/general.js";
import * as correct_module from "./batch_correction.js";

export const step_name = "neighbor_index";

/**
 * This step assembles the neighbor search indices from the PCs (see {@linkplain PcaState}) in preparation for nearest neighbor searches in downstream steps.
 * It wraps the [`buildNeighborSearchIndex`](https://jkanche.com/scran.js/global.html#buildNeighborSearchIndex) function 
 * from [**scran.js**](https://github.com/jkanche/scran.js).
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class NeighborIndexState {
    #correct;
    #parameters;
    #cache;

    constructor(correct, parameters = null, cache = null) {
        if (!(correct instanceof correct_module.BatchCorrectionState)) {
            throw new Error("'correct' should be a BatchCorrectionState object");
        }
        this.#correct = correct;

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

    /**
     * @return {BuildNeighborSearchIndexResults} Index for a nearest-neighbor search,
     * available after running {@linkcode NeighborIndexState#compute compute}.
     */
    fetchIndex() {
        if (!("raw" in this.#cache)) {
            this.#raw_compute(this.#parameters.approximate);
        }
        return this.#cache.raw;
    }

    /**
     * @return {object} Object containing the parameters.
     */
    fetchParameters() {
        return { ...this.#parameters }; // avoid pass-by-reference links.
    }

    /***************************
     ******** Compute **********
     ***************************/

    static defaults() {
        return {
            approximate: true
        };
    }

    #raw_compute(approximate) {
        this.#cache.raw = scran.buildNeighborSearchIndex(this.#correct.fetchCorrected(), { 
            approximate: approximate, 
            numberOfDims: this.#correct.fetchNumberOfDimensions(),
            numberOfCells: this.#correct.fetchNumberOfCells()
        });
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

        if (this.#correct.changed || approximate != this.#parameters.approximate) {
            utils.freeCache(this.#cache.raw);
            this.#raw_compute(approximate);
            this.#parameters.approximate = approximate;
            this.changed = true;
        }

        return;
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
            approximate: phandle.open("approximate", { load: true }).values[0] > 0
        };
    }

    let cache = {};
    return new NeighborIndexState(pca, parameters, cache);
}
