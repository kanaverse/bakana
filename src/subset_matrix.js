import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";

/**
 * This step subsets datasets based on a list of indices to continue with the rest of the analysis
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class SubsetCellsState {
    #qc
    #parameters;
    #cache;

    constructor(qc, parameters = null, cache = null) {
        if (!(qc instanceof qc_module.QualityControlState)) {
            throw new Error("'qc' should be a State object from './quality_control.js'");
        }
        this.#qc = qc;

        this.#parameters = (parameters === null ? {} : parameters);
        this.#cache = (cache === null ? {} : cache);
        this.changed = false;
    }

    free() {
        utils.freeCache(this.#cache.matrix);
    }

    /***************************
     ******** Getters **********
     ***************************/

    fetchSubsetMatrix() {
        if (!("matrix" in this.#cache)) {
            this.#raw_compute(this.#parameters.indices);
        }
        return this.#cache.matrix;
    }

    fetchIndices() {
        return this.#parameters.indices ? this.#parameters.indices : [];
    }


    /***************************
     ******** Compute **********
     ***************************/

     #raw_compute(indices) {
        var mat = this.#qc.fetchFilteredMatrix();

        if (!indices || (Array.isArray(indices) && indices.length == 0)) {
            this.#cache.matrix = mat;
        } else {
            this.#cache.matrix = scran.subsetColumns(mat, indices);
        }

        return;
    }

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     *
     * @return The object is updated with new results.
     */
     compute(indices) {
        this.changed = false;
        if (this.#qc.changed) {
            this.changed = true;
        } 

        if (this.changed) {
            this.#raw_compute(indices);
            this.#parameters.indices = indices;
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
        // Token effort.
        let ghandle = handle.createGroup("subset");
        ghandle.createGroup("parameters"); 
        ghandle.createGroup("results"); 
    }
}

/**************************
 ******** Loading *********
 **************************/

export function unserialize(handle, qc) {
    return {
        state: new SubsetCellsState(qc),
        parameters: {}
    }
}
