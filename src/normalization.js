import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as qc_module from "./quality_control.js";

/**
 * This step performs normalization and log-transformation on the QC-filtered matrix from the {@linkplain QualityControlState}.
 * It wraps the `logNormCounts` function from [**scran.js**](https://github.com/jkanche/scran.js).
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class NormalizationState {
    #qc;
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

    fetchNormalizedMatrix() {
        if (!("matrix" in this.#cache)) {
            this.#raw_compute();
        }
        return this.#cache.matrix;
    }

    /**
     * Extract normalized expression values.
     * @param {number} index - An integer specifying the row index to extract.
     * @return A Float64Array of length equal to the number of (QC-filtered) cells, containing the log-normalized expression values for each cell.
     */
    fetchExpression(index) {
        var mat = this.fetchNormalizedMatrix();
        var buffer = utils.allocateCachedArray(mat.numberOfColumns(), "Float64Array", this.#cache); // re-using the buffer.
        mat.row(index, { buffer: buffer });
        return buffer.slice();
    }

    /***************************
     ******** Compute **********
     ***************************/

    #raw_compute() {
        var mat = this.#qc.fetchFilteredMatrix();
        var buffer = utils.allocateCachedArray(mat.numberOfColumns(), "Float64Array", this.#cache);

        var discards = this.#qc.fetchDiscards();
        var sums = this.#qc.fetchSums({ unsafe: true }); // Better not have any more allocations in between now and filling of size_factors!

        // Reusing the totals computed earlier.
        var size_factors = buffer.array();
        var j = 0;
        discards.array().forEach((x, i) => {
            if (!x) {
                size_factors[j] = sums[i];
                j++;
            }
        });

        if (j != mat.numberOfColumns()) {
            throw "normalization and filtering are not in sync";
        }

        var block = this.#qc.fetchFilteredBlock();

        utils.freeCache(this.#cache.matrix);
        this.#cache.matrix = scran.logNormCounts(mat, { sizeFactors: buffer, block: block });
        return;
    }

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     *
     * @return The object is updated with new results.
     */
    compute() {
        this.changed = false;
        if (this.#qc.changed) {
            this.changed = true;
        } 

        if (this.changed) {
            this.#raw_compute();
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
        let ghandle = handle.createGroup("normalization");
        ghandle.createGroup("parameters"); 
        ghandle.createGroup("results"); 
    }
}

/**************************
 ******** Loading *********
 **************************/

export function unserialize(handle, qc) {
    return {
        state: new NormalizationState(qc),
        parameters: {}
    }
}
