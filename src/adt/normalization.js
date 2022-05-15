import * as scran from "scran.js"; 
import * as utils from "../utils/general.js";
import * as nutils from "../utils/normalization.js";
import * as qc_module from "./quality_control.js";
import * as filter_module from "../cell_filtering.js";

/**
 * This step performs normalization and log-transformation on the QC-filtered ADT matrix from the {@linkplain AdtQualityControlState}.
 * It wraps the `groupedSizeFactors` and `logNormCounts` functions from [**scran.js**](https://github.com/jkanche/scran.js).
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class AdtNormalizationState {
    #qc;
    #filter;
    #parameters;
    #cache;

    constructor(qc, filter, parameters = null, cache = null) {
        if (!(qc instanceof qc_module.AdtQualityControlState)) {
            throw new Error("'qc' should be a AdtQualityControlState object");
        }
        this.#qc = qc;

        if (!(filter instanceof filter_module.CellFilteringState)) {
            throw new Error("'filter' should be a CellFilteringState object");
        }
        this.#filter = filter;

        this.#parameters = (parameters === null ? {} : parameters);
        this.#cache = (cache === null ? {} : cache);
        this.changed = false;

        this.#parameters.target_matrix = "ADT";
    }

    free() {
        utils.freeCache(this.#cache.matrix);
        utils.freeCache(this.#cache.exp_buffer);
        utils.freeCache(this.#cache.total_buffer);
        utils.freeCache(this.#cache.sf_buffer);
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
        var buffer = utils.allocateCachedArray(mat.numberOfColumns(), "Float64Array", this.#cache, "exp_buffer"); // re-using the buffer.
        mat.row(index, { buffer: buffer });
        return buffer.slice();
    }

    /***************************
     ******** Compute **********
     ***************************/

    useRNAMatrix() {
        // For testing only!
        this.#parameters.target_matrix = "RNA";
        return;
    }

    #raw_compute() {
        var mat = this.#filter.fetchFilteredMatrix({ type: this.#parameters.target_matrix });
        var block = this.#filter.fetchFilteredBlock();

        var buffer = this.#cache.sf_buffer;
        if (buffer.length != mat.numberOfColumns()) {
            throw new Error("length of size factor vector should equal number of columns after QC");
        }

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

        if (this.#qc.changed || this.#filter.changed) {
            var mat = this.#filter.fetchFilteredMatrix({ type: this.#parameters.target_matrix });
            var total_buffer = utils.allocateCachedArray(mat.numberOfColumns(), "Float64Array", this.#cache, "total_buffer");
            nutils.subsetSums(this.#qc, this.#filter, total_buffer.array());

            var block = this.#filter.fetchFilteredBlock();
            var sf_buffer = utils.allocateCachedArray(mat.numberOfColumns(), "Float64Array", this.#cache, "sf_buffer");
            scran.quickAdtSizeFactors(mat, { totals: total_buffer, block: block, buffer: sf_buffer });

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
        let ghandle = handle.createGroup("adt-normalization");
        ghandle.createGroup("parameters"); // Token effort.
        let rhandle = ghandle.createGroup("results"); 
        rhande.writeDataSet("size_factors", "Float64", null, this.#cache.sf_buffer);
    }
}

/**************************
 ******** Loading *********
 **************************/

export function unserialize(handle, qc) {
    let ghandle = handle.open("adt-normalization");
    let rhandle = ghandle.open("results");
    let sf = rhandle.open("size_factors", { load: true }).values;

    let cache = {};
    cache.sf_buffer = scran.createFloat64WasmArray(sf.length);
    cache.sf_buffer.set(sf);

    return {
        state: new AdtNormalizationState(qc, null, cache),
        parameters: {}
    }
}
