import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as nutils from "./utils/normalization.js";
import * as qc_module from "./rna_quality_control.js";
import * as filter_module from "./cell_filtering.js";

export const step_name = "rna_normalization";

/**
 * This step performs normalization and log-transformation on the QC-filtered matrix from the {@linkplain QualityControlState}.
 * It wraps the [`logNormCounts`](https://www.jkanche.com/scran.js/global.html#logNormCounts) function
 * from [**scran.js**](https://github.com/jkanche/scran.js).
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class RnaNormalizationState {
    #qc
    #filter;
    #parameters;
    #cache;

    constructor(qc, filter, parameters = null, cache = null) {
        if (!(qc instanceof qc_module.RnaQualityControlState)) {
            throw new Error("'qc' should be a RnaQualityControlState object");
        }
        this.#qc = qc;

        if (!(filter instanceof filter_module.CellFilteringState)) {
            throw new Error("'filter' should be a CellFilteringState object");
        }
        this.#filter = filter;

        this.#parameters = (parameters === null ? {} : parameters);
        this.#cache = (cache === null ? {} : cache);
        this.changed = false;
    }

    free() {
        utils.freeCache(this.#cache.matrix);
        utils.freeCache(this.#cache.sum_buffer);
    }

    /***************************
     ******** Getters **********
     ***************************/

    valid() {
        let filtered = this.#filter.fetchFilteredMatrix();
        return filtered.has("RNA");
    }

    /**
     * @return {external:ScranMatrix} A {@linkplain external:ScranMatrix ScranMatrix} object containing normalized expression values,
     * available after running {@linkcode RnaNormalizationState#compute compute}.
     */
    fetchNormalizedMatrix() {
        if (!("matrix" in this.#cache)) {
            this.#raw_compute();
        }
        return this.#cache.matrix;
    }

    /**
     * @return {Float64WasmArray} Array of length equal to the number of cells, 
     * containing the RNA-derived size factor for each cell.
     * This is available after running {@linkcode RnaNormalizationState#compute compute}.
     */
    fetchSizeFactors() {
        let buff;
        if (this.#cache.sum_buffer) {
            buff = utils.allocateCachedArray(this.#cache.sum_buffer.length, "Float64Array", this.#cache, "centered_buffer");
            scran.centerSizeFactors(this.#cache.sum_buffer, { buffer: buff, block: this.#filter.fetchFilteredBlock() })
        }
        return buff;
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

    #raw_compute() {
        var mat = this.#filter.fetchFilteredMatrix().get("RNA");
        let buffer = nutils.subsetSums(this.#qc, this.#filter, mat, this.#cache, "sum_buffer");

        var block = this.#filter.fetchFilteredBlock();
        utils.freeCache(this.#cache.matrix);
        this.#cache.matrix = scran.logNormCounts(mat, { sizeFactors: buffer, block: block, allowZeros: true });
        return;
    }

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     *
     * @param {object} parameters - Parameter object, equivalent to the `rna_normalization` property of the `parameters` of {@linkcode runAnalysis}.
     *
     * @return The object is updated with new results.
     */
    compute(parameters) {
        this.changed = false;
        if (this.#qc.changed || this.#filter.changed) {
            if (this.valid()) {
                this.changed = true;
            }
        } 

        if (this.changed) {
            this.#raw_compute();
        }
        return;
    }

    static defaults() {
        return {};
    }
}

/**************************
 ******** Loading *********
 **************************/

export function unserialize(handle, qc, filter) {
    return new RnaNormalizationState(qc, filter);
}
