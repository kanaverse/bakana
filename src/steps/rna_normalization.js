import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as nutils from "./utils/normalization.js";
import * as qc_module from "./rna_quality_control.js";
import * as filter_module from "./cell_filtering.js";

export const step_name = "rna_normalization";

/**
 * This step performs normalization and log-transformation on the QC-filtered matrix from the {@linkplain QualityControlState}.
 * It wraps the [`normalizeCounts`](https://kanaverse.github.io/scran.js/global.html#normalizeCounts) function
 * from [**scran.js**](https://github.com/kanaverse/scran.js).
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
        utils.freeCache(this.#cache.sf_buffer);
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
        return this.#cache.sf_buffer;
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
        let raw_sf = nutils.subsetSums(this.#qc, this.#filter, mat);

        let block = this.#filter.fetchFilteredBlock();
        let buffer = utils.allocateCachedArray(raw_sf.length, "Float64Array", this.#cache, "sf_buffer");
        scran.centerSizeFactors(raw_sf, { block: block, buffer: buffer });

        utils.freeCache(this.#cache.matrix);
        this.#cache.matrix = scran.normalizeCounts(mat, { sizeFactors: buffer, allowZeros: true });
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
