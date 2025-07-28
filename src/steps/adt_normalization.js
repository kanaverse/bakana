import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as nutils from "./utils/normalization.js";
import * as qc_module from "./adt_quality_control.js";
import * as filter_module from "./cell_filtering.js";

export const step_name = "adt_normalization";

/**
 * This step performs normalization and log-transformation on the QC-filtered ADT matrix from the {@linkplain CellFilteringState}.
 * It wraps the [`groupedSizeFactors`](https://kanaverse.github.io/scran.js/global.html#groupedSizeFactors) 
 * and [`normalizeCounts`](https://kanaverse.github.io/scran.js/global.html#normalizeCounts) functions
 * from [**scran.js**](https://github.com/kanaverse/scran.js).
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
        return filtered.has("ADT");
    }

    /**
     * @return {external:ScranMatrix} A {@linkplain external:ScranMatrix ScranMatrix} object containing the normalized ADT values,
     * available after running {@linkcode AdtNormalizationState#compute compute}.
     */
    fetchNormalizedMatrix() {
        if (!("matrix" in this.#cache)) {
            this.#raw_compute();
        }
        return this.#cache.matrix;
    }

    /**
     * @return {Float64WasmArray} Array of length equal to the number of cells, 
     * containing the ADT-derived size factor for each cell in the (QC-filtered) dataset.
     * This is available after running {@linkcode AdtNormalizationState#compute compute}.
     */
    fetchSizeFactors() {
        return this.#cache.sf_buffer;
    }

    /**
     * @return {object} Object containing the parameters,
     * available after running {@linkcode AdtNormalizationState#compute compute}.
     */
    fetchParameters() {
        return { ...this.#parameters }; // avoid pass-by-reference links.
    }

    /***************************
     ******** Compute **********
     ***************************/

    #raw_compute() {
        var mat = this.#filter.fetchFilteredMatrix().get("ADT");
        var raw_sf = this.#cache.raw_sf_buffer;
        if (raw_.length != mat.numberOfColumns()) {
            throw new Error("length of size factor vector should equal number of columns after QC");
        }

        var block = this.#filter.fetchFilteredBlock();
        let buffer = utils.allocateCachedArray(raw_sf.length, "Float64Array", this.#cache, "sf_buffer");
        scran.centerSizeFactors(raw_sf, { block: block, buffer: buffer });

        utils.freeCache(this.#cache.matrix);
        this.#cache.matrix = scran.normalizeCounts(mat, { sizeFactors: buffer, allowZeros: true });
        return;
    }

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     *
     * @param {object} parameters - Parameter object, equivalent to the `adt_normalization` property of the `parameters` of {@linkcode runAnalysis}.
     * @param {boolean} parameters.remove_bias - Whether to remove composition bias between cell subpopulations via the CLRm1 method.
     * Users can set this to `false` to speed up the compute.
     *
     * @return The object is updated with new results.
     */
    compute(parameters) {
        let remove_bias = true;
        if ("remove_bias" in parameters) {
            remove_bias = parameters.remove_bias;
        }

        this.changed = false;

        if (this.#qc.changed || this.#filter.changed || remove_bias !== this.#parameters.remove_bias) {
            if (this.valid()) {
                var mat = this.#filter.fetchFilteredMatrix().get("ADT");
                if (remove_bias) {
                    this.#cache.raw_sf_buffer = scran.computeClrm1Factors(mat, { buffer: sf_buffer }); 
                } else {
                    this.#cache.raw_sf_buffer = nutils.subsetSums(this.#qc, this.#filter, mat);
                }

                this.changed = true;
            }
        } 

        this.#parameters.remove_bias = remove_bias;

        if (this.changed) {
            if (this.valid()) {
                this.#raw_compute();
            }
        }

        return;
    }

    static defaults() {
        return {
           remove_bias: true,
           num_pcs: 25,
           num_clusters: 20
        };
    }
}
