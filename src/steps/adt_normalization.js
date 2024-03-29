import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as nutils from "./utils/normalization.js";
import * as qc_module from "./adt_quality_control.js";
import * as filter_module from "./cell_filtering.js";

export const step_name = "adt_normalization";

/**
 * This step performs normalization and log-transformation on the QC-filtered ADT matrix from the {@linkplain CellFilteringState}.
 * It wraps the [`groupedSizeFactors`](https://kanaverse.github.io/scran.js/global.html#groupedSizeFactors) 
 * and [`logNormCounts`](https://kanaverse.github.io/scran.js/global.html#logNormCounts) functions
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
        utils.freeCache(this.#cache.total_buffer);
        utils.freeCache(this.#cache.sf_buffer);
        utils.freeCache(this.#cache.centered_sf_buffer);
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
        let buff;
        if (this.#cache.sf_buffer) {
            buff = utils.allocateCachedArray(this.#cache.sf_buffer.length, "Float64Array", this.#cache, "centered_sf_buffer");
            scran.centerSizeFactors(this.#cache.sf_buffer, { buffer: buff, block: this.#filter.fetchFilteredBlock() })
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
        var mat = this.#filter.fetchFilteredMatrix().get("ADT");
        var block = this.#filter.fetchFilteredBlock();

        var buffer = this.#cache.sf_buffer;
        if (buffer.length != mat.numberOfColumns()) {
            throw new Error("length of size factor vector should equal number of columns after QC");
        }

        utils.freeCache(this.#cache.matrix);
        this.#cache.matrix = scran.logNormCounts(mat, { sizeFactors: buffer, block: block, allowZeros: true });
        return;
    }

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     *
     * @param {object} parameters - Parameter object, equivalent to the `adt_normalization` property of the `parameters` of {@linkcode runAnalysis}.
     * @param {boolean} parameters.remove_bias - Whether to remove composition bias between cell subpopulations.
     * This is done by clustering cells and computing median-based size factors between the average pseudo-cells for each cluster.
     * Users can set this to `false` to speed up the compute.
     * @param {number} parameters.num_pcs - Number of PCs to use for creating a low-dimensional embedding for clustering.
     * Only used if `remove_bias = true`.
     * @param {number} parameters.num_clusters - Number of clusters to create with k-means clustering.
     * Only used if `remove_bias = true`.
     *
     * @return The object is updated with new results.
     */
    compute(parameters) {
        const { num_pcs, num_clusters } = parameters;
        let remove_bias = true;
        if ("remove_bias" in parameters) {
            remove_bias = parameters.remove_bias;
        }

        this.changed = false;

        if (this.#qc.changed || 
            this.#filter.changed || 
            remove_bias !== this.#parameters.remove_bias ||
            (
                remove_bias &&
                (
                    num_pcs !== this.#parameters.num_pcs || 
                    num_clusters != this.#parameters.num_clusters
                ) 
            )
        ) {
            if (this.valid()) {
                var mat = this.#filter.fetchFilteredMatrix().get("ADT");
                let total_buffer = nutils.subsetSums(this.#qc, this.#filter, mat, this.#cache, "total_buffer");
                var block = this.#filter.fetchFilteredBlock();
                var sf_buffer = utils.allocateCachedArray(mat.numberOfColumns(), "Float64Array", this.#cache, "sf_buffer");

                if (remove_bias) {
                    scran.quickAdtSizeFactors(mat, { 
                        totals: total_buffer, 
                        block: block, 
                        buffer: sf_buffer, 
                        numberOfPCs: num_pcs, 
                        numberOfClusters: num_clusters 
                    });
                } else {
                    scran.centerSizeFactors(total_buffer, { buffer: sf_buffer, block: block });
                }

                this.changed = true;
            }

        } 

        this.#parameters.remove_bias = remove_bias;
        this.#parameters.num_pcs = num_pcs;
        this.#parameters.num_clusters = num_clusters;

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
