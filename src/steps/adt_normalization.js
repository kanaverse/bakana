import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as nutils from "./utils/normalization.js";
import * as qc_module from "./adt_quality_control.js";
import * as filter_module from "./cell_filtering.js";

export const step_name = "adt_normalization";

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
    }

    free() {
        utils.freeCache(this.#cache.matrix);
        utils.freeCache(this.#cache.total_buffer);
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
     * @return {?ScranMatrix} A ScranMatrix object containing the normalized ADT values,
     * available after running {@linkcode AdtNormalizationState#compute compute}.
     * Alternatively `null`, if no ADTs are present in the dataset.
     */
    fetchNormalizedMatrix() {
        if (!this.valid()) {
            return null;
        }

        if (!("matrix" in this.#cache)) {
            this.#raw_compute();
        }

        return this.#cache.matrix;
    }

    /**
     * @return {?Float64WasmArray} Array of length equal to the number of cells, 
     * containing the ADT-derived size factor for each cell in the (QC-filtered) dataset.
     * This is available after running {@linkcode AdtNormalizationState#compute compute}.
     * Alternatively `null`, if no ADTs are present in the dataset.
     */
    fetchSizeFactors() {
        if (!this.valid()) {
            return null;
        } else {
            return this.#cache.sf_buffer;
        }
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
     * @param {number} num_pcs - Number of PCs to use for creating a low-dimensional embedding for clustering.
     * @param {number} num_clusters - Number of clusters to create with k-means clustering.
     *
     * @return The object is updated with new results.
     */
    compute(num_pcs, num_clusters) {
        this.changed = false;

        if (this.#qc.changed || this.#filter.changed || num_pcs !== this.#parameters.num_pcs || num_clusters != this.#parameters.num_clusters) {
            if (this.valid()) {
                var mat = this.#filter.fetchFilteredMatrix().get("ADT");
                let total_buffer = nutils.subsetSums(this.#qc, this.#filter, mat, this.#cache, "total_buffer");

                var block = this.#filter.fetchFilteredBlock();
                var sf_buffer = utils.allocateCachedArray(mat.numberOfColumns(), "Float64Array", this.#cache, "sf_buffer");
                scran.quickAdtSizeFactors(mat, { 
                    totals: total_buffer, 
                    block: block, 
                    buffer: sf_buffer, 
                    numberOfPCs: num_pcs, 
                    numberOfClusters: num_clusters 
                });

                this.changed = true;
            }

            this.#parameters.num_pcs = num_pcs;
            this.#parameters.num_clusters = num_clusters;
        } 

        if (this.changed) {
            if (this.valid()) {
                this.#raw_compute();
            }
        }

        return;
    }

    static defaults() {
        return {
           num_pcs: 25,
           num_clusters: 20
        };
    }

    /*************************
     ******** Saving *********
     *************************/

    serialize(handle) {
        let ghandle = handle.createGroup(step_name);
        let phandle = ghandle.createGroup("parameters"); 
        phandle.writeDataSet("num_pcs", "Int32", [], this.#parameters.num_pcs);
        phandle.writeDataSet("num_clusters", "Int32", [], this.#parameters.num_clusters);

        let rhandle = ghandle.createGroup("results"); 
        if (this.valid()) {
            rhandle.writeDataSet("size_factors", "Float64", null, this.#cache.sf_buffer);
        }
    }
}

/**************************
 ******** Loading *********
 **************************/

export function unserialize(handle, qc, filter) {
    let cache = {};
    let parameters = AdtNormalizationState.defaults();
    let output;

    if (step_name in handle.children) {
        let ghandle = handle.open(step_name);

        let phandle = ghandle.open("parameters");
        parameters.num_pcs = phandle.open("num_pcs", { load: true }).values[0];
        parameters.num_clusters = phandle.open("num_clusters", { load: true }).values[0];

        try {
            let rhandle = ghandle.open("results");
            
            if ("size_factors" in rhandle.children) {
                let sf = rhandle.open("size_factors", { load: true }).values;
                cache.sf_buffer = scran.createFloat64WasmArray(sf.length);
                cache.sf_buffer.set(sf);
            }

            output = new AdtNormalizationState(qc, filter, parameters, cache);
        } catch (e) {
            utils.freeCache(cache.sf_buffer);
            utils.freeCache(output);
            throw e;
        }
    } else {
        // Fallback for v1.
        output = new AdtNormalizationState(qc, filter, parameters, cache);
    }

    return output;
}
