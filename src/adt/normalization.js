import * as scran from "scran.js"; 
import * as utils from "../utils/general.js";
import * as nutils from "../utils/normalization.js";
import * as qc_module from "./quality_control.js";
import * as filter_module from "../cell_filtering.js";

export const step_name = "adt_normalization";

/**
 * This step performs normalization and log-transformation on the QC-filtered ADT matrix from the {@linkplain AdtQualityControlState}.
 * It wraps the `groupedSizeFactors` and `logNormCounts` functions from [**scran.js**](https://github.com/jkanche/scran.js).
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class AdtNormalizationState extends nutils.NormalizationStateBase {
    #qc;
    #filter;
    #parameters;
    #cache;

    constructor(qc, filter, parameters = null, cache = null) {
        super();

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

    valid() {
        return this.#filter.hasAvailable(this.#parameters.target_matrix);
    }

    fetchNormalizedMatrix() {
        if (this.valid()) {
            if (!("matrix" in this.#cache)) {
                this.#raw_compute();
            }
            return this.#cache.matrix;
        } else {
            return null;
        }
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
     * @param {number} num_pcs - Number of PCs to use for creating a low-dimensional embedding for clustering.
     * @param {number} num_clusters - Number of clusters to create with k-means clustering.
     *
     * @return The object is updated with new results.
     */
    compute(num_pcs, num_clusters) {
        this.changed = false;

        if (this.#qc.changed || this.#filter.changed || num_pcs !== this.#parameters.num_pcs || num_clusters != this.#parameters.num_clusters) {
            if (this.valid()) {
                var mat = this.#filter.fetchFilteredMatrix({ type: this.#parameters.target_matrix });
                var total_buffer = utils.allocateCachedArray(mat.numberOfColumns(), "Float64Array", this.#cache, "total_buffer");
                nutils.subsetSums(this.#qc, this.#filter, total_buffer.array());

                var block = this.#filter.fetchFilteredBlock();
                var sf_buffer = utils.allocateCachedArray(mat.numberOfColumns(), "Float64Array", this.#cache, "sf_buffer");
                scran.quickAdtSizeFactors(mat, { totals: total_buffer, block: block, buffer: sf_buffer, numberOfPCs: num_pcs, numberOfClusters: num_clusters });
                this.changed = true;
            }
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
    let ghandle = handle.open(step_name);

    let parameters = AdtNormalizationState.defaults();
    let phandle = ghandle.open("parameters");
    parameters.num_pcs = phandle.open("num_pcs", { load: true }).values[0];
    parameters.num_clusters = phandle.open("num_clusters", { load: true }).values[0];

    let output;
    let cache = {};
    try {
        let rhandle = ghandle.open("results");
        
        if ("size_factors" in rhandle.children) {
            let sf = rhandle.open("size_factors", { load: true }).values;
            cache.sf_buffer = scran.createFloat64WasmArray(sf.length);
            cache.sf_buffer.set(sf);
        }

        output = new AdtNormalizationState(qc, filter, cache);
    } catch (e) {
        utils.freeCache(cache.sf_buffer);
        utils.freeCache(output);
        throw e;
    }

    return {
        state: output,
        parameters: {}
    }
}
