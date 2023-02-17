import * as scran from "scran.js";
import * as utils from "./utils/general.js";
import * as filter_module from "./cell_filtering.js";
import * as combine_module from "./combine_embeddings.js";

export const step_name = "batch_correction";

/**
 * Correct for batch effects in PC space based on mutual nearest neighbors.
 * This wraps the [`mnnCorrect`](https://jkanche.com/scran.js/global.html#mnnCorrect) function
 * from [**scran.js**](https://jkanche.com/scran.js).
 * 
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class BatchCorrectionState {
    #filter;
    #combined;
    #parameters;
    #cache;

    constructor(filter, combined, parameters = null, cache = null) {
        if (!(filter instanceof filter_module.CellFilteringState)) {
            throw new Error("'filter' should be a CellFilteringState object");
        }
        this.#filter = filter;

        if (!(combined instanceof combine_module.CombineEmbeddingsState)) {
            throw new Error("'pca' should be a CombineEmbeddingsState object");
        }
        this.#combined = combined;

        this.#parameters = (parameters === null ? {} : parameters);
        this.#cache = (cache === null ? {} : cache);
        this.changed = false;
    }

    free() {
        utils.freeCache(this.#cache.corrected);
    }

    /***************************
     ******** Getters **********
     ***************************/

    /**
     * @return {Float64WasmArray} Buffer containing the batch-corrected embeddings as a column-major dense matrix,
     * where the rows are the dimensions and the columns are the cells.
     * This is available after running {@linkcode BatchCorrectionState#compute compute}.
     */
    fetchCorrected() {
        return this.#cache.corrected;
    }

    /**
     * @return {number} Number of cells in {@linkcode BatchCorrectionState#fetchCorrected fetchCorrected}.
     */
    fetchNumberOfCells() {
        return this.#combined.fetchNumberOfCells();
    }

    /**
     * @return {number} Number of dimensions in {@linkcode BatchCorrectionState#fetchCorrected fetchCorrected}.
     */
    fetchNumberOfDimensions() {
        return this.#combined.fetchNumberOfDimensions();
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

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     *
     * @param {object} parameters - Parameter object, equivalent to the `batch_correction` property of the `parameters` of {@linkcode runAnalysis}.
     * @param {string} parameters.method - The correction method to use.
     * Currently this can be either `"mnn"` or `"none"`.
     * If `"mnn"`, it is recommended that upstream PCA steps (i.e., {@linkplain RnaPcaState} and {@linkplain AdtPcaState}) use `block_method = "weight"`.
     * @param {number} parameters.num_neighbors - Number of neighbors to use during MNN correction.
     * @param {boolean} parameters.approximate - Whether to use an approximate method to identify MNNs.
     *
     * @return The object is updated with new results.
     */
    compute(parameters) {
        let { method, num_neighbors, approximate} = parameters;
        this.changed = false;

        if (this.#filter.changed || this.#combined.changed) {
            this.changed = true;
        }
        let block = this.#filter.fetchFilteredBlock();
        let needs_correction = (method == "mnn" && block !== null);

        if (this.changed || method !== this.#parameters.method || num_neighbors !== this.#parameters.num_neighbors || approximate !== this.#parameters.approximate) { 
            if (needs_correction) {
                let pcs = this.#combined.fetchCombined();
                let corrected = utils.allocateCachedArray(pcs.length, "Float64Array", this.#cache, "corrected");
                scran.mnnCorrect(pcs, block, { 
                    k: num_neighbors, 
                    buffer: corrected, 
                    numberOfCells: this.#combined.fetchNumberOfCells(), 
                    numberOfDims: this.#combined.fetchNumberOfDimensions(), 
                    approximate: approximate 
                });
                this.changed = true;
            }
        }

        if (this.changed) {
            // If no correction is actually required, we shouldn't respond to
            // changes in parameters, because they won't have any effect.
            if (!needs_correction) {
                utils.freeCache(this.#cache.corrected);
                this.#cache.corrected = this.#combined.fetchCombined().view();
            }
        }

        // Updating all parameters, even if they weren't used.
        this.#parameters.method = method;
        this.#parameters.num_neighbors = num_neighbors;
        this.#parameters.approximate = approximate;
        return;
    }

    static defaults() {
        return {
            method: "mnn",
            num_neighbors: 15,
            approximate: true
        };
    }
}

/**************************
 ******** Loading *********
 **************************/

export function unserialize(handle, filter, combined) {
    let cache = {};
    let parameters = BatchCorrectionState.defaults();
    let output;
    
    if (step_name in handle.children) {
        let ghandle = handle.open(step_name);

        let phandle = ghandle.open("parameters"); 
        parameters.method = phandle.open("method", { load: true }).values[0];
        parameters.num_neighbors = phandle.open("num_neighbors", { load: true }).values[0];
        parameters.approximate = phandle.open("approximate", { load: true }).values[0] > 0;

        try {
            let rhandle = ghandle.open("results");

            if ("corrected" in rhandle.children) {
                let corrected = rhandle.open("corrected", { load: true }).values;
                cache.corrected = scran.createFloat64WasmArray(corrected.length);
                cache.corrected.set(corrected);
            } else {
                // Creating a view from the upstream combined state.
                let pcs = combined.fetchCombined();
                cache.corrected = pcs.view();
            }

            output = new BatchCorrectionState(filter, combined, parameters, cache);
        } catch (e) {
            utils.freeCache(cache.corrected);
            utils.freeCache(output);
            throw e;
        }
    } else {
        // Fallback for v1.
        let ghandle = handle.open("pca");

        let rhandle = ghandle.open("results");
        if ("corrected" in rhandle.children) {
            let corrected = rhandle.open("corrected", { load: true }).values;
            let corbuffer = utils.allocateCachedArray(corrected.length, "Float64Array", cache, "corrected");
            corbuffer.set(corrected);
        } else {
            cache.corrected = combined.fetchCombined().view();
        }

        output = new BatchCorrectionState(filter, combined, parameters, cache);
    }

    return output;
}
