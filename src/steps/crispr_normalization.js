import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as filter_module from "./cell_filtering.js";

export const step_name = "crispr_normalization";

/**
 * This step performs normalization and log-transformation on the QC-filtered CRISPR guide count matrix from the {@linkplain CellFilteringState}.
 * It wraps the [`logNormCounts`](https://www.jkanche.com/scran.js/global.html#logNormCounts) functions
 * from [**scran.js**](https://github.com/jkanche/scran.js).
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class CrisprNormalizationState {
    #filter;
    #parameters;
    #cache;

    constructor(filter, parameters = null, cache = null) {
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
    }

    /***************************
     ******** Getters **********
     ***************************/

    valid() {
        let filtered = this.#filter.fetchFilteredMatrix();
        return filtered.has("CRISPR");
    }

    /**
     * @return {ScranMatrix} A ScranMatrix object containing the normalized guide abundance values,
     * available after running {@linkcode AdtNormalizationState#compute compute}.
     */
    fetchNormalizedMatrix() {
        if (!("matrix" in this.#cache)) {
            this.#raw_compute();
        }
        return this.#cache.matrix;
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
        var mat = this.#filter.fetchFilteredMatrix().get("CRISPR");
        var block = this.#filter.fetchFilteredBlock();
        utils.freeCache(this.#cache.matrix);
        this.#cache.matrix = scran.logNormCounts(mat, { block: block, allowZeros: true });
        return;
    }

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     *
     * @return The object is updated with new results.
     */
    compute() {
        this.changed = false;

        if (this.#filter.changed) {
            if (this.valid()) {
                this.#raw_compute();
                this.changed = true;
            }
        }

        return;
    }

    static defaults() {
        return {};
    }

    /*************************
     ******** Saving *********
     *************************/

    serialize(handle) {
        let ghandle = handle.createGroup(step_name);
        let phandle = ghandle.createGroup("parameters"); 
        let rhandle = ghandle.createGroup("results"); 
    }
}

/**************************
 ******** Loading *********
 **************************/

export function unserialize(handle, filter) {
    let cache = {};
    let parameters = {};
    return new CrisprNormalizationState(filter, parameters, cache);
}
