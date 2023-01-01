import * as scran from "scran.js";
import * as utils from "./utils/general.js";
import * as rna_pca_module from "./rna_pca.js";
import * as adt_pca_module from "./adt_pca.js";
import * as crispr_pca_module from "./crispr_pca.js";

export const step_name = "combine_embeddings";

function find_nonzero_upstream_states(pca_states, weights) {
    let tmp = utils.findValidUpstreamStates(pca_states);
    let to_use = [];
    for (const k of tmp) {
        if (weights[k] > 0) {
            to_use.push(k);
        }
    }
    return to_use;
}

/**
 * This step combines multiple embeddings from different modalities into a single matrix for downstream analysis.
 * It wraps the [`scaleByNeighbors`](https://jkanche.com/scran.js/global.html#scaleByNeighbors) function
 * from [**scran.js**](https://jkanche.com/scran.js).
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class CombineEmbeddingsState {
    #pca_states;
    #parameters;
    #cache;

    constructor(pca_states, parameters = null, cache = null) {
        if (!(pca_states.RNA instanceof rna_pca_module.RnaPcaState)) {
            throw new Error("'pca_states.RNA' should be an RnaPcaState object");
        }
        if (!(pca_states.ADT instanceof adt_pca_module.AdtPcaState)) {
            throw new Error("'pca_states.ADT' should be an AdtPcaState object");
        }
        if (!(pca_states.CRISPR instanceof crispr_pca_module.CrisprPcaState)) {
            throw new Error("'pca_states.CRISPR' should be an CrisprPcaState object");
        }
        this.#pca_states = pca_states;

        this.#parameters = (parameters === null ? {} : parameters);
        this.#cache = (cache === null ? {} : cache);
        this.changed = false;
    }

    free() {
        utils.freeCache(this.#cache.combined_buffer);
    }

    /***************************
     ******** Getters **********
     ***************************/

    /**
     * @return {Float64WasmArray} Buffer containing the combined embeddings as a column-major dense matrix,
     * where the rows are the dimensions and the columns are the cells.
     * This is available after running {@linkcode CombineEmbeddingsState#compute compute}.
     */
    fetchCombined() {
        return this.#cache.combined_buffer;
    }

    /**
     * @return {number} Number of cells in {@linkcode CombineEmbeddingsState#fetchCombined fetchCombined},
     * available after running {@linkcode CombineEmbeddingsState#compute compute}.
     */
    fetchNumberOfCells() {
        return this.#cache.num_cells;
    }

    /**
     * @return {number} Number of dimensions in {@linkcode CombineEmbeddingsState#fetchCombined fetchCombined},
     * available after running {@linkcode CombineEmbeddingsState#compute compute}.
     */
    fetchNumberOfDimensions() {
        return this.#cache.total_dims;
    }

    /**
     * @return {object} Object containing the parameters.
     */
    fetchParameters() {
        // Avoid any pass-by-reference activity.
        return { ...this.#parameters };
    }

    /***************************
     ******** Compute **********
     ***************************/

    static defaults() {
        return { 
            rna_weight: 1,
            adt_weight: 1,
            crispr_weight: 0,
            approximate: true
        };
    }

    static createPcsView(cache, upstream) {
        utils.freeCache(cache.combined_buffer);
        cache.combined_buffer = upstream.principalComponents({ copy: "view" }).view();
        cache.num_cells = upstream.numberOfCells();
        cache.total_dims = upstream.numberOfPCs();
    }

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     *
     * @param {number} rna_weight - Relative weight of the RNA embeddings.
     * @param {number} adt_weight - Relative weight of the ADT embeddings.
     * @param {number} crispr_weight - Relative weight of the CRISPR embeddings.
     * @param {boolean} approximate - Whether an approximate nearest neighbor search should be used by `scaleByNeighbors`.
     *
     * @return The object is updated with new results.
     */
    compute(rna_weight, adt_weight, crispr_weight, approximate) {
        this.changed = false;

        for (const v of Object.values(this.#pca_states)) {
            if (v.changed) {
                this.changed = true;
                break;
            }
        }

        if (approximate !== this.#parameters.approximate) {
            this.#parameters.approximate = approximate;
            this.changed = true;
        }

        if (rna_weight !== this.#parameters.rna_weight || adt_weight !== this.#parameters.adt_weight || crispr_weight !== this.#parameters.crispr_weight) {
            this.#parameters.rna_weight = rna_weight;
            this.#parameters.adt_weight = adt_weight;
            this.#parameters.crispr_weight = crispr_weight;
            this.changed = true;
        }

        if (this.changed) { 
            const weights = { RNA: rna_weight, ADT: adt_weight, CRISPR: crispr_weight };
            let to_use = find_nonzero_upstream_states(this.#pca_states, weights);

            if (to_use.length > 1) {
                let weight_arr = to_use.map(x => weights[x]);
                let collected = [];
                let total = 0;
                let ncells = null;

                for (const k of to_use) {
                    let curpcs = this.#pca_states[k].fetchPCs();
                    collected.push(curpcs.principalComponents({ copy: "view" }));
                    if (ncells == null) {
                        ncells = curpcs.numberOfCells();
                    } else if (ncells !== curpcs.numberOfCells()) {
                        throw new Error("number of cells should be consistent across all embeddings");
                    }
                    total += curpcs.numberOfPCs();
                }

                let buffer = utils.allocateCachedArray(ncells * total, "Float64Array", this.#cache, "combined_buffer");
                scran.scaleByNeighbors(collected, ncells, { buffer: buffer, weights: weight_arr, approximate: approximate });
                this.#cache.num_cells = ncells;
                this.#cache.total_dims = total;

            } else {
                // If there's only one embedding, we shouldn't respond to changes
                // in parameters, because they won't have any effect.
                let pcs = this.#pca_states[to_use[0]].fetchPCs();
                this.constructor.createPcsView(this.#cache, pcs);
            }
        }

        // Updating all parameters anyway. This requires us to take ownership
        // of 'weights' to avoid pass-by-reference shenanigans.
        return;
    }

    /*************************
     ******** Saving *********
     *************************/

    serialize(handle) {
        let ghandle = handle.createGroup(step_name);

        {
            let phandle = ghandle.createGroup("parameters"); 
            phandle.writeDataSet("approximate", "Uint8", [], this.#parameters.approximate);
            phandle.writeDataSet("rna_weight", "Float64", [], this.#parameters.rna_weight);
            phandle.writeDataSet("adt_weight", "Float64", [], this.#parameters.adt_weight);
            phandle.writeDataSet("crispr_weight", "Float64", [], this.#parameters.crispr_weight);
        }

        {
            let rhandle = ghandle.createGroup("results");
            let pcs = this.fetchCombined();
            if (pcs.owner === null) {
                // If it's not a view, we save it; otherwise we assume
                // that we can recover it from the upstream PCA states.
                rhandle.writeDataSet(
                    "combined", 
                    "Float64", 
                    [this.fetchNumberOfCells(), this.fetchNumberOfDimensions()], // remember, it's transposed.
                    pcs
                ); 
            }
        }
    }
}

export function unserialize(handle, pca_states) {
    let cache = {};
    let parameters = CombineEmbeddingsState.defaults();
    let output;

    try {
        if (step_name in handle.children) {
            let ghandle = handle.open(step_name);

            {
                let phandle = ghandle.open("parameters");
                parameters.approximate = phandle.open("approximate", { load: true }).values[0] > 0;

                if ('weights' in phandle.children) { // v2.0
                    let whandle = phandle.open("weights");
                    if ("RNA" in whandle.children) {
                        parameters.rna_weight = whandle.open("RNA", { load: true }).values[0];
                    }
                    if ("ADT" in whandle.children) {
                        parameters.adt_weight = whandle.open("ADT", { load: true }).values[0];
                    }
                    if ("CRISPR" in whandle.children) {
                        parameters.crispr_weight = whandle.open("CRISPR", { load: true }).values[0];
                    }
                } else {
                    parameters.rna_weight = phandle.open("rna_weight", { load: true }).values[0];
                    parameters.adt_weight = phandle.open("adt_weight", { load: true }).values[0];
                    parameters.crispr_weight = phandle.open("crispr_weight", { load: true }).values[0];
                }
            }

            let rhandle = ghandle.open("results");

            if ("combined" in rhandle.children) {
                let phandle = rhandle.open("combined", { load: true });
                cache.num_cells = phandle.shape[0];
                cache.total_dims = phandle.shape[1];

                let vals = phandle.values;
                cache.combined_buffer = scran.createFloat64WasmArray(vals.length);
                cache.combined_buffer.set(vals);
            }
        }

        if (!("combined_buffer" in cache)) {
            // This only happens if there was only one upstream PCA state; in which case, 
            // we figure out which upstream PCA state contains the PC vector
            // and create a view on it so that our fetchPCs() works properly.
            // (v1 and earlier also implicitly falls in this category.)
            const weights = { RNA: parameters.rna_weight, ADT: parameters.adt_weight, CRISPR: parameters.crispr_weight };
            let to_use = find_nonzero_upstream_states(pca_states, weights);

            if (to_use.length != 1) {
                throw new Error("only one upstream PCA state should be valid with non-zero weight if 'combined' is not available");
            }

            let pcs = pca_states[to_use[0]].fetchPCs();
            CombineEmbeddingsState.createPcsView(cache, pcs);
        }

        output = new CombineEmbeddingsState(pca_states, parameters, cache);
    } catch (e) {
        utils.freeCache(cache.combined_buffer);
        utils.freeCache(output);
        throw e;
    }

    return output;
}
