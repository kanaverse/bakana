import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as qcutils from "./utils/quality_control.js";
import * as inputs_module from "./inputs.js";

export const step_name = "cell_filtering";

function findUsefulUpstreamStates(states, msg) {
    let is_valid = utils.findValidUpstreamStates(states, msg);
    let useful = [];
    for (const v of is_valid) {
        if (!states[v].skipped()) {
            useful.push(v);
        }
    }
    return useful;
}

/**
 * This step filters the count matrices to remove low-quality cells.
 * It wraps the `filterCells` function from [**scran.js**](https://github.com/jkanche/scran.js).
 * For multi-modal datasets, this combines the quality calls from all modalities; a cell is removed if it is considered low-quality in any individual modality.
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class CellFilteringState {
    #inputs;
    #qc_states;
    #cache;
    #parameters;

    constructor(inputs, qc_states, parameters = null, cache = null) {
        if (!(inputs instanceof inputs_module.InputsState)) {
            throw new Error("'inputs' should be an InputsState");
        }
        this.#inputs = inputs;

        for (const v of Object.values(qc_states)) {
            if (!(v instanceof qcutils.QualityControlStateBase)) {
                throw new Error("'qc_states' should contain QualityControlStateBase objects");
            }
        }
        this.#qc_states = qc_states;

        this.#parameters = (parameters === null ? {} : parameters);
        this.#cache = (cache === null ? {} : cache);
        this.changed = false;
    }

    free() {
        utils.freeCache(this.#cache.block_buffer);
        utils.freeCache(this.#cache.discard_buffer);
        utils.freeCache(this.#cache.matrix);
    }

    /***************************
     ******** Getters **********
     ***************************/

    hasAvailable(type) {
        if ("matrix" in this.#cache) {
            return this.#cache.matrix.has(type);
        } else {
            return this.#inputs.hasAvailable(type);
        }
    }

    fetchFilteredMatrix({ type = "RNA" } = {}) {
        if (!("matrix" in this.#cache)) {
            this.#raw_compute_matrix();
        }
        return this.#cache.matrix.get(type);
    }

    /**
     * Fetch the filtered vector of block assignments.
     *
     * @return An Int32WasmArray of length equal to the number of cells after filtering, containing the block assignment for each cell.
     *
     * Alternatively `null` if no blocks are present in the dataset.
     */
    fetchFilteredBlock() {
        if (!("block_buffer" in this.#cache)) {
            this.#raw_compute_block();
        }
        return this.#cache.block_buffer;
    }

    fetchDiscards() {
        if ("discard_buffer" in this.#cache) {
            return this.#cache.discard_buffer;
        } else {
            return null;
        }
    }

    /**
     * Fetch an annotation for the cells remaining after QC filtering.
     *
     * @param {string} col - Name of the annotation field of interest.
     *
     * @return {Array|TypedArray} Array of length equal to the number of filtered cells, containing the requested annotations.
     */
    fetchFilteredAnnotations(col) { 
        let vec = this.#inputs.fetchAnnotations(col);
        if ("discard_buffer" in this.#cache) {
            var discard = this.#cache.discard_buffer.array();
            return vec.filter((x, i) => !discard[i]);
        } else {
            return vec;
        }
    }

    /***************************
     ******** Compute **********
     ***************************/

    #raw_compute_matrix() {
        utils.freeCache(this.#cache.matrix);
        this.#cache.matrix = new scran.MultiMatrix;
        let available = this.#inputs.listAvailableTypes();
        for (const a of available) {
            let src = this.#inputs.fetchCountMatrix({ type: a });

            let sub;
            if ("discard_buffer" in this.#cache) {
                sub = scran.filterCells(src, this.#cache.discard_buffer);
            } else {
                sub = src.clone();
            }

            this.#cache.matrix.add(a, sub);
        }
    }

    #raw_compute_block() {
        utils.freeCache(this.#cache.block_buffer);

        let block = this.#inputs.fetchBlock();
        if (block !== null) {
            if ("discard_buffer" in this.#cache) {
                // Filtering on the block.
                let bcache = utils.allocateCachedArray(this.#cache.matrix.numberOfColumns(), "Int32Array", this.#cache, "block_buffer");
                scran.filterBlock(block, this.#cache.discard_buffer, { buffer: bcache });
            } else {
                this.#cache.block_buffer = block.view();
            }
        } else {
            this.#cache.block_buffer = null;
        }
    }

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     * Each argument is taken from the property of the same name in the `cell_filtering` property of the `parameters` of {@linkcode runAnalysis}.
     *
     * @return The object is updated with the new results.
     */
    compute() {
        this.changed = false;

        // Checking upstreams.
        if (this.#inputs.changed) {
            this.changed = true;
        }
        for (const x of Object.values(this.#qc_states)) {
            if (x.changed) {
                this.changed = true;
            }
        }

        if (this.changed) {
            let to_use = findUsefulUpstreamStates(this.#qc_states, "QC");

            if (to_use.length > 0) {
                let disc_buffer;
                let first = this.#qc_states[to_use[0]].fetchDiscards();

                if (to_use.length > 1) {
                    // A discard signal in any modality causes the cell to be removed. 
                    disc_buffer = utils.allocateCachedArray(first.length, "Uint8Array", this.#cache, "discard_buffer");
                    let disc_arr = disc_buffer.array();
                    disc_arr.fill(0);

                    for (const u of to_use) {
                        this.#qc_states[u].fetchDiscards().forEach((y, i) => { disc_arr[i] |= y; });
                    }
                } else {
                    // If there's only one valid modality, we just create a view on it
                    // to avoid unnecessary duplication.
                    disc_buffer = first.view();
                    utils.freeCache(this.#cache.discard_buffer);
                    this.#cache.discard_buffer = disc_buffer;
                }

            } else {
                // Deleting this so that serialization will behave correctly.
                utils.freeCache(this.#cache.discard_buffer);
                delete this.#cache.discard_buffer;
            }
        }
    }

    static defaults() {
        return {};
    }

    /***************************
     ******** Results **********
     ***************************/

    /**
     * Obtain a summary of the state, typically for display on a UI like **kana**.
     *
     * @return {Object} An object containing:
     *
     * - `retained`: the number of cells retained after filtering out low-quality cells.
     */
    summary() {
        let remaining = 0;
        if ("discard_buffer" in this.#cache) {
            this.#cache.discard_buffer.forEach(x => { remaining += (x == 0); });
        } else {
            let available = this.#inputs.hasAvailable();
            remaining = this.#inputs.fetchCountMatrix(available[0]).numberOfColumns();
        }
        return { "retained": remaining };
    }

    /*************************
     ******** Saving *********
     *************************/

    serialize(handle) {
        let ghandle = handle.createGroup(step_name);
        ghandle.createGroup("parameters"); // for tradition

        let rhandle = ghandle.createGroup("results"); 

        // If it's not present, then none of the upstream QC states applied
        // any filtering, so we don't really have anything to do here.
        if ("discard_buffer" in this.#cache) {
            let disc = this.#cache.discard_buffer;

            // If it's present and it's not a view, we save it; otherwise, we
            // skip it, under the assumption that only one of the upstream QC
            // states contains the discard vector.
            if (disc.owner === null) {
                rhandle.writeDataSet("discards", "Uint8", null, disc);
            }

        }
    }
}

export function unserialize(handle, inputs, qc_states) {
    let parameters = CellFilteringState.defaults();
    let cache = {};
    let output;

    try {
        if (step_name in handle.children) {
            let ghandle = handle.open(step_name);

            let rhandle = ghandle.open("results");
            if ("discards" in rhandle.children) {
                let discards = rhandle.open("discards", { load: true }).values; 
                cache.discard_buffer = scran.createUint8WasmArray(discards.length);
                cache.discard_buffer.set(discards);
            }
        } 

        if (!("discard_buffer" in cache)) {
            let to_use = findUsefulUpstreamStates(qc_states, "QC");

            if (to_use.length == 1) {
                // We figure out which upstream QC state contains the discard vector
                // and create a view on it so that our discard_buffer checks work properly.
                // (v1 and earlier also implicitly falls in this category.)
                let use_state = qc_states[to_use[0]];
                cache.discard_buffer = use_state.fetchDiscards().view();
            } else if (to_use.length == 0) {
                // No-op; we don't need to define discard_buffer.
                ;
            } else {
                throw new Error("no more than one upstream QC state should be valid if 'discards' is not available");
            }
        }

        output = new CellFilteringState(inputs, qc_states, parameters, cache);
    } catch (e) {
        utils.freeCache(cache.discard_buffer);
        utils.freeCache(output);
        throw e;
    }

    return {
        state: output,
        parameters: { ...parameters } // make a copy to avoid pass-by-reference links with state's internal parameters
    };
}
