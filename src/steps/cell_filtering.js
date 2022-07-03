import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as qcutils from "./utils/quality_control.js";
import * as inputs_module from "./inputs.js";

export const step_name = "cell_filtering";

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
            this.#raw_compute();
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
            this.#raw_compute();
        }
        return this.#cache.block_buffer;
    }

    fetchDiscards() {
        return this.#cache.discard_buffer;
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
        var discard = this.#cache.discard_buffer.array();
        return vec.filter((x, i) => !discard[i]);
    }

    /***************************
     ******** Compute **********
     ***************************/

    #raw_compute() {
        let to_use = utils.findValidUpstreamStates(this.#qc_states, "QC");
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

        // Subsetting all modalities.
        utils.freeCache(this.#cache.matrix);
        this.#cache.matrix = new scran.MultiMatrix;
        let available = this.#inputs.listAvailableTypes();
        for (const a of available) {
            let src = this.#inputs.fetchCountMatrix({ type: a });
            let sub = scran.filterCells(src, disc_buffer);
            this.#cache.matrix.add(a, sub);
        }

        // Filtering on the block.
        let block = this.#inputs.fetchBlock();
        if (block !== null) {
            let bcache = utils.allocateCachedArray(this.#cache.matrix.numberOfColumns(), "Int32Array", this.#cache, "block_buffer");
            scran.filterBlock(block, this.#cache.discard_buffer, { buffer: bcache });
        } else {
            utils.freeCache(this.#cache.block_buffer);
            this.#cache.block_buffer = null;
        }

        return;
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
            if (x.valid() && x.changed) {
                this.changed = true;
            }
        }

        if (this.changed) {
            this.#raw_compute();
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
        this.#cache.discard_buffer.forEach(x => { remaining += (x == 0); });
        return { "retained": remaining };
    }

    /*************************
     ******** Saving *********
     *************************/

    serialize(handle) {
        let ghandle = handle.createGroup(step_name);
        ghandle.createGroup("parameters"); // for tradition

        let rhandle = ghandle.createGroup("results"); 
        let disc = this.fetchDiscards();
        if (disc.owner === null) {
            // If it's not a view, we save it; otherwise, we skip it, under the
            // assumption that only one of the upstream QC states contains the
            // discard vector.
            rhandle.writeDataSet("discards", "Uint8", null, disc);
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
            // This only happens if there was only one upstream QC state; in which case, 
            // we figure out which upstream QC state contains the discard vector
            // and create a view on it so that our fetchDiscards() works properly.
            // (v1 and earlier also implicitly falls in this category.)

            let to_use = utils.findValidUpstreamStates(qc_states, "QC");
            if (to_use.length != 1) {
                throw new Error("only one upstream QC state should be valid if 'discards' is not available");
            }
            let use_state = qc_states[to_use[0]];
            cache.discard_buffer = use_state.fetchDiscards().view();
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
