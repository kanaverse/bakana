import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as rutils from "./utils/reader.js";
import * as qcutils from "./utils/quality_control.js";
import * as inputs_module from "./inputs.js";

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
        return this.#cache.matrix.has(type);
    }

    fetchFilteredMatrix({ type = "RNA" } = {}) {
        if (!("matrix" in this.#cache)) {
            this.#apply_filters();
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
        if (!("blocked" in this.#cache)) {
            this.#apply_filters();
        }
        if (this.#cache.blocked) {
            return this.#cache.block_buffer;
        } else {
            return null;
        }
    }

    fetchDiscards() {
        return this.#cache.discard_buffer;
    }

    /**
     * Fetch an annotation for the cells remaining after QC filtering.
     *
     * @param {string} col - Name of the annotation field of interest.
     *
     * @return An object similar to that returned by {@linkcode InputsState#fetchAnnotations InputsState.fetchAnnotations},
     * after filtering the vectors to only contain information for the remaining cells.
     * - For `type: "factor"`, the array in `index` will be filtered.
     * - For `type: "array"`, the array in `values` will be filtered.
     */
    fetchFilteredAnnotations(col) { 
        let vec = this.#inputs.fetchAnnotations(col);
        var discard = this.#cache.discard_buffer.array();
        let filterfun = (x, i) => !discard[i];
        if (vec.type === "factor") {
            vec.index = vec.index.filter(filterfun);
        } else {
            vec.values = vec.values.filter(filterfun);
        }
        return vec;
    }

    /***************************
     ******** Compute **********
     ***************************/

    #apply_filters() {
        let available = this.#inputs.listAvailableTypes();

        // A discard signal in any modality causes the cell to be removed.
        let mat = this.#inputs.fetchCountMatrix(available[0]);
        let disc_buffer = utils.allocateCachedArray(mat.numberOfColumns(), "Uint8Array", this.#cache, "discard_buffer");
        let disc_arr = disc_buffer.array();
        disc_arr.fill(0);

        for (const x of Object.values(this.#qc_states)) {
            if (x.valid()) {
                x.fetchDiscards().forEach((y, i) => { disc_arr[i] |= y; });
            }
        }

        // Subsetting all modalities.
        utils.freeCache(this.#cache.matrix);
        this.#cache.matrix = new rutils.MultiMatrix;
        for (const a of available) {
            let src = this.#inputs.fetchCountMatrix({ mode: a });
            let sub = scran.filterCells(src, disc_buffer);
            this.#cache.matrix.add(a, sub);
        }

        // Also subsetting the blocking factor, if any.
        let block = this.#inputs.fetchBlock();
        this.#cache.blocked = (block !== null);
        if (this.#cache.blocked) {
            let bcache = utils.allocateCachedArray(this.#cache.matrix.numberOfColumns(), "Int32Array", this.#cache, "block_buffer");
            let bcache_arr = bcache.array();
            let block_arr = block.array();
            let j = 0;
            for (let i = 0; i < block_arr.length; i++) {
                if (disc_arr[i] == 0) {
                    bcache_arr[j] = block_arr[i];
                    j++;
                }
            }
        }

        return;
    }

    compute() {
        if (this.#inputs.changed) {
            this.changed = true;
        }

        for (const x of Object.values(this.#qc_states)) {
            if (x.valid() && x.changed) {
                this.changed = true;
            }
        }

        if (this.changed) {
            this.#apply_filters();
        }
    }

    /***************************
     ******** Results **********
     ***************************/

    summary() {
        let remaining = 0;
        this.#cache.discard_buffer.forEach(x => { remaining += (x == 0); });
        return { "retained": remaining };
    }
}
