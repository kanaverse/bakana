import * as scran from "scran.js"; 
import * as utils from "../utils/general.js";
import * as inputs_module from "../inputs.js";
import * as qcutils from "../utils/quality_control.js";

/**
 * This step applies quality control on the ADT count matrix.
 * This wraps `computePerCellAdtQcMetrics` and related functions from [**scran.js**](https://github.com/jkanche/scran.js).
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class AdtQualityControlState {
    #inputs;
    #cache;
    #parameters;

    constructor(inputs, parameters = null, cache = null) {
        if (!(inputs instanceof inputs_module.InputsState)) {
            throw new Error("'inputs' should be a State object from './inputs.js'");
        }
        this.#inputs = inputs;

        this.#parameters = (parameters === null ? {} : parameters);
        this.#cache = (cache === null ? {} : cache);
        this.changed = false;

        this.#parameters.__use_main_matrix__ = false;
    }

    free() {
        utils.freeCache(this.#cache.metrics);
        utils.freeCache(this.#cache.filters);
        utils.freeCache(this.#cache.block_buffer);
        utils.freeCache(this.#cache.metrics_buffer);
    }

    /***************************
     ******** Getters **********
     ***************************/

    fetchDiscards() {
        return this.#cache.filters.discardOverall({ copy: "view" });
    }

    /***************************
     ******** Compute **********
     ***************************/

    useMainMatrix() {
        this.#parameters.__use_main_matrix__ = true;
        return;
    }

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     * Each argument is taken from the property of the same name in the `quality_control` property of the `parameters` of {@linkcode runAnalysis}.
     *
     * @param {string} igg_prefix - Prefix of the identifiers for isotype controls.
     * @param {number} nmads - Number of MADs to use for automatically selecting the filter threshold for each metric.
     * @param {number} min_detected_drop - Minimum proportional drop in the number of detected features before a cell is to be considered low-quality.
     *
     * @return The object is updated with the new results.
     */
    compute(igg_prefix, nmads, min_detected_drop) {
        this.changed = false;

        if (this.#inputs.changed || igg_prefix !== this.#parameters.igg_prefix) {
            utils.freeCache(this.#cache.metrics);

            if (this.#inputs.hasAdt() || this.#parameters.__use_main_matrix__) {
                var mat = (this.#parameters.__use_main_matrix__ ? this.#inputs.fetchCountMatrix() : this.#inputs.fetchAdtCountMatrix());
                var gene_info = (this.#parameters.__use_main_matrix__ ? this.#inputs.fetchGenes() : this.#inputs.fetchAdtGenes());

                // Finding the prefix.
                var subsets = utils.allocateCachedArray(mat.numberOfRows(), "Uint8Array", this.#cache, "metrics_buffer");
                subsets.fill(0);
                var sub_arr = subsets.array();
                var lower_igg = igg_prefix.toLowerCase();
                for (const [key, val] of Object.entries(gene_info)) {
                    val.forEach((x, i) => { 
                        if (x.toLowerCase().startsWith(lower_igg)) {
                            sub_arr[i] = 1;                        
                        }
                    });
                }

                this.#cache.metrics = scran.computePerCellAdtQcMetrics(mat, [subsets]);
            } else {
                delete this.#cache.metrics;
            }

            this.#parameters.igg_prefix = igg_prefix;
            this.changed = true;
        }

        if (this.changed || nmads !== this.#parameters.nmads || min_detected_drop !== this.#parameters.min_detected_drop) {
            utils.freeCache(this.#cache.filters);

            if ("metrics" in this.#cache) {
                let block = this.#inputs.fetchBlock();
                this.#cache.filters = scran.computePerCellAdtQcFilters(this.#cache.metrics, { numberOfMADs: nmads, minDetectedDrop: min_detected_drop, block: block });
            } else {
                delete this.#cache.filters;
            }

            this.#parameters.nmads = nmads;
            this.#parameters.min_detected_drop = min_detected_drop;
            this.changed = true;
        }

        return;
    }

    /***************************
     ******** Results **********
     ***************************/

    #format_metrics({ copy = true } = {}) {
        return {
            sums: this.#cache.metrics.sums({ copy: copy }),
            detected: this.#cache.metrics.detected({ copy: copy }),
            igg_total: this.#cache.metrics.subsetTotals(0, { copy: copy })
        };
    }

    #format_thresholds({ copy = true } = {}) {
        return {
            detected: this.#cache.filters.thresholdsDetected({ copy: copy }),
            igg_total: this.#cache.filters.thresholdsSubsetTotals(0, { copy: copy })
        }
    }

    #valid() {
        return "metrics" in this.#cache;
    }

    /**
     * Obtain a summary of the state, typically for display on a UI like **kana**.
     *
     * @return An object containing:
     *
     * - `data`: an object containing one property for each sample.
     *   Each property is itself an object containing `sums`, `detected` and `igg_total`,
     *   which are TypedArrays containing the relevant QC metrics for all cells in that sample.
     * - `thresholds`: an object containing one property for each sample.
     *   Each property is itself an object containing `detected` and `igg_total`,
     *   which are numbers containing the thresholds on the corresponding QC metrics for that sample.
     * - `retained`: the number of cells remaining after QC filtering.
     */
    summary() {
        var output = null;

        if (this.#valid()) {
            output = {};

            var blocks = this.#inputs.fetchBlockLevels();
            if (blocks === null) {
                blocks = [ "default" ];
                output.data = { default: this.#format_metrics() };
            } else {
                let metrics = this.#format_metrics({ copy: "view" });
                let bids = this.#inputs.fetchBlock();
                output.data = qcutils.splitMetricsByBlock(metrics, blocks, bids);
            }

            let listed = this.#format_thresholds({ copy: "view" });
            output.thresholds = qcutils.splitThresholdsByBlock(listed, blocks);

            // We don't use sums for filtering but we do report it in the metrics,
            // so we just add some NaNs to the thresholds for consistency.
            for (const [k, v] of Object.entries(output.thresholds)) {
                v.sums = NaN;
            }
        }

        return output;
    }

    /*************************
     ******** Saving *********
     *************************/

    serialize(handle) {
        let ghandle = handle.createGroup("adt-quality_control");

        {
            let phandle = ghandle.createGroup("parameters"); 
            phandle.writeDataSet("igg_prefix", "String", [], this.#parameters.igg_prefix);
            phandle.writeDataSet("nmads", "Float64", [], this.#parameters.nmads);
            phandle.writeDataSet("min_detected_drop", "Float64", [], this.#parameters.min_detected_drop);
        }

        {
            let rhandle = ghandle.createGroup("results"); 

            if (this.#valid()) {
                {
                    let mhandle = rhandle.createGroup("metrics");
                    let data = this.#format_metrics({ copy: "view" });
                    mhandle.writeDataSet("sums", "Float64", null, data.sums)
                    mhandle.writeDataSet("detected", "Int32", null, data.detected);
                    mhandle.writeDataSet("igg_total", "Float64", null, data.igg_total);
                }

                {
                    let thandle = rhandle.createGroup("thresholds");
                    let thresholds = this.#format_thresholds({ copy: "hdf5" }); 
                    for (const x of [ "detected", "igg_total" ]) {
                        let current = thresholds[x];
                        thandle.writeDataSet(x, "Float64", null, current);
                    }
                }

                let disc = this.fetchDiscards();
                rhandle.writeDataSet("discards", "Uint8", null, disc);
            }
        }
    }
}

/**************************
 ******** Loading *********
 **************************/

class AdtQCFiltersMimic {
    constructor(detected, igg_total, discards) {
        this.detected_ = detected;
        this.igg_total_ = igg_total;
        this.discards = scran.createUint8WasmArray(discards.length);
        this.discards.set(discards);
    }

    thresholdsDetected({ copy }) {
        return utils.mimicGetter(this.detected_, copy);
    }

    thresholdsSubsetTotals(index, { copy }) {
        if (index != 0) {
            throw "only 'index = 0' is supported for mimics";
        }
        return utils.mimicGetter(this.igg_total_, copy);
    }

    discardOverall({ copy }) {
        return utils.mimicGetter(this.discards, copy);
    }

    free() {
        this.discards.free();
    }
}

export function unserialize(handle, inputs) {
    let ghandle = handle.open("adt-quality_control");

    let parameters = {};
    {
        let phandle = ghandle.open("parameters"); 
        parameters = {
            mito_prefix: phandle.open("igg_prefix", { load: true }).values[0],
            nmads: phandle.open("nmads", { load: true }).values[0],
            min_detected_drop: phandle.open("min_detected_drop", { load: true }).values[0]
        }
    }

    let cache = {};
    {
        let rhandle = ghandle.open("results");

        if ("metrics" in rhandle.children) {
            let mhandle = rhandle.open("metrics");
            let detected = mhandle.open("detected", { load: true }).values;
            cache.metrics = scran.emptyPerCellAdtQcMetricsResults(detected.length, 1);
            cache.metrics.detected({ copy: false }).set(detected);
            let igg_total = mhandle.open("igg_total", { load: true }).values;
            cache.metrics.subsetTotals(0, { copy: false }).set(igg_total);

            let thandle = rhandle.open("thresholds");
            let thresholds_detected = thandle.open("detected", { load: true }).values;
            let thresholds_igg_total = thandle.open("igg_total", { load: true }).values;

            let discards = rhandle.open("discards", { load: true }).values; 
            cache.filters = new QCFiltersMimic(
                thresholds_detected,
                thresholds_igg_total,
                discards
            );
        }
    }

    return {
        state: new AdtQualityControlState(inputs, parameters, cache),
        parameters: { ...parameters }
    };
}
