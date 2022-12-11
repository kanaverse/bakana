import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as inputs_module from "./inputs.js";
import * as qcutils from "./utils/quality_control.js";

export const step_name = "adt_quality_control";

/**
 * This step applies quality control on the ADT count matrix.
 * Specifically, it computes the QC metrics and filtering thresholds, wrapping `computePerCellAdtQcMetrics` and `computePerCellAdtQcFilters` from [**scran.js**](https://github.com/jkanche/scran.js).
 * Note that the actual filtering is done by {@linkcode CellFilteringState}.
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class AdtQualityControlState extends qcutils.QualityControlStateBase {
    #inputs;
    #cache;
    #parameters;

    constructor(inputs, parameters = null, cache = null) {
        super();

        if (!(inputs instanceof inputs_module.InputsState)) {
            throw new Error("'inputs' should be a State object from './inputs.js'");
        }
        this.#inputs = inputs;

        this.#parameters = (parameters === null ? {} : parameters);
        this.#cache = (cache === null ? {} : cache);
        this.changed = false;

        this.#parameters.target_matrix = "ADT";
    }

    free() {
        utils.freeCache(this.#cache.metrics);
        utils.freeCache(this.#cache.filters);
        utils.freeCache(this.#cache.metrics_buffer);
    }

    /***************************
     ******** Getters **********
     ***************************/

    valid() {
        let input = this.#inputs.fetchCountMatrix();
        return input.has(this.#parameters.target_matrix);
    }

    skipped() {
        return this.#parameters.skip;
    }

    /**
     * @return {object} Object containing the parameters.
     */
    fetchParameters() {
        let output = { ...this.#parameters }; // avoid pass-by-reference links.
        delete output.target_matrix;
        return output;
    }

    /**
     * @return {?PerCellAdtQcFiltersResults} Result of filtering on the ADT-derived QC metrics.
     * This is available after running {@linkcode AdtQualityControlState#compute compute}.
     * Alternatively `null`, if no ADTs are available or if this step was skipped.
     */
    fetchFilters() {
        if (this.valid() && !this.skipped()) {
            return this.#cache.filters;
        } else {
            return null;
        }
    }

    /**
     * @return {?Uint8WasmArray} Buffer containing the discard vector of length equal to the number of cells,
     * where each element is truthy if the corresponding cell is to be discarded.
     * This is available after running {@linkcode AdtQualityControlState#compute compute}.
     * Alternatively `null`, if no ADTs are available or if this step was skipped.
     */
    fetchDiscards() {
        if (this.valid() && !this.skipped()) {
            return this.#cache.filters.discardOverall({ copy: "view" });
        } else {
            return null;
        }
    }

    /**
     * @return {?PerCellAdtQcMetricsResults} ADT-derived QC metrics,
     * available after running {@linkcode AdtQualityControlState#compute compute}.
     * Alternatively `null`, if no ADTs are available or if this step was skipped.
     */
    fetchMetrics() {
        if (this.valid() && !this.skipped()) {
            return this.#cache.metrics;
        } else {
            return null;
        }
    }

    /***************************
     ******** Compute **********
     ***************************/

    useRNAMatrix() {
        // For testing only!
        this.#parameters.target_matrix = "RNA";
        return;
    }

    static defaults() {
        return {
            skip: false,
            igg_prefix: "IgG",
            nmads: 3,
            min_detected_drop: 0.1
        };
    }

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     * Each argument is taken from the property of the same name in the `adt_quality_control` property of the `parameters` of {@linkcode runAnalysis}.
     *
     * @param {boolean} skip - Whether to skip the calculation of quality control statistics.
     * @param {string} igg_prefix - Prefix of the identifiers for isotype controls.
     * @param {number} nmads - Number of MADs to use for automatically selecting the filter threshold for each metric.
     * @param {number} min_detected_drop - Minimum proportional drop in the number of detected features before a cell is to be considered low-quality.
     *
     * @return The object is updated with the new results.
     */
    compute(skip, igg_prefix, nmads, min_detected_drop) {
        this.changed = false;

        // If the metrics or filters aren't available and we're not skipping
        // this step, then we need to run this.
        let skip_impossible = skip || !this.valid();
        let unskip_metrics = (!skip_impossible && !("metrics" in this.#cache));
        let unskip_filters = (!skip_impossible && !("filters" in this.#cache));

        if (this.#inputs.changed || igg_prefix !== this.#parameters.igg_prefix || unskip_metrics) {
            utils.freeCache(this.#cache.metrics);

            if (skip_impossible) {
                // Delete anything existing, as it won't be valid as other
                // things have changed upstream of us. This ensures that we
                // can re-run this step later via unskip_metrics = true.
                delete this.#cache.metrics;
            } else {
                var mat = this.#inputs.fetchCountMatrix().get(this.#parameters.target_matrix);
                var gene_info = this.#inputs.fetchFeatureAnnotations()[this.#parameters.target_matrix];

                // Finding the prefix.
                var subsets = utils.allocateCachedArray(mat.numberOfRows(), "Uint8Array", this.#cache, "metrics_buffer");
                subsets.fill(0);
                var sub_arr = subsets.array();
                var lower_igg = igg_prefix.toLowerCase();
                for (const key of gene_info.columnNames()) {
                    let val = gene_info.column(key);
                    val.forEach((x, i) => { 
                        if (x.toLowerCase().startsWith(lower_igg)) {
                            sub_arr[i] = 1;                        
                        }
                    });
                }

                this.#cache.metrics = scran.computePerCellAdtQcMetrics(mat, [subsets]);
            }

            this.#parameters.igg_prefix = igg_prefix;
            this.changed = true;
        }

        if (this.changed || nmads !== this.#parameters.nmads || min_detected_drop !== this.#parameters.min_detected_drop || unskip_filters) {
            utils.freeCache(this.#cache.filters);

            if (skip_impossible) {
                delete this.#cache.filters;
            } else {
                let block = this.#inputs.fetchBlock();
                this.#cache.filters = scran.computePerCellAdtQcFilters(this.#cache.metrics, { numberOfMADs: nmads, minDetectedDrop: min_detected_drop, block: block });
            }

            this.#parameters.nmads = nmads;
            this.#parameters.min_detected_drop = min_detected_drop;
            this.changed = true;
        }

        if (this.#parameters.skip !== skip) {
            this.changed = true;
            this.#parameters.skip = skip;
        } else if (this.#parameters.skip && skip) {
            this.changed = false; // force this to be true, no matter what.
        }

        return;
    }

    /*************************
     ******** Saving *********
     *************************/

    serialize(handle) {
        let ghandle = handle.createGroup(step_name);

        {
            let phandle = ghandle.createGroup("parameters"); 
            phandle.writeDataSet("skip", "Uint8", [], Number(this.#parameters.skip));
            phandle.writeDataSet("igg_prefix", "String", [], this.#parameters.igg_prefix);
            phandle.writeDataSet("nmads", "Float64", [], this.#parameters.nmads);
            phandle.writeDataSet("min_detected_drop", "Float64", [], this.#parameters.min_detected_drop);
        }

        {
            let rhandle = ghandle.createGroup("results"); 

            if ("metrics" in this.#cache) { // if skip=true, metrics may not be computed... but if they are, we save them anyway.
                let mhandle = rhandle.createGroup("metrics");
                mhandle.writeDataSet("sums", "Float64", null, this.#cache.metrics.sums({ copy: "view" }));
                mhandle.writeDataSet("detected", "Int32", null, this.#cache.metrics.detected({ copy: "view" }));
                mhandle.writeDataSet("igg_total", "Float64", null, this.#cache.metrics.subsetTotals(0, { copy: "view" }));
            }

            if ("filters" in this.#cache) { // if skip=true, thresholds may not be reported... but if they are, we save them anyway.
                let thandle = rhandle.createGroup("thresholds");
                thandle.writeDataSet("detected", "Float64", null, this.#cache.filters.thresholdsDetected({ copy: "hdf5" }));
                thandle.writeDataSet("igg_total", "Float64", null, this.#cache.filters.thresholdsSubsetTotals(0, { copy: "hdf5" }));

                let disc = this.fetchDiscards();
                rhandle.writeDataSet("discards", "Uint8", null, disc);
            }
        }
    }
}

/**************************
 ******** Loading *********
 **************************/

class AdtQcFiltersMimic {
    constructor(detected, igg_total, discards) {
        this.detected_ = detected;
        this.igg_total_ = igg_total;
        try {
            this.discards = scran.createUint8WasmArray(discards.length);
            this.discards.set(discards);
        } catch (e) {
            utils.freeCache(this.discards);
            throw e;
        }
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
    let cache = {};
    let parameters = AdtQualityControlState.defaults();
    let output;

    if (step_name in handle.children) {
        let ghandle = handle.open(step_name);

        let phandle = ghandle.open("parameters"); 
        parameters.igg_prefix = phandle.open("igg_prefix", { load: true }).values[0];
        parameters.nmads = phandle.open("nmads", { load: true }).values[0];
        parameters.min_detected_drop = phandle.open("min_detected_drop", { load: true }).values[0];

        if ("skip" in phandle.children) {
            parameters.skip = phandle.open("skip", { load: true }).values[0] > 0;
        }

        try {
            let rhandle = ghandle.open("results");

            if ("metrics" in rhandle.children) { // if skip=true or valid() is false, QC metrics may not be reported.
                let mhandle = rhandle.open("metrics");

                let detected = mhandle.open("detected", { load: true }).values;
                cache.metrics = scran.emptyPerCellAdtQcMetricsResults(detected.length, 1);
                cache.metrics.detected({ copy: false }).set(detected);

                let sums = mhandle.open("sums", { load: true }).values;
                cache.metrics.sums({ copy: false }).set(sums);
                let igg_total = mhandle.open("igg_total", { load: true }).values;
                cache.metrics.subsetTotals(0, { copy: false }).set(igg_total);
            }

            if ("thresholds" in rhandle.children) { // if skip=true or valid() is false, QC thresholds may not be reported.
                let thandle = rhandle.open("thresholds");
                let thresholds_detected = thandle.open("detected", { load: true }).values;
                let thresholds_igg_total = thandle.open("igg_total", { load: true }).values;

                let discards = rhandle.open("discards", { load: true }).values; 
                cache.filters = new AdtQcFiltersMimic(
                    thresholds_detected,
                    thresholds_igg_total,
                    discards
                );
            }

            output = new AdtQualityControlState(inputs, parameters, cache);
        } catch (e) {
            utils.freeCache(cache.metrics);
            utils.freeCache(cache.filters)
            utils.freeCache(output);
            throw e;
        }
    } else {
        // Fallback for v1.
        output = new AdtQualityControlState(inputs, parameters, cache);
    }

    return output;
}
