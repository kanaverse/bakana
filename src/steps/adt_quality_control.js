import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as inputs_module from "./inputs.js";

export const step_name = "adt_quality_control";

/**
 * This step applies quality control on the ADT count matrix.
 * Specifically, it computes the QC metrics and filtering thresholds, wrapping `computePerCellAdtQcMetrics` and `computePerCellAdtQcFilters` from [**scran.js**](https://github.com/jkanche/scran.js).
 * Note that the actual filtering is done by {@linkcode CellFilteringState}.
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
    }

    free() {
        utils.freeCache(this.#cache.metrics);
        utils.freeCache(this.#cache.filters);
        utils.freeCache(this.#cache.metrics_buffer);
        utils.freeCache(this.#cache.discard_buffer);
    }

    /***************************
     ******** Getters **********
     ***************************/

    valid() {
        let input = this.#inputs.fetchCountMatrix();
        return input.has("ADT");
    }

    /**
     * @return {object} Object containing the parameters.
     */
    fetchParameters() {
        return { ...this.#parameters }; // avoid pass-by-reference links.
    }

    /**
     * @return {?PerCellAdtQcFiltersResults} Result of filtering on the ADT-derived QC metrics.
     * This is available after running {@linkcode AdtQualityControlState#compute compute}.
     */
    fetchFilters() {
        return this.#cache.filters;
    }

    /**
     * @return {Uint8WasmArray} Buffer containing the discard vector of length equal to the number of cells,
     * where each element is truthy if the corresponding cell is to be discarded.
     * This is available after running {@linkcode AdtQualityControlState#compute compute}.
     */
    fetchDiscards() {
        return this.#cache.discard_buffer;
    }

    /**
     * @return {PerCellAdtQcMetricsResults} ADT-derived QC metrics,
     * available after running {@linkcode AdtQualityControlState#compute compute}.
     */
    fetchMetrics() {
        return this.#cache.metrics;
    }

    /***************************
     ******** Compute **********
     ***************************/

    static defaults() {
        return {
            igg_prefix: "IgG",
            nmads: 3,
            min_detected_drop: 0.1
        };
    }

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     * Each argument is taken from the property of the same name in the `adt_quality_control` property of the `parameters` of {@linkcode runAnalysis}.
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

            var mat = this.#inputs.fetchCountMatrix().get("ADT");
            var gene_info = this.#inputs.fetchFeatureAnnotations()["ADT"];

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
            this.#parameters.igg_prefix = igg_prefix;
            this.changed = true;
        }

        if (this.changed || nmads !== this.#parameters.nmads || min_detected_drop !== this.#parameters.min_detected_drop) {
            utils.freeCache(this.#cache.filters);

            let block = this.#inputs.fetchBlock();
            this.#cache.filters = scran.computePerCellAdtQcFilters(this.#cache.metrics, { numberOfMADs: nmads, minDetectedDrop: min_detected_drop, block: block });
            var discard = utils.allocateCachedArray(mat.numberOfRows(), "Uint8Array", this.#cache, "discard_buffer");
            this.#cache.filters.filter(this.#cache.metrics, { block: block, buffer: discard });

            this.#parameters.nmads = nmads;
            this.#parameters.min_detected_drop = min_detected_drop;
            this.changed = true;
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
            phandle.writeDataSet("igg_prefix", "String", [], this.#parameters.igg_prefix);
            phandle.writeDataSet("nmads", "Float64", [], this.#parameters.nmads);
            phandle.writeDataSet("min_detected_drop", "Float64", [], this.#parameters.min_detected_drop);
        }

        {
            let rhandle = ghandle.createGroup("results"); 

            let mhandle = rhandle.createGroup("metrics");
            mhandle.writeDataSet("sums", "Float64", null, this.#cache.metrics.sums({ copy: "view" }));
            mhandle.writeDataSet("detected", "Int32", null, this.#cache.metrics.detected({ copy: "view" }));
            mhandle.writeDataSet("igg_total", "Float64", null, this.#cache.metrics.subsetTotals(0, { copy: "view" }));

            let thandle = rhandle.createGroup("thresholds");
            thandle.writeDataSet("detected", "Float64", null, this.#cache.filters.thresholdsDetected({ copy: "view" }));
            thandle.writeDataSet("igg_total", "Float64", null, this.#cache.filters.thresholdsSubsetTotals(0, { copy: "view" }));

            rhandle.writeDataSet("discards", "Uint8", null, this.#cache.discard_buffer);
        }
    }
}

/**************************
 ******** Loading *********
 **************************/

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

        try {
            let rhandle = ghandle.open("results");

            if ("metrics" in rhandle.children) { // if skip=true or valid() is false, QC metrics may not be reported.
                let mhandle = rhandle.open("metrics");

                let detected = mhandle.open("detected", { load: true }).values;
                cache.metrics = scran.emptyPerCellAdtQcMetricsResults(detected.length, 1);
                cache.metrics.detected({ fillable: true }).set(detected);

                let sums = mhandle.open("sums", { load: true }).values;
                cache.metrics.sums({ fillable: true }).set(sums);
                let igg_total = mhandle.open("igg_total", { load: true }).values;
                cache.metrics.subsetTotals(0, { fillable: true }).set(igg_total);
            }

            if ("thresholds" in rhandle.children) { // if skip=true or valid() is false, QC thresholds may not be reported.
                let discards = rhandle.open("discards", { load: true }).values; 
                cache.discard_buffer = scran.createUint8WasmArray(discards.length);
                cache.discard_buffer.set(discards);

                let thandle = rhandle.open("thresholds");
                let thresholds_detected = thandle.open("detected", { load: true }).values;
                let thresholds_igg_total = thandle.open("igg_total", { load: true }).values;

                cache.filters = scran.emptyPerCellAdtQcFiltersResults(1, thresholds_detected.length);
                cache.filters.thresholdsDetected({ fillable: true }).set(thresholds_detected);
                cache.filters.thresholdsSubsetTotals(0, { fillable: true }).set(thresholds_igg_total);
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
