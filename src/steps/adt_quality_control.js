import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as inputs_module from "./inputs.js";

export const step_name = "adt_quality_control";

/**
 * Results of computing per-cell ADT-derived QC metrics,
 * see [here](https://kanaverse.github.io/scran.js/PerCellAdtQcMetricsResults.html) for details.
 *
 * @external PerCellAdtQcMetricsResults
 */

/**
 * Suggested filters for the ADT-derived QC metrics,
 * see [here](https://kanaverse.github.io/scran.js/SuggestAdtQcFiltersResults.html) for details.
 *
 * @external SuggestAdtQcFiltersResults
 */

/**
 * This step applies quality control on the ADT count matrix.
 * Specifically, it computes the QC metrics and filtering thresholds, 
 * wrapping the [`perCellAdtQcMetrics`](https://kanaverse.github.io/scran.js/global.html#perCellAdtQcMetrics)
 * and [`suggestAdtQcFilters`](https://kanaverse.github.io/scran.js/global.html#suggestAdtQcFilters) functions
 * from [**scran.js**](https://github.com/kanaverse/scran.js).
 * Note that the actual filtering is done by {@linkplain CellFilteringState}.
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
     * @return {external:SuggestAdtQcFiltersResults} Result of filtering on the ADT-derived QC metrics.
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
     * @return {external:PerCellAdtQcMetricsResults} ADT-derived QC metrics,
     * available after running {@linkcode AdtQualityControlState#compute compute}.
     */
    fetchMetrics() {
        return this.#cache.metrics;
    }

    /****************************
     ******** Defaults **********
     ****************************/

    static defaults() {
        return {
            automatic: true,
            tag_id_column: null,
            igg_prefix: "IgG",
            nmads: 3,
            min_detected_drop: 0.1
        };
    }

    static configureFeatureParameters(lower_igg, annotations) {
        let counter = val => {
            let n = 0;
            val.forEach(x => {
                if (x.toLowerCase().startsWith(lower_igg)) {
                    n++;
                }
            });
            return n;
        };

        let best_key = null;
        let best = 0;

        let rn = annotations.rowNames();
        if (rn !== null) {
            best = counter(rn);
        }

        for (const key of annotations.columnNames()) {
            let latest = counter(annotations.column(key));
            if (latest > best) {
                best_key = key;
                best = latest;
            }
        }

        return best_key;
    }

    /***************************
     ******** Compute **********
     ***************************/

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     * 
     * @param {object} parameters - Parameter object, equivalent to the `adt_quality_control` property of the `parameters` of {@linkcode runAnalysis}.
     * @param {boolean} parameters.automatic - Automatically choose feature-based parameters based on the feature annotations. 
     * Specifically, `tag_id_column` is set to the column with the most matches to `igg_prefix`.
     * @param {?(string|number)} parameters.tag_id_column - Name or index of the column of the feature annotations that contains the tag identifiers.
     * If `null`, the row names are used.
     * Ignored if `automatic = true`.
     * @param {?string} parameters.igg_prefix - Prefix of the identifiers for isotype controls.
     * If `null`, no prefix-based identification is performed.
     * @param {number} parameters.nmads - Number of MADs to use for automatically selecting the filter threshold for each metric.
     * @param {number} parameters.min_detected_drop - Minimum proportional drop in the number of detected features before a cell is to be considered low-quality.
     *
     * @return The object is updated with the new results.
     */
    compute(parameters) {
        let { igg_prefix, nmads, min_detected_drop } = parameters;
        this.changed = false;

        let automatic;
        let tag_id_column; 
        if ("automatic" in parameters) {
            automatic = parameters.automatic;
            tag_id_column = parameters.tag_id_column;
        } else {
            automatic = true;
            tag_id_column = null;
        }

        if (
            this.#inputs.changed || 
            automatic !== this.#parameters.automatic ||
            igg_prefix !== this.#parameters.igg_prefix ||
            (!automatic && tag_id_column !== this.#parameters.tag_id_column)
        ) {
            utils.freeCache(this.#cache.metrics);

            if (this.valid()) {
                var tag_info = this.#inputs.fetchFeatureAnnotations()["ADT"];
                var subsets = utils.allocateCachedArray(tag_info.numberOfRows(), "Uint8Array", this.#cache, "metrics_buffer");
                subsets.fill(0);

                if (igg_prefix !== null) {
                    var lower_igg = igg_prefix.toLowerCase();
                    let key = tag_id_column;
                    if (automatic) {
                        key = AdtQualityControlState.configureFeatureParameters(lower_igg, tag_info);
                    }

                    let val = (key == null ? tag_info.rowNames() : tag_info.column(key));
                    if (val !== null) {
                        var sub_arr = subsets.array();
                        val.forEach((x, i) => { 
                            if (x.toLowerCase().startsWith(lower_igg)) {
                                sub_arr[i] = 1;                        
                            }
                        });
                    }
                }

                var mat = this.#inputs.fetchCountMatrix().get("ADT");
                this.#cache.metrics = scran.perCellAdtQcMetrics(mat, [subsets]);
                this.changed = true;
            } else {
                delete this.#cache.metrics;
            }
        }

        this.#parameters.automatic = automatic;
        this.#parameters.tag_id_column = tag_id_column;
        this.#parameters.igg_prefix = igg_prefix;

        if (this.changed || nmads !== this.#parameters.nmads || min_detected_drop !== this.#parameters.min_detected_drop) {
            utils.freeCache(this.#cache.filters);

            if (this.valid()) {
                let block = this.#inputs.fetchBlock();
                this.#cache.filters = scran.suggestAdtQcFilters(this.#cache.metrics, { numberOfMADs: nmads, minDetectedDrop: min_detected_drop, block: block });
                var discard = utils.allocateCachedArray(this.#cache.metrics.numberOfCells(), "Uint8Array", this.#cache, "discard_buffer");
                this.#cache.filters.filter(this.#cache.metrics, { block: block, buffer: discard });
                this.changed = true;
            } else {
                delete this.#cache.filters;
            }

            this.#parameters.nmads = nmads;
            this.#parameters.min_detected_drop = min_detected_drop;
        }

        return;
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

                cache.filters = scran.emptySuggestAdtQcFiltersResults(1, thresholds_detected.length);
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
