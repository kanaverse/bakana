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
        utils.freeCache(this.#cache.keep_buffer);
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
     * @return {Uint8WasmArray} Buffer containing a vector of length equal to the number of cells,
     * where each element is truthy if the corresponding cell is to be retained after filtering.
     * This is available after running {@linkcode AdtQualityControlState#compute compute}.
     */
    fetchKeep() {
        return this.#cache.keep_buffer;
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

    /**
     * @return {object} Object containing default parameters,
     * see the `parameters` argument in {@linkcode AdtQualityControlState#compute compute} for details.
     */
    static defaults() {
        return {
            guess_ids: true,
            tag_id_column: null,
            igg_prefix: "IgG",

            filter_strategy: "automatic", 
            nmads: 3,
            min_detected_drop: 0.1,

            detected_threshold: 0,
            igg_threshold: 1
        };
    }

    static #configureFeatureParameters(lower_igg, annotations) {
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
     * @param {boolean} [parameters.guess_ids] - Automatically choose feature-based parameters based on the feature annotations. 
     * Specifically, `tag_id_column` is set to the column with the most matches to `igg_prefix`.
     * @param {?(string|number)} [parameters.tag_id_column] - Name or index of the column of the feature annotations that contains the tag identifiers.
     * If `null`, the row names are used.
     * Ignored if `guess_ids = true`.
     * @param {?string} [parameters.igg_prefix]  - Prefix of the identifiers for isotype controls.
     * If `null`, no prefix-based identification is performed.
     * @param {string} [parameters.filter_strategy] - Strategy for defining a filter threshold for the QC metrics.
     * This can be `"automatic"` or `"manual"`.
     * @param {number} [parameters.nmads] - Number of MADs to use for automatically selecting the filter threshold for each metric.
     * Only used when `filter_strategy = "automatic"`.
     * @param {number} [parameters.min_detected_drop] - Minimum proportional drop in the number of detected features before a cell is to be considered low-quality.
     * Only used when `filter_strategy = "automatic"`.
     * @param {number} [parameters.detected_threshold] - Manual threshold on the detected number of features for each cell.
     * Cells are only retained if the detected number is equal to or greater than this threshold.
     * Only used when `filter_strategy = "manual"`.
     * @param {number} [parameters.igg_threshold] - Manual threshold on the isotype control totals for each cell.
     * Cells are only retained if their totals are less than or equal to this threshold.
     * Only used when `filter_strategy = "manual"`.
     *
     * @return The object is updated with the new results.
     */
    compute(parameters) {
        parameters = utils.defaultizeParameters(parameters, AdtQualityControlState.defaults(), [ "automatic" ]);
        this.changed = false;

        // Some back-compatibility here.
        if (typeof parameters.guess_ids === "undefined") {
            if ("automatic" in parameters) {
                parameters.guess_ids = parameters.automatic;
            } else {
                parameters.guess_ids = true;
            }
        }

        if (
            this.#inputs.changed || 
            parameters.guess_ids !== this.#parameters.guess_ids ||
            parameters.igg_prefix !== this.#parameters.igg_prefix ||
            (!parameters.guess_ids && parameters.tag_id_column !== this.#parameters.tag_id_column)
        ) {
            utils.freeCache(this.#cache.metrics);

            if (this.valid()) {
                var tag_info = this.#inputs.fetchFeatureAnnotations()["ADT"];
                var subsets = utils.allocateCachedArray(tag_info.numberOfRows(), "Uint8Array", this.#cache, "metrics_buffer");
                subsets.fill(0);

                if (parameters.igg_prefix !== null) {
                    var lower_igg = parameters.igg_prefix.toLowerCase();
                    let key = parameters.tag_id_column;
                    if (parameters.guess_ids) {
                        key = AdtQualityControlState.#configureFeatureParameters(lower_igg, tag_info);
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

        if (this.changed || 
            parameters.filter_strategy !== this.#parameters.filter_strategy ||
            parameters.nmads !== this.#parameters.nmads || 
            parameters.min_detected_drop !== this.#parameters.min_detected_drop ||
            parameters.detected_threshold !== this.#parameters.detected_threshold ||
            parameters.igg_threshold !== this.#parameters.igg_threshold
        ) {
            utils.freeCache(this.#cache.filters);

            if (this.valid()) {
                let block = this.#inputs.fetchBlock();

                if (parameters.filter_strategy === "automatic") {
                    this.#cache.filters = scran.suggestAdtQcFilters(this.#cache.metrics, { numberOfMADs: parameters.nmads, block: block });
                } else if (parameters.filter_strategy === "manual") {
                    let block_levels = this.#inputs.fetchBlockLevels();
                    this.#cache.filters = scran.emptySuggestAdtQcFiltersResults(1, block_levels === null ? 1 : block_levels.length);
                    this.#cache.filters.detected({ copy: false }).fill(parameters.detected_threshold);
                    this.#cache.filters.subsetSum(0, { copy: false }).fill(parameters.igg_threshold);
                } else {
                    throw new Error("unknown ADT QC filtering strategy '" + filter_strategy + "'");
                }

                var keep = utils.allocateCachedArray(this.#cache.metrics.numberOfCells(), "Uint8Array", this.#cache, "keep_buffer");
                this.#cache.filters.filter(this.#cache.metrics, { block: block, buffer: keep });
                this.changed = true;
            } else {
                delete this.#cache.filters;
            }
        }

        this.#parameters = parameters;
        return;
    }
}
