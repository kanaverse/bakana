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

            filter_strategy: "automatic", 
            nmads: 3,
            min_detected_drop: 0.1,

            detected_threshold: 0,
            igg_threshold: 1
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
     * @param {string} parameters.filter_strategy - Strategy for defining a filter threshold for the QC metrics.
     * This can be `"automatic"` or `"manual"`.
     * @param {number} parameters.nmads - Number of MADs to use for automatically selecting the filter threshold for each metric.
     * Only used when `filter_strategy = "automatic"`.
     * @param {number} parameters.min_detected_drop - Minimum proportional drop in the number of detected features before a cell is to be considered low-quality.
     * Only used when `filter_strategy = "automatic"`.
     * @param {number} parameters.detected_threshold - Manual threshold on the detected number of features for each cell.
     * Cells are only retained if the detected number is equal to or greater than this threshold.
     * Only used when `filter_strategy = "manual"`.
     * @param {number} parameters.igg_threshold - Manual threshold on the isotype control totals for each cell.
     * Cells are only retained if their totals are less than or equal to this threshold.
     * Only used when `filter_strategy = "manual"`.
     *
     * @return The object is updated with the new results.
     */
    compute(parameters) {
        let { igg_prefix, filter_strategy, nmads, min_detected_drop, detected_threshold, igg_threshold  } = parameters;
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

        if (this.changed || 
            filter_strategy !== this.#parameters.filter_strategy ||
            nmads !== this.#parameters.nmads || 
            min_detected_drop !== this.#parameters.min_detected_drop ||
            detected_threshold !== this.#parameters.detected_threshold ||
            igg_threshold !== this.#parameters.igg_threshold
        ) {
            utils.freeCache(this.#cache.filters);

            if (this.valid()) {
                let block = this.#inputs.fetchBlock();

                if (filter_strategy === "automatic") {
                    this.#cache.filters = scran.suggestAdtQcFilters(this.#cache.metrics, { numberOfMADs: nmads, block: block });
                } else if (filter_strategy === "manual") {
                    let block_levels = this.#inputs.fetchBlockLevels();
                    this.#cache.filters = scran.emptySuggestAdtQcFiltersResults(1, block_levels === null ? 1 : block_levels.length);
                    this.#cache.filters.thresholdsDetected({ copy: false }).fill(detected_threshold);
                    this.#cache.filters.thresholdsSubsetTotals(0, { copy: false }).fill(igg_threshold);
                } else {
                    throw new Error("unknown ADT QC filtering strategy '" + filter_strategy + "'");
                }

                var discard = utils.allocateCachedArray(this.#cache.metrics.numberOfCells(), "Uint8Array", this.#cache, "discard_buffer");
                this.#cache.filters.filter(this.#cache.metrics, { block: block, buffer: discard });
                this.changed = true;
            } else {
                delete this.#cache.filters;
            }

            this.#parameters.filter_strategy = filter_strategy;
            this.#parameters.nmads = nmads;
            this.#parameters.min_detected_drop = min_detected_drop;
            this.#parameters.detected_threshold = detected_threshold;
            this.#parameters.igg_threshold = igg_threshold;
        }

        return;
    }
}
