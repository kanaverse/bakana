import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as inputs_module from "./inputs.js";

export const step_name = "crispr_quality_control";

/**
 * Results of computing per-cell CRISPR-derived QC metrics,
 * see [here](https://kanaverse.github.io/scran.js/PerCellCrisprQcMetricsResults.html) for details.
 *
 * @external PerCellCrisprQcMetricsResults
 */

/**
 * Suggested filters for the CRISPR-derived QC metrics,
 * see [here](https://kanaverse.github.io/scran.js/SuggestCrisprQcFiltersResults.html) for details.
 *
 * @external SuggestCrisprQcFiltersResults
 */

/**
 * This step applies quality control on the CRISPR guide count matrix.
 * Specifically, it computes the QC metrics and filtering thresholds, 
 * wrapping the [`perCellCrisprQcMetrics`](https://kanaverse.github.io/scran.js/global.html#perCellCrisprQcMetrics)
 * and [`suggestCrisprQcFilters`](https://kanaverse.github.io/scran.js/global.html#suggestCrisprQcFilters) functions
 * from [**scran.js**](https://github.com/kanaverse/scran.js).
 * Note that the actual filtering is done by {@linkplain CellFilteringState}.
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class CrisprQualityControlState {
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
        return input.has("CRISPR");
    }

    /**
     * @return {object} Object containing the parameters.
     */
    fetchParameters() {
        return { ...this.#parameters }; // avoid pass-by-reference links.
    }

    /**
     * @return {external:SuggestCrisprQcFiltersResults} Result of filtering on the CRISPR-derived QC metrics.
     * This is available after running {@linkcode CrisprQualityControlState#compute compute}.
     */
    fetchFilters() {
        return this.#cache.filters;
    }

    /**
     * @return {Uint8WasmArray} Buffer containing the discard vector of length equal to the number of cells,
     * where each element is truthy if the corresponding cell is to be discarded.
     * This is available after running {@linkcode CrisprQualityControlState#compute compute}.
     */
    fetchDiscards() {
        return this.#cache.discard_buffer;
    }

    /**
     * @return {external:PerCellCrisprQcMetricsResults} CRISPR-derived QC metrics,
     * available after running {@linkcode CrisprQualityControlState#compute compute}.
     */
    fetchMetrics() {
        return this.#cache.metrics;
    }

    /***************************
     ******** Compute **********
     ***************************/

    static defaults() {
        return {
            filter_strategy: "automatic",
            nmads: 3,

            max_threshold: 0
        };
    }

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     *
     * @param {object} parameters - Parameter object, equivalent to the `crispr_quality_control` property of the `parameters` of {@linkcode runAnalysis}.
     * @param {string} parameters.filter_strategy - Strategy for defining a filter threshold for the QC metrics.
     * This can be `"automatic"` or `"manual"`.
     * @param {number} parameters.nmads - Number of MADs to use for automatically selecting the filter threshold on the maximum count. 
     * Only used when `filter_strategy = "automatic"`.
     * @param {number} parameters.max_threshold - Manual threshold on the maximum count for each cell.
     * Cells are only retained if their maximums are greater than or equal to this threshold.
     * Only used when `filter_strategy = "manual"`.
     *
     * @return The object is updated with the new results.
     */
    compute(parameters) {
        let { filter_strategy, nmads, max_threshold } = parameters;
        this.changed = false;

        if (this.#inputs.changed) {
            utils.freeCache(this.#cache.metrics);

            if (this.valid()) {
                var mat = this.#inputs.fetchCountMatrix().get("CRISPR");
                this.#cache.metrics = scran.perCellCrisprQcMetrics(mat);
                this.changed = true;
            } else {
                delete this.#cache.metrics;
            }
        }

        if (this.changed || 
            filter_strategy !== this.#parameters.filter_strategy ||
            nmads !== this.#parameters.nmads ||
            max_threshold !== this.#parameters.max_threshold
        ) {
            utils.freeCache(this.#cache.filters);

            if (this.valid()) {
                let block = this.#inputs.fetchBlock();

                if (filter_strategy === "automatic") {
                    this.#cache.filters = scran.suggestCrisprQcFilters(this.#cache.metrics, { numberOfMADs: nmads, block: block });
                } else if (filter_strategy === "manual") {
                    let block_levels = this.#inputs.fetchBlockLevels();
                    this.#cache.filters = scran.emptySuggestCrisprQcFiltersResults(block_levels === null ? 1 : block_levels.length);
                    this.#cache.filters.thresholdsMaxCount({ copy: false }).fill(max_threshold);
                } else {
                    throw new Error("unknown CRISPR QC filter strategy '" + filter_strategy + "'");
                }

                var discard = utils.allocateCachedArray(this.#cache.metrics.numberOfCells(), "Uint8Array", this.#cache, "discard_buffer");
                this.#cache.filters.filter(this.#cache.metrics, { block: block, buffer: discard });
                this.changed = true;
            } else {
                delete this.#cache.filters;
            }

            this.#parameters.filter_strategy = filter_strategy;
            this.#parameters.nmads = nmads;
            this.#parameters.max_threshold = max_threshold;
        }

        return;
    }
}
