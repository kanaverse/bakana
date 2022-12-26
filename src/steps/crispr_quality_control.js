import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as inputs_module from "./inputs.js";

export const step_name = "crispr_quality_control";

/**
 * Results of computing per-cell CRISPR-derived QC metrics,
 * see [here](https://www.jkanche.com/scran.js/PerCellCrisprQcMetricsResults.html) for details.
 *
 * @external PerCellCrisprQcMetricsResults
 */

/**
 * Suggested filters for the CRISPR-derived QC metrics,
 * see [here](https://www.jkanche.com/scran.js/SuggestCrisprQcFiltersResults.html) for details.
 *
 * @external SuggestCrisprQcFiltersResults
 */

/**
 * This step applies quality control on the CRISPR guide count matrix.
 * Specifically, it computes the QC metrics and filtering thresholds, 
 * wrapping the [`perCellCrisprQcMetrics`](https://www.jkanche.com/scran.js/global.html#perCellCrisprQcMetrics)
 * and [`suggestCrisprQcFilters`](https://www.jkanche.com/scran.js/global.html#suggestCrisprQcFilters) functions
 * from [**scran.js**](https://github.com/jkanche/scran.js).
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
            nmads: 3
        };
    }

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     * Each argument is taken from the property of the same name in the `crispr_quality_control` property of the `parameters` of {@linkcode runAnalysis}.
     *
     * @param {number} nmads - Number of MADs to use for automatically selecting the filter threshold on the maximum count. 
     *
     * @return The object is updated with the new results.
     */
    compute(igg_prefix, nmads, min_detected_drop) {
        this.changed = false;

        if (this.#inputs.changed || igg_prefix !== this.#parameters.igg_prefix) {
            utils.freeCache(this.#cache.metrics);

            if (this.valid()) {
                var mat = this.#inputs.fetchCountMatrix().get("CRISPR");
                this.#cache.metrics = scran.perCellCrisprQcMetrics(mat);
                this.changed = true;
            } else {
                delete this.#cache.metrics;
            }
        }

        if (this.changed || nmads !== this.#parameters.nmads) {
            utils.freeCache(this.#cache.filters);

            if (this.valid()) {
                let block = this.#inputs.fetchBlock();
                this.#cache.filters = scran.suggestCrisprQcFilters(this.#cache.metrics, { numberOfMADs: nmads, block: block });
                var discard = utils.allocateCachedArray(this.#cache.metrics.numberOfCells(), "Uint8Array", this.#cache, "discard_buffer");
                this.#cache.filters.filter(this.#cache.metrics, { block: block, buffer: discard });
                this.changed = true;
            } else {
                delete this.#cache.filters;
            }

            this.#parameters.nmads = nmads;
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
            phandle.writeDataSet("nmads", "Float64", [], this.#parameters.nmads);
        }

        {
            let rhandle = ghandle.createGroup("results"); 
            
            if (this.valid()) {
                let mhandle = rhandle.createGroup("metrics");
                mhandle.writeDataSet("sums", "Float64", null, this.#cache.metrics.sums({ copy: "view" }));
                mhandle.writeDataSet("detected", "Int32", null, this.#cache.metrics.detected({ copy: "view" }));
                mhandle.writeDataSet("max_proportion", "Float64", null, this.#cache.metrics.maxProportions({ copy: "view" }));
                mhandle.writeDataSet("max_index", "Float64", null, this.#cache.metrics.maxIndex({ copy: "view" }));

                let thandle = rhandle.createGroup("thresholds");
                thandle.writeDataSet("max_count", "Float64", null, this.#cache.filters.thresholdsMaxCount({ copy: "view" }));

                rhandle.writeDataSet("discards", "Uint8", null, this.#cache.discard_buffer);
            }
        }
    }
}

/**************************
 ******** Loading *********
 **************************/

export function unserialize(handle, inputs) {
    let cache = {};
    let parameters = CrisprQualityControlState.defaults();
    let output;

    if (step_name in handle.children) {
        let ghandle = handle.open(step_name);

        let phandle = ghandle.open("parameters"); 
        parameters.nmads = phandle.open("nmads", { load: true }).values[0];

        try {
            let rhandle = ghandle.open("results");

            if ("metrics" in rhandle.children) { // if skip=true or valid() is false, QC metrics may not be reported.
                let mhandle = rhandle.open("metrics");

                let detected = mhandle.open("detected", { load: true }).values;
                cache.metrics = scran.emptyPerCellCrisprQcMetricsResults(detected.length);
                cache.metrics.detected({ fillable: true }).set(detected);

                let sums = mhandle.open("sums", { load: true }).values;
                cache.metrics.sums({ fillable: true }).set(sums);

                let max_prop = mhandle.open("max_proportion", { load: true }).values;
                cache.metrics.maxProportions({ fillable: true }).set(max_prop);

                let max_index = mhandle.open("max_index", { load: true }).values;
                cache.metrics.maxIndex({ fillable: true }).set(max_index);
            }

            if ("thresholds" in rhandle.children) { // if skip=true or valid() is false, QC thresholds may not be reported.
                let discards = rhandle.open("discards", { load: true }).values; 
                cache.discard_buffer = scran.createUint8WasmArray(discards.length);
                cache.discard_buffer.set(discards);

                let thandle = rhandle.open("thresholds");
                let thresholds_max_count = thandle.open("max_count", { load: true }).values;

                cache.filters = scran.emptySuggestCrisprQcFiltersResults(thresholds_max_count.length);
                cache.filters.thresholdsMaxCount({ fillable: true }).set(thresholds_max_count);
            }

            output = new CrisprQualityControlState(inputs, parameters, cache);
        } catch (e) {
            utils.freeCache(cache.metrics);
            utils.freeCache(cache.filters)
            utils.freeCache(output);
            throw e;
        }
    }

    return output;
}
