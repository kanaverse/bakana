import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import { mito } from "./mito.js";
import * as inputs_module from "./inputs.js";

export const step_name = "rna_quality_control";

/**
 * Results of computing per-cell RNA-derived QC metrics,
 * see [here](https://www.jkanche.com/scran.js/PerCellRnaQcMetricsResults.html) for details.
 *
 * @external PerCellRnaQcMetricsResults
 */

/**
 * Suggested filters for the RNA-derived QC metrics,
 * see [here](https://www.jkanche.com/scran.js/SuggestRnaQcFiltersResults.html) for details.
 *
 * @external SuggestRnaQcFiltersResults
 */

/**
 * This step applies quality control on the RNA count matrix.
 * Specifically, it computes the QC metrics and filtering thresholds, 
 * wrapping the [`perCellRnaQcMetrics`](https://www.jkanche.com/scran.js/global.html#perCellRnaQcMetrics)
 * and [`suggestRnaQcFilters`](https://www.jkanche.com/scran.js/global.html#suggestRnaQcFilters) functions
 * from [**scran.js**](https://github.com/jkanche/scran.js).
 * Note that the actual filtering is done by {@linkplain CellFilteringState}.
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class RnaQualityControlState {
    #inputs;
    #cache;
    #parameters;

    constructor(inputs, parameters = null, cache = null) {
        if (!(inputs instanceof inputs_module.InputsState)) {
            throw new Error("'inputs' should be an InputsState object");
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
        return input.has("RNA");
    }

    /**
     * @return {object} Object containing the parameters.
     */
    fetchParameters() {
        return { ...this.#parameters }; // avoid pass-by-reference links.
    }

    /**
     * @return {Uint8WasmArray} Buffer containing the discard vector of length equal to the number of cells,
     * where each element is truthy if the corresponding cell is to be discarded.
     */
    fetchDiscards() {
        return this.#cache.discard_buffer;
    }

    /**
     * @return {external:SuggestRnaQcFiltersResults} Result of filtering on the RNA-derived QC metrics.
     */
    fetchFilters() {
        return this.#cache.filters;
    }

    /**
     * @return {external:PerCellRnaQcMetricsResults} RNA-derived QC metrics.
     */
    fetchMetrics() {
        return this.#cache.metrics;
    }

    /***************************
     ******** Compute **********
     ***************************/

    static defaults () {
        return {
            use_mito_default: true,
            mito_prefix: "mt-",
            nmads: 3
        };
    }

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     * Each argument is taken from the property of the same name in the `quality_control` property of the `parameters` of {@linkcode runAnalysis}.
     *
     * @param {boolean} use_mito_default - Whether the internal mitochondrial gene lists should be used.
     * @param {string} mito_prefix - Prefix of the identifiers for mitochondrial genes, when `use_mito_default = false`.
     * @param {number} nmads - Number of MADs to use for automatically selecting the filter threshold for each metric.
     *
     * @return The object is updated with the new results.
     */
    compute(use_mito_default, mito_prefix, nmads) {
        this.changed = false;

        if (this.#inputs.changed || use_mito_default !== this.#parameters.use_mito_default || mito_prefix !== this.#parameters.mito_prefix) {
            utils.freeCache(this.#cache.metrics);

            if (this.valid()) {
                var mat = this.#inputs.fetchCountMatrix().get("RNA");
                var gene_info = this.#inputs.fetchFeatureAnnotations()["RNA"];

                // Finding the prefix.
                var subsets = utils.allocateCachedArray(mat.numberOfRows(), "Uint8Array", this.#cache, "metrics_buffer");
                subsets.fill(0);
                var sub_arr = subsets.array();

                // TODO: use the guessed features to narrow the Ensembl/symbol search.
                for (const key of gene_info.columnNames()) {
                    let val = gene_info.column(key);
                    if (use_mito_default) {
                        val.forEach((x, i) => {
                            if (mito.symbol.has(x) || mito.ensembl.has(x)) {
                                sub_arr[i] = 1;
                            }
                        });
                    } else {
                        var lower_mito = mito_prefix.toLowerCase();
                        val.forEach((x, i) => {
                            if(x.toLowerCase().startsWith(lower_mito)) {
                                sub_arr[i] = 1;
                            }
                        });
                    }
                }

                this.#cache.metrics = scran.perCellRnaQcMetrics(mat, [subsets]);
                this.changed = true;
            } else {
                delete this.#cache.metrics;
            }

            this.#parameters.use_mito_default = use_mito_default;
            this.#parameters.mito_prefix = mito_prefix;
        }

        if (this.changed || nmads !== this.#parameters.nmads) {
            utils.freeCache(this.#cache.filters);

            if (this.valid()) {
                let block = this.#inputs.fetchBlock();
                this.#cache.filters = scran.suggestRnaQcFilters(this.#cache.metrics, { numberOfMADs: nmads, block: block });
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
        let ghandle = handle.createGroup("quality_control");

        {
            let phandle = ghandle.createGroup("parameters"); 
            phandle.writeDataSet("use_mito_default", "Uint8", [], Number(this.#parameters.use_mito_default));
            phandle.writeDataSet("mito_prefix", "String", [], this.#parameters.mito_prefix);
            phandle.writeDataSet("nmads", "Float64", [], this.#parameters.nmads);
            phandle.writeDataSet("skip", "Uint8", [], 0); // back-compatibility only.
        }

        {
            let rhandle = ghandle.createGroup("results"); 

            if (this.valid()) {
                let mhandle = rhandle.createGroup("metrics");
                mhandle.writeDataSet("sums", "Float64", null, this.#cache.metrics.sums({ copy: "view" }));
                mhandle.writeDataSet("detected", "Int32", null, this.#cache.metrics.detected({ copy: "view" }));
                mhandle.writeDataSet("proportion", "Float64", null, this.#cache.metrics.subsetProportions(0, { copy: "view" }));

                let thandle = rhandle.createGroup("thresholds");
                thandle.writeDataSet("sums", "Float64", null, this.#cache.filters.thresholdsSums({ copy: "view" }));
                thandle.writeDataSet("detected", "Float64", null, this.#cache.filters.thresholdsDetected({ copy: "view" }));
                thandle.writeDataSet("proportion", "Float64", null, this.#cache.filters.thresholdsSubsetProportions(0, { copy: "view" }));

                rhandle.writeDataSet("discards", "Uint8", null, this.#cache.discard_buffer);
            }
        }
    }
}

/**************************
 ******** Loading *********
 **************************/

export function unserialize(handle, inputs) {
    let ghandle = handle.open("quality_control");

    let parameters = RnaQualityControlState.defaults(); 
    {
        let phandle = ghandle.open("parameters"); 
        parameters.use_mito_default = phandle.open("use_mito_default", { load: true }).values[0] > 0;
        parameters.mito_prefix = phandle.open("mito_prefix", { load: true }).values[0];
        parameters.nmads = phandle.open("nmads", { load: true }).values[0];
    }

    let output;
    let cache = {};
    try {
        let rhandle = ghandle.open("results");

        if ("metrics" in rhandle.children) { // QC metrics may not be reported if skipped.
            let mhandle = rhandle.open("metrics");
            let sums = mhandle.open("sums", { load: true }).values;

            cache.metrics = scran.emptyPerCellRnaQcMetricsResults(sums.length, 1);
            cache.metrics.sums({ fillable: true }).set(sums);

            let detected = mhandle.open("detected", { load: true }).values;
            cache.metrics.detected({ fillable: true }).set(detected);
            let proportions = mhandle.open("proportion", { load: true }).values;
            cache.metrics.subsetProportions(0, { fillable: true }).set(proportions);
        }

        if ("thresholds" in rhandle.children) { // if skip=true, QC thresholds may not be reported.
            let discards = rhandle.open("discards", { load: true }).values; 
            cache.discard_buffer = scran.createUint8WasmArray(discards.length);
            cache.discard_buffer.set(discards);

            let thandle = rhandle.open("thresholds");
            let thresholds_sums = thandle.open("sums", { load: true }).values;
            let thresholds_detected = thandle.open("detected", { load: true }).values;
            let thresholds_proportion = thandle.open("proportion", { load: true }).values;

            cache.filters = scran.emptySuggestRnaQcFiltersResults(1, thresholds_sums.length);
            cache.filters.thresholdsSums({ fillable: true }).set(thresholds_sums);
            cache.filters.thresholdsDetected({ fillable: true }).set(thresholds_detected);
            cache.filters.thresholdsSubsetProportions(0, { fillable: true }).set(thresholds_proportion);
        }

        output = new RnaQualityControlState(inputs, parameters, cache);
    } catch (e) {
        utils.freeCache(cache.metrics);
        utils.freeCache(cache.filters)
        utils.freeCache(output);
        throw e;
    }

    return output;
}
