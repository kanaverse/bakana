import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import { mito } from "./mito.js";
import * as inputs_module from "./inputs.js";
import * as rutils from "../readers/index.js";

var downloadFun = async (url) => {
    let resp = await fetch(url);
    if (!resp.ok) {
        throw new Error("failed to fetch content at " + url + "(" + resp.status + ")");
    }
    return new Uint8Array(await resp.arrayBuffer());
};

const baseUrl = "https://github.com/kanaverse/kana-special-features/releases/download/v1.0.0";

export const step_name = "rna_quality_control";

let mito_lists = {};
const species_to_taxa = { mouse: 10090, human: 9606 };

/**
 * Results of computing per-cell RNA-derived QC metrics,
 * see [here](https://kanaverse.github.io/scran.js/PerCellRnaQcMetricsResults.html) for details.
 *
 * @external PerCellRnaQcMetricsResults
 */

/**
 * Suggested filters for the RNA-derived QC metrics,
 * see [here](https://kanaverse.github.io/scran.js/SuggestRnaQcFiltersResults.html) for details.
 *
 * @external SuggestRnaQcFiltersResults
 */

/**
 * This step applies quality control on the RNA count matrix.
 * Specifically, it computes the QC metrics and filtering thresholds, 
 * wrapping the [`perCellRnaQcMetrics`](https://kanaverse.github.io/scran.js/global.html#perCellRnaQcMetrics)
 * and [`suggestRnaQcFilters`](https://kanaverse.github.io/scran.js/global.html#suggestRnaQcFilters) functions
 * from [**scran.js**](https://github.com/kanaverse/scran.js).
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
            dataset_id_column: null
            reference_mito_species: null,
            reference_mito_features: "ENSEMBL",
            mito_prefix: "mt-",
            nmads: 3
        };
    }

    async #acquire_reference(species, feature_type) {
        if (typeof species == "string") {
            if (species in species_to_taxa) {
                species = species_to_taxa[species];
            } else {
                throw new Error("unknown species '" + species + "' for extracting mitochondrial lists");
            }
        }

        let target = String(species) + "-mito-" + feature_type.toLowerCase() + ".txt.gz";
        if (!(target in mito_lists)) {
            let contents = await downloadFun(baseUrl + target);
            let lines = rutils.readLines2(contents, { compression: "gz" });
            mito_lists[target] = new Set(lines.map(x => x[0]));
        }

        return mito_lists[target];
    }

    #back_compatibility_hack(use_mito_default) {
        let genes = this.#inputs.fetchFeatureAnnotations()["RNA"];
        let best_key = null;
        let best = { type: "symbol", species: "human", confidence: 0 };

        if (genes.rowNames() !== null) {
            let col = genes.rowNames();
            let val = scran.guessFeatures(col);
            if (val.confidence > best.confidence && (use_mito_default || val.type == "symbol")) {
                best = val;
            }
        }

        for (const key of genes.columnNames()) {
            let col = genes.column(key);
            let val = scran.guessFeatures(col);
            if (val.confidence > best.confidence && (use_mito_default || val.type == "symbol")) {
                best = val;
                best_key = key;
            }
        }

        return {
            dataset_id_column: best_key,
            reference_mito_species: best.species,
            reference_mito_features: best.type.toUpperCase()
        };
    }

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     *
     * @param {object} parameters - Parameter object, equivalent to the `rna_quality_control` property of the `parameters` of {@linkcode runAnalysis}.
     * @param {?(string|number)} parameters.dataset_id_column - Name or index of the column of the feature annotations that contains the feature identifiers for the RNA modality.
     * If `null`, the row names are used.
     * @param {?(string|number)} parameters.reference_mito_species - Species for which to look up the mitochondrial features.
     * This can either be the common name (e.g., `"mouse"`, `"human"`) or a taxonomy ID (e.g. 9606, 10090).
     * If `null`, no look-up is performed.
     * @param {string} parameters.reference_mito_features - Name of the feature type in the reference list of mitochondrial features.
     * This can be `"ENSEMBL"`, `"SYMBOL"`, or `"ENTREZ"`.
     * @param {string} parameters.mito_prefix - Case-insensitive prefix to use to identify mitochondrial genes from the dataset.
     * Only used when `reference_mito_species = null`; in such cases, `dataset_id_column` should point to symbols.
     * @param {number} parameters.nmads - Number of MADs to use for automatically selecting the filter threshold for each metric.
     *
     * @return The object is updated with the new results.
     * @async
     */
    async compute(parameters) {
        let { mito_prefix, nmads } = parameters;
        let dataset_id_column;
        let reference_mito_species;
        let reference_mito_features;

        // Some back-compatibility here.
        if ("dataset_id_column" in parameters) {
            dataset_id_column = parameters.dataset_id_column;
            reference_mito_species = parameters.reference_mito_species;
            reference_mito_features = parameters.reference_mito_features;
        } else {
            let backcomp = this.#back_compatibility_hack(parameters.use_mito_default);
            dataset_id_column = backcomp.dataset_id_column;
            reference_mito_species = backcomp.reference_mito_species;
            reference_mito_features = backcomp.reference_mito_features;
        }

        this.changed = false;

        if (this.#inputs.changed || use_mito_default !== this.#parameters.use_mito_default || mito_prefix !== this.#parameters.mito_prefix) {
            utils.freeCache(this.#cache.metrics);

            if (this.valid()) {
                var gene_info = this.#inputs.fetchFeatureAnnotations()["RNA"];
                let val = (dataset_id_column == null ? gene_info.rowNames() : gene_info.column(dataset_id_column));
                var subsets = utils.allocateCachedArray(gene_info.numberOfRows(), "Uint8Array", this.#cache, "metrics_buffer");
                subsets.fill(0);

                if (val !== null) {
                    if (reference_mito_species !== null) {
                        let lists = await this.#acquire_reference(reference_mito_species, reference_mito_features);
                        var sub_arr = subsets.array();
                        val.forEach((x, i) => {
                            if (lists.has(x)) {
                                sub_arr[i] = 1;
                            }
                        });
                    } else {
                        var lower_mito = mito_prefix.toLowerCase();
                        var sub_arr = subsets.array();
                        val.forEach((x, i) => {
                            if(x.toLowerCase().startsWith(lower_mito)) {
                                sub_arr[i] = 1;
                            }
                        });
                    }
                }

                var mat = this.#inputs.fetchCountMatrix().get("RNA");
                this.#cache.metrics = scran.perCellRnaQcMetrics(mat, [subsets]);
                this.changed = true;
            } else {
                delete this.#cache.metrics;
            }

            this.#parameters.dataset_id_column = dataset_id_column;
            this.#parameters.reference_mito_species = reference_mito_species;
            this.#parameters.reference_mito_features = reference_mito_features;
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
}

/**************************
 ******** Loading *********
 **************************/

export function unserialize(handle, inputs) {
    let ghandle = handle.open("rna_quality_control" in handle.children ? "rna_quality_control" : "quality_control");

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

/**************************
 ******** Setters *********
 **************************/

/**
 * Specify a function to download the reference mitochondrial gene lists.
 *
 * @param {function} fun - Function that accepts a single string containing a URL and returns any value that can be used in the {@linkplain SimpleFile} constructor.
 * This is most typically a Uint8Array of that URL's contents, but it can also be a path to a locally cached file on Node.js.
 *
 * @return `fun` is set as the global downloader for this step. 
 * The _previous_ value of the downloader is returned.
 */
export function setRnaQualityControlDownloadFun(fun) {
    let previous = downloadFun;
    downloadFun = fun;
    return previous;
}
