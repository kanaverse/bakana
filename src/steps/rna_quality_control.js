import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
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

const mito_lists = {};

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
    #automatic;

    constructor(inputs, parameters = null, cache = null) {
        if (!(inputs instanceof inputs_module.InputsState)) {
            throw new Error("'inputs' should be an InputsState object");
        }
        this.#inputs = inputs;

        this.#parameters = (parameters === null ? {} : parameters);
        this.#cache = (cache === null ? {} : cache);
        this.#automatic = false;
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

    /****************************
     ******** Defaults **********
     ****************************/

    static defaults () {
        return {
            automatic: true,
            gene_id_column: null,
            use_reference_mito: true,
            species: [],
            gene_id_type: "ENSEMBL",
            mito_prefix: "mt-",
            nmads: 3
        };
    }

    static configureFeatureParameters(use_reference_mito, guesses) {
        let best_key = null;
        let best = { type: "symbol", species: "human", confidence: 0 };

        if ("row_names" in guesses) {
            let val = guesses.row_names;
            if (val.confidence > best.confidence && (use_reference_mito || val.type == "symbol")) {
                best = val;
            }
        }

        for (const [key, val] of Object.entries(guesses.columns)) {
            if (val.confidence > best.confidence && (use_reference_mito || val.type == "symbol")) {
                best = val;
                best_key = key;
            }
        }

        return {
            gene_id_column: best_key,
            species: [best.species],
            gene_id_type: best.type.toUpperCase()
        };
    }

    /**
     * Array of strings containing the taxonomy IDs for species where mitochondrial gene lists are available.
     * @type {Array}
     */
    static mitochondriaSpecies = [ "9606", "10090", "6239" ];

    /***************************
     ******** Compute **********
     ***************************/

    async #acquire_reference(species, feature_type) {
        let output = new Set;

        for (const s0 of species) {
            let s = utils.toTaxonomy(s0);
            let target = s + "-mito-" + feature_type.toLowerCase() + ".txt.gz";
            if (!(target in mito_lists)) {
                let contents = await downloadFun(baseUrl + "/" + target);
                let lines = await rutils.readLines2(contents, { compression: "gz" });
                mito_lists[target] = lines;
            }

            mito_lists[target].forEach(x => { output.add(x); });
        }

        return output;
    }

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     *
     * @param {object} parameters - Parameter object, equivalent to the `rna_quality_control` property of the `parameters` of {@linkcode runAnalysis}.
     * @param {boolean} parameters.automatic - Automatically choose feature-based parameters based on the feature annotation for the RNA modality.
     * If set to `true`, the following logic is applied:
     *
     * - If `use_reference_mito = true`, the annotation column that best matches human/mouse Ensembl/symbols is set as `gene_id_column`.
     *   Based on the identified species and feature type, `species` and `gene_id_type` are also set.
     * - If `use_reference_mito = false`, the annotation column that best matches human/mouse symbols is set as `gene_id_column`.
     *
     * @param {?(string|number)} parameters.gene_id_column - Name or index of the column of the feature annotations that contains the gene identifiers for the RNA modality.
     * If `null`, the row names are used.
     * Ignored if `automatic = true`.
     * @param {boolean} parameters.use_reference_mito - Whether to use the reference lists of mitochondrial genes.
     * If `false`, mitochondrial genes are instead identified from their prefix.
     * @param {Array} parameters.species - Array of strings specifying zero, one or more species to use to obtain a reference list of mitochondrial genes.
     * Each entry can either be the common name (e.g., `"mouse"`, `"human"`) or a taxonomy ID (e.g. 9606, 10090; see {@linkcode RnaQualityControlState#mitochondriaSpecies mitochondriaSpecies}).
     * Ignored if `automatic = true`.
     * @param {string} parameters.gene_id_type - Name of the feature type in the reference list of mitochondrial genes.
     * This can be any one of `"ENSEMBL"`, `"SYMBOL"`, or `"ENTREZ"`.
     * Ignored if `automatic = true`.
     * @param {?string} parameters.mito_prefix - Case-insensitive prefix to use to identify mitochondrial genes from the dataset.
     * Only used when `use_reference_mito = false`; in such cases, `gene_id_column` should point to symbols.
     * If `null`, no prefix-based identification is performed.
     * @param {number} parameters.nmads - Number of MADs to use for automatically selecting the filter threshold for each metric.
     *
     * @return The object is updated with the new results.
     * @async
     */
    async compute(parameters) {
        let { mito_prefix, nmads } = parameters;
        let automatic;
        let use_reference_mito;
        let gene_id_column;
        let species;
        let gene_id_type;

        // Some back-compatibility here.
        if ("use_reference_mito" in parameters) {
            automatic = parameters.automatic;
            use_reference_mito = parameters.use_reference_mito;
            gene_id_column = parameters.gene_id_column;
            species = parameters.species;
            gene_id_type = parameters.gene_id_type;
        } else {
            automatic = true;
            use_reference_mito = parameters.use_mito_default;
            let def = RnaQualityControlState.defaults();
            gene_id_column = def.gene_id_column;
            species = def.species;
            gene_id_type = def.gene_id_type;
        }

        this.changed = false;

        if (
            this.#inputs.changed || 
            automatic !== this.#parameters.automatic ||
            use_reference_mito !== this.#parameters.use_reference_mito || 
            (
                !automatic && 
                (
                    gene_id_column !== this.#parameters.gene_id_column || 
                    (!use_reference_mito && mito_prefix !== this.#parameters.mito_prefix) ||
                    (
                        use_reference_mito && 
                        (
                            utils.changedParameters(species, this.#parameters.species) || 
                            gene_id_type !== this.#parameters.gene_id_type
                        )
                    )
                )
            ) 
        ) {
            utils.freeCache(this.#cache.metrics);

            if (this.valid()) {
                let gene_id_column2 = gene_id_column;
                let species2 = species;
                let gene_id_type2 = gene_id_type;

                if (automatic) {
                    let guesses = this.#inputs.guessRnaFeatureTypes();
                    let backcomp = RnaQualityControlState.configureFeatureParameters(use_reference_mito, guesses);
                    gene_id_column2 = backcomp.gene_id_column;
                    species2 = backcomp.species;
                    gene_id_type2 = backcomp.gene_id_type;
                }

                var gene_info = this.#inputs.fetchFeatureAnnotations()["RNA"];
                let val = (gene_id_column2 == null ? gene_info.rowNames() : gene_info.column(gene_id_column2));
                var subsets = utils.allocateCachedArray(gene_info.numberOfRows(), "Uint8Array", this.#cache, "metrics_buffer");
                subsets.fill(0);

                if (val !== null) {
                    if (use_reference_mito) {
                        let lists = await this.#acquire_reference(species2, gene_id_type2);
                        var sub_arr = subsets.array();
                        val.forEach((x, i) => {
                            if (lists.has(x)) {
                                sub_arr[i] = 1;
                            }
                        });
                    } else if (mito_prefix !== null) {
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
        }

        this.#parameters.automatic = automatic;
        this.#parameters.gene_id_column = gene_id_column;
        this.#parameters.use_reference_mito = use_reference_mito;
        this.#parameters.species = species;
        this.#parameters.gene_id_type = gene_id_type;
        this.#parameters.mito_prefix = mito_prefix;

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
export function setRnaQualityControlDownload(fun) {
    let previous = downloadFun;
    downloadFun = fun;
    return previous;
}
