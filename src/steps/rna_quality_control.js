import * as scran from "scran.js"; 
import * as bioc from "bioconductor";
import * as utils from "./utils/general.js";
import * as inputs_module from "./inputs.js";
import * as rutils from "../readers/index.js";

const baseUrl = "https://github.com/kanaverse/kana-special-features/releases/download/v1.0.0";

export const step_name = "rna_quality_control";

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
        utils.freeCache(this.#cache.keep_buffer);
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
        let output = { ...this.#parameters }; // avoid pass-by-reference links.
        output.species = bioc.CLONE(output.species);
        return output;
    }

    /**
     * @return {Uint8WasmArray} Buffer containing a vector of length equal to the number of cells,
     * where each element is truthy if the corresponding cell is to be retained after filtering.
     * This is available after running {@linkcode RnaQualityControlState#compute compute}.
     */
    fetchKeep() {
        return this.#cache.keep_buffer;
    }

    /**
     * @return {external:SuggestRnaQcFiltersResults} Result of filtering on the RNA-derived QC metrics.
     * This is available after running {@linkcode RnaQualityControlState#compute compute}.
     */
    fetchFilters() {
        return this.#cache.filters;
    }

    /**
     * @return {external:PerCellRnaQcMetricsResults} RNA-derived QC metrics,
     * available after running {@linkcode RnaQualityControlState#compute compute}.
     */
    fetchMetrics() {
        return this.#cache.metrics;
    }

    /****************************
     ******** Defaults **********
     ****************************/

    /**
     * @return {object} Object containing default parameters,
     * see the `parameters` argument in {@linkcode RnaQualityControlState#compute compute} for details.
     */
    static defaults () {
        return {
            guess_ids: true,
            gene_id_column: null,
            use_reference_mito: true,
            species: [],
            gene_id_type: "ENSEMBL",
            mito_prefix: "mt-",

            filter_strategy: "automatic", 
            nmads: 3,

            sum_threshold: 100,
            detected_threshold: 100,
            mito_threshold: 0.2
        };
    }

    static configureFeatureParameters(use_reference_mito, guesses) {
        let best_key = null;
        let best = { type: "symbol", species: "9606", confidence: 0 };

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
    static mitochondriaSpecies = [ 
        "9606",  // Mouse
        "10090", // Human
        "6239",  // C. elegans
        "10116", // Rat
        "9541",  // M. fascicularis
        "7227",  // Fly
        "7955",  // Zebrafish
        "9598"   // Chimp
    ];

    /***************************
     ******** Remotes **********
     ***************************/

    async #acquire_reference(species, feature_type) {
        let output = new Set;
        let mito_lists = RnaQualityControlState.#mito_lists;

        for (const s of species) {
            let target = s + "-mito-" + feature_type.toLowerCase() + ".txt.gz";
            if (!(target in mito_lists)) {
                let contents = await RnaQualityControlState.#downloadFun(baseUrl + "/" + target);
                let lines = await rutils.readLines2(contents, { compression: "gz" });
                mito_lists[target] = lines;
            }

            mito_lists[target].forEach(x => { output.add(x); });
        }

        return output;
    }

    static #mito_lists = {};

    /**
     * Flush all cached lists of mitochondrial genes.
     *
     * By default, {@linkcode RnaQualityControlState#compute compute} will cache the mitochondrial gene lists in a static member for re-use across {@linkplain RnaQualityControlState} instances.
     * These cached lists are not tied to any single instance and will not be removed by garbage collectors or by {@linkcode freeAnalysis}.
     * Rather, this function should be called to release the relevant memory.
     */
    static flush() {
        RnaQualityControlState.#mito_lists = {};
        return;
    }

    static #downloadFun = utils.defaultDownload;

    /**
     * Specify a function to download the reference mitochondrial gene lists.
     *
     * @param {function} fun - Function that accepts a single string containing a URL and returns any value that can be used in the {@linkplain SimpleFile} constructor.
     * This is most typically a Uint8Array of that URL's contents, but it can also be a path to a locally cached file on Node.js.
     *
     * @return `fun` is set as the global downloader for this step. 
     * The _previous_ value of the downloader is returned.
     */
    static setDownload(fun) {
        let previous = RnaQualityControlState.#downloadFun;
        RnaQualityControlState.#downloadFun = fun;
        return previous;
    }

    /***************************
     ******** Compute **********
     ***************************/

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     *
     * @param {object} parameters - Parameter object, equivalent to the `rna_quality_control` property of the `parameters` of {@linkcode runAnalysis}.
     * @param {boolean} [parameters.guess_ids] - Automatically choose feature-based parameters based on the feature annotation for the RNA modality.
     * If set to `true`, the following logic is applied:
     *
     * - If `use_reference_mito = true`, the annotation column that best matches human/mouse Ensembl/symbols is set as `gene_id_column`.
     *   Based on the identified species and feature type, `species` and `gene_id_type` are also set.
     * - If `use_reference_mito = false`, the annotation column that best matches human/mouse symbols is set as `gene_id_column`.
     *
     * @param {?(string|number)} [parameters.gene_id_column] - Name or index of the column of the feature annotations that contains the gene identifiers for the RNA modality.
     * If `null`, the row names are used.
     * Ignored if `guess_ids = true`.
     * @param {boolean} [parameters.use_reference_mito] - Whether to use the reference lists of mitochondrial genes.
     * If `false`, mitochondrial genes are instead identified from their prefix.
     * @param {Array} [parameters.species] - Array of strings specifying zero, one or more species to use to obtain a reference list of mitochondrial genes.
     * Each entry should be a taxonomy ID (e.g. `"9606"`, `"10090"`) as specified in {@linkcode RnaQualityControlState#mitochondriaSpecies mitochondriaSpecies}).
     * Ignored if `guess_ids = true`.
     * @param {string} [parameters.gene_id_type] - Name of the feature type in the reference list of mitochondrial genes.
     * This can be any one of `"ENSEMBL"`, `"SYMBOL"`, or `"ENTREZ"`.
     * Ignored if `guess_ids = true`.
     * @param {?string} [parameters.mito_prefix] - Case-insensitive prefix to use to identify mitochondrial genes from the dataset.
     * Only used when `use_reference_mito = false`; in such cases, `gene_id_column` should point to symbols.
     * If `null`, no prefix-based identification is performed.
     * @param {string} [parameters.filter_strategy] - Strategy for defining a filter threshold for the QC metrics.
     * This can be `"automatic"` or `"manual"`.
     * @param {number} [parameters.nmads] - Number of MADs to use for automatically selecting the filter threshold for each metric.
     * Only used when `filter_strategy = "automatic"`.
     * @param {number} [parameters.sum_threshold] - Manual threshold on the sum of counts for each cell.
     * Cells are only retained if their sums are equal to or greater than this threshold.
     * Only used when `filter_strategy = "manual"`.
     * @param {number} [parameters.detected_threshold] - Manual threshold on the detected number of features for each cell.
     * Cells are only retained if the detected number is equal to or greater than this threshold.
     * Only used when `filter_strategy = "manual"`.
     * @param {number} [parameters.mito_threshold] - Manual threshold on the mitochondrial proportion for each cell.
     * Cells are only retained if their totals are less than or equal to this threshold.
     * Only used when `filter_strategy = "manual"`.
     *
     * @return The object is updated with the new results.
     * @async
     */
    async compute(parameters) {
        parameters = utils.defaultizeParameters(parameters, RnaQualityControlState.defaults(), ["automatic"]);
        this.changed = false;

        // Some back-compatibility here.
        if (typeof parameters.guess_ids == "undefined") {
            if ("automatic" in parameters) {
                parameters.guess_ids = parameters.automatic;
            } else {
                parameters.guess_ids = true;
            }
        }

        if (
            this.#inputs.changed || 
            parameters.guess_ids !== this.#parameters.guess_ids ||
            parameters.use_reference_mito !== this.#parameters.use_reference_mito || 
            (
                !parameters.guess_ids && 
                (
                    parameters.gene_id_column !== this.#parameters.gene_id_column || 
                    (!parameters.use_reference_mito && parameters.mito_prefix !== this.#parameters.mito_prefix) ||
                    (
                        parameters.use_reference_mito && 
                        (
                            utils.changedParameters(parameters.species, this.#parameters.species) || 
                            parameters.gene_id_type !== this.#parameters.gene_id_type
                        )
                    )
                )
            ) 
        ) {
            utils.freeCache(this.#cache.metrics);

            if (this.valid()) {
                let gene_id_column2 = parameters.gene_id_column;
                let species2 = parameters.species;
                let gene_id_type2 = parameters.gene_id_type;

                if (parameters.guess_ids) {
                    let guesses = this.#inputs.guessRnaFeatureTypes();
                    let backcomp = RnaQualityControlState.configureFeatureParameters(parameters.use_reference_mito, guesses);
                    gene_id_column2 = backcomp.gene_id_column;
                    species2 = backcomp.species;
                    gene_id_type2 = backcomp.gene_id_type;
                }

                var gene_info = this.#inputs.fetchFeatureAnnotations()["RNA"];
                let val = (gene_id_column2 == null ? gene_info.rowNames() : gene_info.column(gene_id_column2));
                var subsets = utils.allocateCachedArray(gene_info.numberOfRows(), "Uint8Array", this.#cache, "metrics_buffer");
                subsets.fill(0);

                if (val !== null) {
                    if (parameters.use_reference_mito) {
                        let lists = await this.#acquire_reference(species2, gene_id_type2);
                        var sub_arr = subsets.array();
                        val.forEach((x, i) => {
                            if (lists.has(x)) {
                                sub_arr[i] = 1;
                            }
                        });
                    } else if (parameters.mito_prefix !== null) {
                        var lower_mito = parameters.mito_prefix.toLowerCase();
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

        if (this.changed || 
            parameters.filter_strategy !== this.#parameters.filter_strategy ||
            parameters.nmads !== this.#parameters.nmads ||
            parameters.sum_threshold !== this.#parameters.sum_threshold ||
            parameters.detected_threshold !== this.#parameters.detected_threshold ||
            parameters.mito_threshold !== this.#parameters.mito_threshold
        ) {
            utils.freeCache(this.#cache.filters);

            if (this.valid()) {
                let block = this.#inputs.fetchBlock();

                if (parameters.filter_strategy === "automatic") {
                    this.#cache.filters = scran.suggestRnaQcFilters(this.#cache.metrics, { numberOfMADs: parameters.nmads, block: block });
                } else if (parameters.filter_strategy === "manual") {
                    let block_levels = this.#inputs.fetchBlockLevels();
                    this.#cache.filters = scran.emptySuggestRnaQcFiltersResults(1, block_levels === null ? 1 : block_levels.length);
                    this.#cache.filters.sum({ copy: false }).fill(parameters.sum_threshold);
                    this.#cache.filters.detected({ copy: false }).fill(parameters.detected_threshold);
                    this.#cache.filters.subsetProportion(0, { copy: false }).fill(parameters.mito_threshold);
                } else {
                    throw new Error("unknown RNA QC filtering strategy '" + filter_strategy + "'");
                }

                var keep = utils.allocateCachedArray(this.#cache.metrics.numberOfCells(), "Uint8Array", this.#cache, "keep_buffer");
                this.#cache.filters.filter(this.#cache.metrics, { block: block, buffer: keep });
                this.changed = true;
            } else {
                delete this.#cache.filters;
            }
        }

        this.#parameters = parameters;
        this.#parameters.species = bioc.CLONE(this.#parameters.species); // avoid pass-by-reference behavior.
        return;
    }
}
