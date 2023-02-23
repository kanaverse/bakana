import * as scran from "scran.js";
import * as bioc from "bioconductor";
import * as utils from "./utils/general.js";
import * as mutils from "./utils/markers.js";
import * as rutils from "../readers/index.js";
import * as inputs_module from "./inputs.js";
import * as filter_module from "./cell_filtering.js";
import * as norm_module from "./rna_normalization.js";
import * as markers_module from "./marker_detection.js";

var downloadFun = async (url) => {
    let resp = await fetch(url);
    if (!resp.ok) {
        throw new Error("failed to fetch content at " + url + "(" + resp.status + ")");
    }
    return new Uint8Array(await resp.arrayBuffer());
};

const base = "https://github.com/LTLA/kana-feature-sets/releases/download/v1.0.0";

export const step_name = "feature_set_enrichment";

/**
 * This step tests for enrichment of particular feature sets in the set of top marker genes,
 * based on marker rankings from {@linkplain MarkerDetectionState}.
 * It wraps the [`testFeatureSetEnrichment`](https://kanaverse.github.io/scran.js/global.html#testFeatureSetEnrichment) 
 * and [`scoreFeatureSet`](https://kanaverse.github.io/scran.js/global.html#scoreFeatureSet) functions
 * from [**scran.js**](https://github.com/kanaverse/scran.js).
 * 
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class FeatureSetEnrichmentState {
    #inputs;
    #filter;
    #normalized;
    #markers;
    #parameters;
    #cache;

    constructor(inputs, filter, normalized, markers, parameters = null, cache = null) {
        if (!(inputs instanceof inputs_module.InputsState)) {
            throw new Error("'inputs' should be a State object from './inputs.js'");
        }
        this.#inputs = inputs;

        if (!(filter instanceof filter_module.CellFilteringState)) {
            throw new Error("'filter' should be a CellFilteringState object");
        }
        this.#filter = filter;

        if (!(normalized instanceof norm_module.RnaNormalizationState)) {
            throw new Error("'normalized' should be a RnaNormalizationState object from './rna_normalization.js'");
        }
        this.#normalized = normalized;

        if (!(markers instanceof markers_module.MarkerDetectionState)) {
            throw new Error("'markers' should be a State object from './marker_detection.js'");
        }
        this.#markers = markers;

        this.#parameters = (parameters === null ? {} : parameters);
        this.#cache = (cache === null ? {} : cache);
        this.#cache.adhoc_results = {};

        this.changed = false;
    }

    free() {
        this.#cache = {};
        return; // nothing extra to free here.
    }

    /***************************
     ******** Getters **********
     ***************************/

    valid() {
        let mat = this.#inputs.fetchCountMatrix();
        return mat.has("RNA");
    }

    /**
     * @return {object} Object where each entry corresponds to a feature set collection.
     * Each value is itself an object containing:
     *
     * - `names`: Array of strings containing the names of the feature sets in the collection.
     * - `descriptions`: Array of strings containing the set descriptions.
     * - `sizes`: Int32Array containing the set sizes.
     * - `universe`: number of features in the universe.
     */
    fetchCollectionDetails() {
        let collected = this.#cache.prepared;
        let output = {};
        for (const [k, v] of Object.entries(collected)) {
            output[k] = { names: v.names, descriptions: v.descriptions, sizes: v.sizes, universe: v.indices.length };
        }
        return output;
    }

    /**
     * @param {number} group - Cluster index of interest.
     * @param {string} effect_size - Effect size to use for ranking.
     * This should be one of `"cohen"`, `"auc"`, `"lfc"` or `"delta_detected"`.
     * @param {string} summary - Summary statistic to use for ranking.
     * This should be one of `"min"`, `"mean"` or `"min_rank"`.
     *
     * @return {object} Object where each entry corresponds to a feature set collection.
     * Each value is itself an object containing:
     *
     * - `counts`: Int32Array containing the number of markers present in each set.
     * - `pvalues`: Float64Array containing the enrichment p-values for each set.
     * - `num_markers`: number of markers selected for testing.
     */
    fetchGroupResults(group, effect_size, summary) {
        let key = String(group) + ":" + effect_size + ":" + summary;

        if (!(key in this.#cache.adhoc_results)) {
            let res = this.#markers.fetchResults()["RNA"];
            let top_markers = this.#parameters.top_markers;
            let collections = this.#cache.prepared;

            if (effect_size == "delta_detected") {
                effect_size = "deltaDetected";
            }

            // Avoid picking down-regulated genes in the marker set.
            let min_threshold = effect_size == "auc" ? 0.5 : 0;

            // Larger is better except for 'min_rank'.
            let use_largest = effect_size !== "min_rank"; 
            let sumidx = mutils.summaries2int[summary];

            let output = {};
            for (const [name, info] of Object.entries(collections)) {
                let stats = res[effect_size](group, { summary: sumidx, copy: false });
                let curstats = bioc.SLICE(stats, info.indices);
                let threshold = scran.computeTopThreshold(curstats, top_markers, { largest: use_largest });
                let in_set = [];

                if (use_largest) {
                    if (threshold < min_threshold) {
                        threshold = min_threshold;
                    }
                    curstats.forEach((x, i) => {
                        if (x >= threshold) {
                            in_set.push(i);
                        }
                    });
                } else {
                    curstats.forEach((x, i) => {
                        if (x <= threshold) {
                            in_set.push(i);
                        }
                    });
                }

                let computed = scran.testFeatureSetEnrichment(in_set, info.sets, curstats.length);
                output[name] = { 
                    counts: computed.count, 
                    pvalues: computed.pvalue,
                    num_markers: in_set.length
                };
            }

            this.#cache.adhoc_results[key] = output;
        }

        return this.#cache.adhoc_results[key];
    }

    /**
     * @param {string} collection - Name of the collection.
     * @param {number} set_index - Index of a feature set inside the specified collection.
     *
     * @return {Int32Array} Array containing the row indices of the RNA count matrix corresponding to the genes in the specified set.
     */
    fetchFeatureSetIndices(collection, set_index) {
        let current = this.#cache.prepared[collection];
        let set = current.sets[set_index];
        let indices = current.indices;
        return bioc.SLICE(indices, set);
    }

    /**
     * @param {string} collection - Name of the collection.
     * @param {number} set_index - Index of a feature set inside the specified collection.
     *
     * @return {Object} Object containing:
     *
     * - `indices`: Int32Array containing the row indices of the genes in the set, relative to the RNA count matrix.
     * - `weights`: Float64Array containing the weights of each gene in the set.
     * - `scores`: Float64Array containing the feature set score for each cell.
     */
    fetchPerCellScores(collection, set_index) {
        let indices = this.fetchFeatureSetIndices(collection, set_index);
        // console.log(bioc.SLICE(this.#inputs.fetchFeatureAnnotations().RNA.column("id"), indices));

        let mat = this.#normalized.fetchNormalizedMatrix();
        let features = utils.allocateCachedArray(mat.numberOfRows(), "Uint8Array", this.#cache, "set_buffer");
        features.fill(0);
        let farr = features.array();
        indices.forEach(x => { farr[x] = 1; }); 

        return scran.scoreFeatureSet(mat, features, { block: this.#filter.fetchFilteredBlock() });
    }

    /**
     * @return {object} Object containing the parameters.
     */
    fetchParameters() {
        return { ...this.#parameters };
    }

    /****************************
     ******** Defaults **********
     ****************************/

    static defaults() {
        return {
            collections: [],
            automatic: true,
            gene_id_column: null, 
            gene_id_type: "ENSEMBL", 
            top_markers: 100
        };
    }

    static configureFeatureParameters(guesses) {
        let best_key = null;
        let best = { type: "symbol", species: "human", confidence: 0 };

        if ("row_names" in guesses) {
            let val = guesses.row_names;
            if (val.confidence > best.confidence) {
                best = val;
            }
        }

        for (const [key, val] of Object.entries(guesses.columns)) {
            if (val.confidence > best.confidence) {
                best = val;
                best_key = key;
            }
        }

        return {
            gene_id_column: best_key,
            gene_id_type: best.type.toUpperCase(),
            species: [best.species]
        };
    }

    /**
     * Available feature set collections for each species.
     * Each key is a taxonomy ID and each value is an array of names of feature set collections for that species.
     * @type {object}
     */
    static availableCollections = {
        "10090": [ "mouse-GO" ],
        "9606": [ "human-GO" ],
        "6239": [ "worm-GO" ],
        "10116": [ "rat-GO" ],
        "7227": [ "fly-GO" ],
        "7955": [ "zebrafish-GO" ],
        "9598": [ "chimp-GO" ]
    };

    /***************************
     ******** Remotes **********
     ***************************/

    async #load_collection(name) {
        let loaded = FeatureSetEnrichmentState.#all_loaded;
        if (!(name in loaded)) {
            let suffixes = [
                "features.csv.gz",
                "sets.txt.gz"
            ];

            let contents = await Promise.all(
                suffixes.map(
                    async suffix => {
                        let full = name + "_" + suffix;
                        let b = await downloadFun(base + "/" + full);
                        return new rutils.SimpleFile(b, { name: full })
                    }
                )
            );

            let gene_info = {};
            {
                let genes = await rutils.readTable2(contents[0].content(), { delim: "," });
                let headers = genes.shift();
                for (const x of headers) {
                    gene_info[x] = [];
                }
                for (const line of genes) {
                    for (var i = 0; i < headers.length; i++) {
                        let curfield = line[i];
                        gene_info[headers[i]].push(curfield == "" ? null : curfield);
                    }
                }
            }

            let set_members = [];
            let set_name = [];
            let set_description = [];
            {
                let features = await rutils.readLines2(contents[1].content());
                for (const line of features) {
                    let values = line.split("\t");
                    set_name.push(values[0]);
                    set_description.push(values[1]);

                    let last = Number(values[2]);
                    let members = [last];
                    for (var i = 3; i < values.length; i++) {
                        let latest = Number(values[i]) + last;
                        members.push(latest);
                        last = latest;
                    }
                    set_members.push(members);
                }
            }

            loaded[name] = {
                features: gene_info,
                members: set_members,
                names: set_name,
                descriptions: set_description
            };
        }

        return loaded[name];
    }

    static #all_loaded = {};

    /**
     * Flush all cached feature set collections.
     *
     * By default, {@linkcode FeatureSetEnrichmentState#compute compute} will cache the feature set collections in a static member for re-use across {@linkplain FeatureSetEnrichmentState} instances.
     * These cached collections are not tied to any single instance and will not be removed by garbage collectors or by {@linkcode freeAnalysis}.
     * Rather, this function should be called to release the relevant memory.
     */
    static flush() {
        FeatureSetEnrichmentState.#all_loaded = {};
        return;
    }

    /***************************
     ******** Compute **********
     ***************************/

    #remap_collection(name, data_id, ref_id) {
        let feats = this.#inputs.fetchFeatureAnnotations()["RNA"];

        let data_id_col;
        if (data_id == null) {
            data_id_col = feats.rowNames();
            if (data_id_col == null) {
                throw new Error("no row names available in the feature annotations");
            }
        } else {
            data_id_col = feats.column(data_id);
        }

        let loaded = FeatureSetEnrichmentState.#all_loaded[name];
        if (!(ref_id in loaded.features)){
            throw new Error("no column '" + ref_id + "' in the feature set annotations");
        }

        let output = scran.remapFeatureSets(data_id_col, loaded.features[ref_id], loaded.members);
        let sizes = new Int32Array(output.sets.length);
        output.sets.forEach((x, i) => { sizes[i] = x.length });
        output.sizes = sizes;

        return output;
    }

    async #prepare_collections(collections, species, gene_id_column, gene_id_type) {
        let allowable = new Set;
        for (const s of species) {
            if (s in FeatureSetEnrichmentState.availableCollections) {
                FeatureSetEnrichmentState.availableCollections[s].forEach(x => allowable.add(x));
            }
        }

        let collected = {};
        for (const x of collections) {
            if (!allowable.has(x)) {
                continue;
            }

            let loaded = await this.#load_collection(x);
            let mapped = this.#remap_collection(x, gene_id_column, gene_id_type);

            collected[x] = {
                names: loaded.names,
                descriptions: loaded.descriptions,
                indices: mapped.target_indices,
                sets: mapped.sets,
                sizes: mapped.sizes
            };
        }

        this.#cache.prepared = collected;
        return;
    }

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     *
     * @param {object} parameters - Parameter object, equivalent to the `feature_set_enrichment` property of the `parameters` of {@linkcode runAnalysis}.
     * @param {Array} parameters.collections - Array of strings containing the names of collections to be tested, see {@linkcode FeatureSetEnrichmentState#availableCollections availableCollections}.
     * @param {boolean} parameters.automatic - Automatically choose feature-based parameters based on the feature annotation for the RNA modality.
     * If `true`, the column of the annotation that best matches human/mouse Ensembl/symbols is identified and used to set `species`, `gene_id_column`, `gene_id_type`.
     * @param {Array} parameters.species - Array of strings specifying zero, one or more species involved in this dataset.
     * Each entry should be a taxonomy ID (e.g. `"9606"`, `"10090"`) as specified in {@linkcode FeatureSetEnrichmentState#availableCollections availableCollections}.
     * This is used internally to filter `collections` to the entries relevant to these species. 
     * Ignored if `automatic = true`.
     * @param {?(string|number)} parameters.gene_id_column - Name or index of the column of the RNA entry of {@linkcode InputsState#fetchFeatureAnnotations InputsState.fetchFeatureAnnotations} containing the identity of each gene. 
     * If `null`, identifiers are taken from the row names.
     * Ignored if `automatic = true`.
     * @param {string} parameters.gene_id_type - Type of feature identifier in `gene_id_column`.
     * This should be one of `"ENSEMBL"`, `"SYMBOL"` or `"ENTREZ"`
     * Ignored if `automatic = true`.
     * @param {number} parameters.top_markers - Number of top markers to use when testing for enrichment.
     *
     * @return The state is updated with new results.
     */
    async compute(parameters) {
        let { collections, automatic, species, gene_id_column, gene_id_type, top_markers } = parameters;
        this.changed = false;
        if (this.#inputs.changed) {
            this.#cache.mapped = {};
            this.changed = true;
        }

        if (this.valid()) {
            if (
                utils.changedParameters(this.#parameters.collections, collections) ||
                automatic !== this.#parameters.automatic ||
                (
                    !automatic && 
                    (
                        this.#parameters.gene_id_column !== gene_id_column || 
                        this.#parameters.gene_id_type !== gene_id_type ||
                        utils.changedParameters(this.#parameters.species, species)
                    )
                )
            ) {
                let gene_id_column2 = gene_id_column;
                let gene_id_type2 = gene_id_type;
                let species2 = species;

                if (automatic) {
                    let guesses = this.#inputs.guessRnaFeatureTypes();
                    let auto = FeatureSetEnrichmentState.configureFeatureParameters(guesses);
                    gene_id_column2 = auto.gene_id_column;
                    gene_id_type2 = auto.gene_id_type;
                    species2 = auto.species;
                }

                await this.#prepare_collections(collections, species2, gene_id_column2, gene_id_type2);
                this.changed = true;
            }

            if (this.changed || this.#markers.changed || top_markers !== this.#parameters.top_markers) {
                this.#cache.adhoc_results = {};
                this.changed = true;
            }

        }

        this.#parameters.automatic = automatic;
        this.#parameters.species = species;
        this.#parameters.collections = collections;
        this.#parameters.gene_id_column = gene_id_column;
        this.#parameters.gene_id_type = gene_id_type;
        this.#parameters.top_markers = top_markers;
        return;
    }
}

/**************************
 ******** Loading *********
 **************************/

export function unserialize(handle, inputs, filter, normalized, markers) {
    let parameters = {};
    let cache = {};

    // Protect against old analysis states that don't have cell_labelling.
    if ("feature_set_enrichment" in handle.children) {
        let ghandle = handle.open("feature_set_enrichment");
        
        {
            let phandle = ghandle.open("parameters");
            parameters.collections = phandle.open("collections", { load: true }).values;
            for (const k of [ "gene_id_column", "gene_id_type", "top_markers" ]) {
                parameters[k] = phandle.open(k, { load: true }).values[0];
            }
        }
    }

    return new FeatureSetEnrichmentState(inputs, filter, normalized, markers, parameters, cache);
}

/**************************
 ******** Setters *********
 **************************/

/**
 * Specify a function to download feature sets.
 *
 * @param {function} fun - Function that accepts a single string containing a URL and returns any value that can be used in the {@linkplain SimpleFile} constructor.
 * This is most typically a Uint8Array of that URL's contents, but it can also be a path to a locally cached file on Node.js.
 *
 * @return `fun` is set as the global downloader for this step.
 * The _previous_ value of the downloader is returned.
 */
export function setFeatureSetEnrichmentDownload(fun) {
    let previous = downloadFun;
    downloadFun = fun;
    return previous;
}
