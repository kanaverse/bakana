import * as scran from "scran.js";
import * as bioc from "bioconductor";
import * as gesel from "gesel";

import * as utils from "./utils/general.js";
import * as mutils from "./utils/markers.js";
import * as rutils from "../readers/index.js";
import * as inputs_module from "./inputs.js";
import * as filter_module from "./cell_filtering.js";
import * as norm_module from "./rna_normalization.js";
import * as markers_module from "./marker_detection.js";

export const step_name = "feature_set_enrichment";

/****************************
 ******** Internals *********
 ****************************/

class FeatureSetManager {
    #cache;

    constructor() {
        this.#cache = {};
    }

    free() {
        utils.freeCache(this.#cache.set_buffer);
        this.#cache = {};
    }

    static flush() {
        // TODO: call a gesel flush() function.
        return;
    }

    static setDownload(fun) {
        console.warn("'FeatureSetState.setDownload' is a no-op, uses 'gesel.setReferenceDownload' and 'gesel.setGeneDownload' instead");
        return;
    }

    async prepare(feats, species, gene_id_column, gene_id_type) {
        let gene_offset = 0;
        let set_offset = 0;
        let collection_offset = 0;

        // Finding our genes.
        let data_id_col;
        if (data_id == null) {
            data_id_col = feats.rowNames();
            if (data_id_col == null) {
                throw new Error("no row names available in the feature annotations");
            }
        } else {
            data_id_col = feats.column(data_id);
        }

        let gene_universe = [];
        let set_names = [];
        let set_descriptions = [];
        let set_sets = [];
        let set_indices = [];
        let set_collections = [];
        let remapped = [];
        let collection_names = [];
        let collection_descriptions = [];
        let collection_species = [];

        for (const spec of species) {
            let gene_mapping = gesel.searchGenes(spec, data_id_col, { types: [ gene_id_type.toLowerCase() });
            let reverse_mapping = new Map;
            for (var i = 0; i < gene_mapping.length; i++) {
                let current = gene_mapping[i];
                if (current.length == 0) {
                    continue;
                }

                let gene = i + gene_offset;
                universe.push(gene);
                for (const gesel_gene of current) {
                    let found = reverse_mapping.get(gesel_gene);
                    if (typeof found == "undefined") {
                        found = [];
                        reverse_mapping.set(gesel_gene, found);
                    }
                    found.push(gene);
                }
            }

            // Formatting the details for each set. This includes reindexing
            // the gesel gene indices to refer to indices in the input features.
            let all_sets = await gesel.fetchAllSets(spec);
            let all_sets2genes = await gesel.fetchGenesForAllSets(spec);

            for (var i = 0; i < N; i++) {
                set_names.push(all_sets[i].name);
                set_names.push(all_sets[i].description);
                set_names.push(all_sets[i].collection + collection_offset);

                let subset = [];
                for (const gesel_gene of all_sets2genes[i]) {
                    let found = reverse_mapping.get(gesel_gene);
                    if (typeof found !== "undefined") {
                        for (const common_index of found) {
                            subset.push(gene);
                        }
                    }
                }

                set_indices.push(new Int32Array(subset));
                set_sizes.push(subset.length);
            }

            // Reindexing the reverse mapping to get indices to the input features.
            let all_genes2sets = await gesel.fetchSetsForAllGenes(spec);
            for (var i = 0; i < gene_mapping.length; i++) {
                let current = new Set;
                for (const gesel_gene of gene_mapping[i]) {
                    for (const set of all_genes2sets[gesel_gene]) {
                        current.add(set);
                    }
                }
                remapped.push((new Int32Array(current)).sort());
            }

            // Sticking the collection details somewhere.
            let all_collections = await gesel.fetchAllCollections(spec);
            for (var i = 0; i < all_collections.length; i++) {
                collection_names.push(all_collections[i].name);
                collection_descriptions.push(all_collections[i].description);
                collection_species.push(spec);
            }
        }

        this.#cache.universe = new Int32Array(gene_universe);

        this.#cache.sets = {
            names: set_names,
            descriptions: set_descriptions,
            sets: set_indices,
            sizes: new Int32Array(set_sizes),
            collections: new Int32Array(set_collections)
        };

        this.#cache.collections = {
            name: collection_names,
            descriptions: collection_descriptions,
            species: collection_species
        };

        this.#cache.mapping_to_sets = remapped;

        return;
    }

    fetchCollectionDetails() {
        return this.#cache.collections;
    }

    fetchSetDetails() {
        return { 
            names: this.#cache.sets.names,
            descriptions: this.#cache.sets.descriptions,
            sizes: this.#cache.sets.sizes,
            collections: this.#cache.sets.collections
        };
    }

    fetchUniverseSize() {
        return this.#cache.universe.length;
    }

    computeEnrichment(group, effect_size, summary, markers, top_markers) {
        if (effect_size == "delta_detected") {
            effect_size = "deltaDetected";
        }

        // Avoid picking down-regulated genes in the marker set.
        let min_threshold = effect_size == "auc" ? 0.5 : 0;

        // Larger is better except for 'min_rank'.
        let use_largest = effect_size !== "min_rank"; 
        let sumidx = mutils.summaries2int[summary];

        let stats = markers[effect_size](group, { summary: sumidx, copy: false });
        let curstats = bioc.SLICE(stats, this.#cache.universe);
        let threshold = scran.computeTopThreshold(curstats, top_markers, { largest: use_largest });

        let set_counts = new Map;
        let num_top = 0;
        let add = sub_index => {
            let gene = this.#cache.universe[sub_index];
            for (const set of this.#cache.mapping_to_sets[gene]) {
                let found = set_counts.get(set);
                if (typeof found == "undefined") {
                    found = 0;
                }
                set_counts.set(set, found + 1);
            }
            num_top++;
        };

        if (use_largest) {
            if (threshold < min_threshold) {
                threshold = min_threshold;
            }
            curstats.forEach((x, i) => {
                if (x >= threshold) {
                    add(i);
                }
            });
        } else {
            curstats.forEach((x, i) => {
                if (x <= threshold) {
                    add(i);
                }
            });
        }

        let indices = new Int32Array(set_counts.size);
        let set_ids = new Int32Array(set_counts.size);
        let counts = new Int32Array(set_counts.size);
        let pvalues = new Float64Array(set_counts.size);
        let counter = 0;
        for (const [k, v] of set_counts) {
            indices[counter] = counter;
            set_ids[counter] = k;
            counts[counter] = v;
            pvalues[counter] = gesel.testEnrichment(v, num_top, this.#cache.sets[k].size, this.#cache.universe.length);
        }

        // Sorting by p-value.
        indices.sort((a, b) => pvalues[a] - pvalues[b]);
        return {
            set_ids: bioc.SLICE(set_ids, indices),
            counts: bioc.SLICE(counts, indices),
            pvalues: bioc.SLICE(pvalues, indices),
            num_markers: num_top
        };
    }

    fetchFeatureSetIndices(set_id) {
        return this.#cache.sets.sets[set_id];
    }

    computePerCellScores(set_id, normalized, block) {
        let indices = this.fetchFeatureSetIndices(set_id);
        // console.log(bioc.SLICE(this.#inputs.fetchFeatureAnnotations().RNA.column("id"), indices));

        let features = utils.allocateCachedArray(normalized.numberOfRows(), "Uint8Array", this.#cache, "set_buffer");
        features.fill(0);
        let farr = features.array();
        indices.forEach(x => { farr[x] = 1; }); 

        return scran.scoreFeatureSet(normalized, features, { block: block });
    }
}

// More internal functions, mostly related to wrangling with parameters.

function _configureFeatureParameters(guesses) {
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

async function _buildCollections(old_parameters, manager, collections, automatic, species, gene_id_column, gene_id_type, annofun, guessfun) {
    if (
        utils.changedParameters(old_parameters.collections, collections) ||
        automatic !== old_parameters.automatic ||
        (
            !automatic && 
            (
                old_parameters.gene_id_column !== gene_id_column || 
                old_parameters.gene_id_type !== gene_id_type ||
                utils.changedParameters(old_parameters.species, species)
            )
        )
    ) {
        let gene_id_column2 = gene_id_column;
        let gene_id_type2 = gene_id_type;
        let species2 = species;

        if (automatic) {
            let auto = _configureFeatureParameters(guessfun());
            gene_id_column2 = auto.gene_id_column;
            gene_id_type2 = auto.gene_id_type;
            species2 = auto.species;
        }

        await manager.prepare(annofun(), collections, species2, gene_id_column2, gene_id_type2);
        return true;
    }

    return false;
}

function _transplantParameters(parameters, collections, automatic, species, gene_id_column, gene_id_type, top_markers) {
    parameters.automatic = automatic;
    parameters.species = bioc.CLONE(species); // make a copy to avoid pass-by-ref behavior.
    parameters.collections = bioc.CLONE(collections);
    parameters.gene_id_column = gene_id_column;
    parameters.gene_id_type = gene_id_type;
    parameters.top_markers = top_markers;
}

function _fetchParameters(parameters) {
    // Avoid pass-by-reference behavior.
    let out = { ...parameters };
    out.species = bioc.CLONE(out.species);
    out.collections = bioc.CLONE(out.collections);
    return out;
}

/************************
 ******** State *********
 ************************/

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

    #parameters;
    #manager;

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

        this.#parameters = (parameters === null ? {} : parameters);
        this.#manager = new FeatureSetManager;
        this.changed = false;
    }

    /**
     * Frees all resources associated with this instance.
     */
    free() {
        this.#manager.free();
        return; 
    }

    valid() {
        let mat = this.#inputs.fetchCountMatrix();
        return mat.has("RNA");
    }

    /**
     * @return {object} Object with the following properties:
     *
     * - `names`: Array of strings of length equal to the number of feature set collections, containing the names of the collections.
     * - `descriptions`: Array of strings of length equal to `names`, containing the descriptions for all collections.
     * - `species`: Array of strings of length equal to `names`, containing the taxonomy IDs for all collections.
     */
    fetchCollectionDetails() {
        return this.#manager.fetchCollectionDetails();
    }

    /**
     * @return {object} Object with the following properties:
     *
     * - `names`: Array of strings of length equal to the number of feature sets across all collections, containing the names of those sets.
     * - `descriptions`: Array of strings of length equal to `names`, containing the set descriptions.
     * - `sizes`: Int32Array of length equal to `names`, containing the set sizes.
     *   Each set's size is defined as the number of features in the dataset that are successfully mapped to a member of the set.
     * - `collections`: Int32Array of length equal to `names`, specifying the collection to which the set belongs.
     *   This is interpreted as the index of the arrays in {@linkcode fetchCollectionDetails}.
     */
    fetchSetDetails() {
        return this.#manager.fetchSetDetails();
    }

    /**
     * @return {number} Number of features from the input dataset that were successfully mapped to at least one gene in the reference database.
     */
    fetchUniverseSize() {
        return this.#manager.fetchUniverseSize();
    }

    /**
     * @param {external:ScoreMarkersResults} markers - Arbitrary marker detection results for an RNA modality, with the same order and identity of genes as from the upstream {@linkplain InputsState}.
     * This is most typically the output from {@linkcode MarkerDetectionState#fetchResults MarkerDetectionState.fetchResults} or equivalents from {@linkplain CustomSelectionsState}.
     * @param {number} group - Index of the group of interest inside `markers`.
     * @param {string} effect_size - Effect size to use for ranking.
     * This should be one of `"cohen"`, `"auc"`, `"lfc"` or `"delta_detected"`.
     * @param {string} summary - Summary statistic to use for ranking.
     * This should be one of `"min"`, `"mean"` or `"min_rank"`.
     *
     * @return {object} Object containing the following properties:
     *
     * - `set_ids`: Int32Array of length equal to the number of sets, containing the set IDs.
     *   Each entry is an index into the arrays returned by {@linkcode fetchSetDetails}.
     * - `counts`: Int32Array of length equal to `set_ids`, containing the number of markers present in each set.
     * - `pvalues`: Float64Array of length equal to `counts`, containing the enrichment p-values for each set.
     * - `num_markers`: number of markers selected for testing.
     */
    computeEnrichment(markers, group, effect_size, summary) {
        return this.#manager.computeEnrichment(group, effect_size, summary, markers, this.#parameters.top_markers);
    }

    /**
     * @param {number} set_id - Feature set ID, defined as an index into the arrays returned by {@linkcode fetchSetDetails}.
     *
     * @return {Int32Array} Array containing the row indices of the RNA count matrix corresponding to the genes in the specified set.
     */
    fetchFeatureSetIndices(set_id) {
        return this.#manager.fetchFeatureSetIndices(set_id);
    }

    /**
     * @param {number} set_id - Feature set ID, defined as an index into the arrays returned by {@linkcode fetchSetDetails}.
     *
     * @return {Object} Object containing:
     *
     * - `indices`: Int32Array containing the row indices of the genes in the set, relative to the RNA count matrix.
     * - `weights`: Float64Array containing the weights of each gene in the set.
     * - `scores`: Float64Array containing the feature set score for each cell.
     */
    computePerCellScores(set_id) {
        return this.#manager.computePerCellScores(set_id, this.#normalized.fetchNormalizedMatrix(), this.#filter.fetchFilteredBlock());
    }

    // Soft-deprecated.
    fetchPerCellScores(collection, set_index) {
        return this.computePerCellScores(collection, set_index);
    }

    /**
     * @return {object} Object containing the parameters.
     */
    fetchParameters() {
        return _fetchParameters(this.#parameters);
    }

    /****************************
     ******** Defaults **********
     ****************************/

    /**
     * @return {object} Default parameters that may be modified and fed into {@linkcode FeatureSetEnrichmentState#compute compute}.
     */
    static defaults() {
        return {
            skip: false,
            automatic: true,
            species: [],
            gene_id_column: null, 
            gene_id_type: "ENSEMBL", 
            top_markers: 100
        };
    }

    /***************************
     ******** Remotes **********
     ***************************/

    static flush() {
        return;
    }

    static setDownload(fun) {
        return FeatureSetManager.setDownload(fun);
    }

    /***************************
     ******** Compute **********
     ***************************/

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     *
     * @param {object} parameters - Parameter object, equivalent to the `feature_set_enrichment` property of the `parameters` of {@linkcode runAnalysis}.
     * @param {boolean} parameters.skip - Whether to skip the preparation of feature set collections.
     * If `true`, none of the other methods (e.g., {@linkcode computeEnrichment}, {@linkcode computePerCellScores}) should be called.
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
        this.changed = false;
        if (this.#inputs.changed) {
            this.changed = true;
        }

        let { skip, automatic, species, gene_id_column, gene_id_type, top_markers } = parameters;
        if (skip !== this.#parameters.skip) {
            this.changed = true;
        }

        if (this.valid() && !skip) {
            let modified = await _buildCollections(
                this.#parameters, 
                this.#manager,
                collections, 
                automatic, 
                species, 
                gene_id_column, 
                gene_id_type, 
                () => this.#inputs.fetchFeatureAnnotations()["RNA"],
                () => this.#inputs.guessRnaFeatureTypes()
            );
            if (modified) {
                this.changed = true;
            }

            if (top_markers !== this.#parameters.top_markers) {
                this.changed = true;
            }
        }

        _transplantParameters(this.#parameters, collections, automatic, species, gene_id_column, gene_id_type, top_markers);
        this.#parameters.skip = skip;
        return;
    }
}

/*****************************
 ******** Standalone *********
 *****************************/

/**
 * Standalone version of {@linkplain FeatureSetEnrichmentState} that provides the same functionality outside of {@linkcode runAnalysis}.
 * Users can supply their own annotation and marker results to compute the enrichment statistics for each group.
 * Users are also responsible for ensuring that the lifetime of the supplied objects exceeds that of the constructed instance,
 * i.e., the Wasm-related `free()` methods of the inputs are not called while the FeatureSetEnrichmentInstance is still in operation.
 */
export class FeatureSetEnrichmentStandalone {
    #annotations;
    #guesses;

    #normalized;
    #block;
    #backmap;

    #parameters;
    #manager;

    /**
     * @param {external:DataFrame} annotations - A {@linkplain external:DataFrame DataFrame} of per-gene annotations, where each row corresponds to a gene.
     * @param {object} [options={}] - Optional parameters.
     * @param {?(external:ScranMatrix)} [options.normalized=null] - A {@linkcode external:ScranMatrix ScranMatrix} of log-normalized expression values,
     * to be used in {@linkcode FeatureSetEnrichmentStandalone#computePerCellScores FeatureSetEnrichmentStandalone.computePerCellScores}.
     * Each row corresponds to a gene in the same order as `annotations`. 
     * @param {?(Array|TypedArray)} [options.block=null] - Array of length equal to the number of columns in `normalized`, containing the block assignments for each column. 
     * If `null`, all columns are assigned to the same block.
     */
    constructor(annotations, { normalized = null, block = null } = {}) {
        this.#annotations = annotations;
        this.#guesses = null;

        this.#normalized = null;
        this.#block = null;
        this.#backmap = null;

        if (normalized !== null) {
            if (normalized.numberOfRows() !== this.#annotations.numberOfRows()) {
                throw new Error("number of rows of 'annotations' and 'normalized' should be identical");
            }

            if (block !== null) {
                if (normalized.numberOfColumns() !== block.length) {
                    throw new Error("number of columns of 'normalized' should equal the length of 'block'");
                }

                let dump = utils.subsetInvalidFactors([ block ]);
                if (dump.retain !== null) {
                    this.#normalized = scran.subsetColumns(normalized, dump.retain);
                    this.#backmap = dump.retain;
                } else {
                    this.#normalized = normalized.clone();
                }

                this.#block = dump.arrays[0].ids;
            } else {
                this.#normalized = normalized.clone();
            }
        }

        this.#parameters = FeatureSetEnrichmentState.defaults();
        this.#manager = new FeatureSetManager; // empty Manager is valid for the defaults, which have no collections anyway.
    }

    #guessFeatureTypes() {
        if (this.#guesses == null) {
            this.#guesses = utils.guessFeatureTypes(this.#annotations);
        }
        return this.#guesses;
    }

    // Testing functions to check that the sanitization worked correctly.
    _peekMatrices() {
        return this.#normalized;
    }

    _peekBlock() {
        return this.#block;
    }

    /**
     * Frees all resources associated with this instance.
     */
    free() {
        scran.free(this.#block);
        scran.free(this.#normalized);
        this.#manager.free();
        return; // nothing extra to free here.
    }

    /**
     * If this method is not called, the parameters default to those in {@linkcode FeatureSetEnrichmentState#defaults FeatureSetEnrichmentState.defaults}.
     *
     * @param {object} parameters - Parameter object, see the argument of the same name in {@linkcode FeatureSetEnrichmentState#compute FeatureSetEnrichmentState.compute} for more details.
     *
     * @return The state is updated with new parameters.
     * @async
     */
    async setParameters(parameters) {
        let { automatic, species, gene_id_column, gene_id_type, top_markers } = parameters;

        await _buildCollections(
            this.#parameters,
            this.#manager,
            collections, 
            automatic, 
            species, 
            gene_id_column, 
            gene_id_type, 
            () => this.#annotations,
            () => this.#guessFeatureTypes()
        );

        _transplantParameters(this.#parameters, collections, automatic, species, gene_id_column, gene_id_type, top_markers);
    }

    /**
     * @return {object} Object containing the details about the available feature set collections,
     * see {@linkcode FeatureSetEnrichmentState#fetchCollectionDetails FeatureSetEnrichmentState.fetchCollectionDetails} for more details.
     */
    fetchCollectionDetails() {
        return this.#manager.fetchCollectionDetails();
    }

    /**
     * @return {object} Object containing the details about the available feature sets,
     * see {@linkcode FeatureSetEnrichmentState#fetchSetDetails FeatureSetEnrichmentState.fetchSetDetails} for more details.
     */
    fetchSetDetails() {
        return this.#manager.fetchSetDetails();
    }

    /**
     * @return {number} Number of features from the input dataset that were successfully mapped to at least one gene in the reference database.
     */
    fetchUniverseSize() {
        return this.#manager.fetchUniverseSize();
    }

    /**
     * @param {external:ScoreMarkersResults} markers - Marker detection results for an RNA modality.
     * @param {number} group - Group index of interest.
     * @param {string} effect_size - Effect size to use for ranking.
     * This should be one of `"cohen"`, `"auc"`, `"lfc"` or `"delta_detected"`.
     * @param {string} summary - Summary statistic to use for ranking.
     * This should be one of `"min"`, `"mean"` or `"min_rank"`.
     *
     * @return {object} Object containing statistics for the enrichment of the top marker genes in each feature set.
     * See {@linkcode FeatureSetEnrichmentState#computeEnrichment FeatureSetEnrichmentState.computeEnrichment} for more details.
     */
    computeEnrichment(markers, group, effect_size, summary) {
        return this.#manager.computeEnrichment(group, effect_size, summary, markers, this.#parameters.top_markers);
    }

    /**
     * @param {number} set_id - Feature set ID, defined as an index into the arrays returned by {@linkcode fetchSetDetails}.
     *
     * @return {Int32Array} Array containing the row indices of the RNA count matrix corresponding to the genes in the specified set.
     */
    fetchFeatureSetIndices(set_id) {
        return this.#manager.fetchFeatureSetIndices(set_id);
    }

    /**
     * @return {object} Object containing the parameters.
     */
    fetchParameters() {
        return _fetchParameters(this.#parameters);
    }

    /**
     * @param {number} set_id - Feature set ID, defined as an index into the arrays returned by {@linkcode fetchSetDetails}.
     *
     * @return {Object} Object containing the per-cell scores for the feature set activity.
     * See {@linkcode FeatureSetEnrichmentState#computePerCellScores FeatureSetEnrichmentState.computePerCellScores} for more details.
     */
    computePerCellScores(set_id) {
        if (this.#normalized == null) {
            throw new Error("no normalized matrix supplied in constructor");
        }

        let output = this.#manager.computePerCellScores(set_id, this.#normalized, this.#block);

        if (this.#backmap !== null) {
            let backfilled = new Float64Array(output.scores.length);
            backfilled.fill(Number.NaN);
            this.#backmap.forEach((x, i) => {
                backfilled[x] = output.scores[i];
            });
            output.scores = backfilled;
        }

        return output;
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
