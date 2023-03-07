import * as scran from "scran.js";
import * as utils from "./utils/general.js";
import * as markers from "./utils/markers.js";
import * as filter_module from "./cell_filtering.js";
import * as rna_norm_module from "./rna_normalization.js";
import * as adt_norm_module from "./adt_normalization.js";
import * as crispr_norm_module from "./crispr_normalization.js";
import * as enrich_module from "./feature_set_enrichment.js";

export const step_name = "custom_selections";

/************************
 ****** Internals *******
 ************************/

class SelectionManager {
    constructor(selections = null, cache = null) {
        this._selections = (selections == null ? {} : selections);
        this._cache = (cache == null ? { results: {} } : cache);
    }

    #liberate(i) {
        for (const [k, v] of Object.entries(this._cache.results[i].raw)) {
            v.free();                                                
        }
    }

    free() {
        utils.freeCache(this._cache.buffer);
        delete this._cache.buffer;

        for (const k of Object.keys(this._cache.results)) {
            this.#liberate(k);
        }
        this._cache.results = {};        

        markers.freeVersusResults(this._cache.versus);
        delete this._cache.versus;
    }

    addSelection(id, selection, to_use, matfun, block, copy, lfc_threshold, compute_auc) {
        let mat = matfun(to_use[0]);
        let ncells = mat.numberOfColumns();
        utils.checkIndices(selection, ncells);

        // Assumes that we have at least one cell in and outside the selection!
        var buffer = utils.allocateCachedArray(ncells, "Int32Array", this._cache);
        buffer.fill(0);
        var tmp = buffer.array();
        selection.forEach(element => { tmp[element] = 1; });

        let res = {};
        for (const k of to_use) {
            let mat = matfun(k);
            res[k] = scran.scoreMarkers(mat, buffer, { block: block, lfcThreshold: lfc_threshold, computeAuc: compute_auc }); 
        }
              
        // Removing previous results, if there were any.
        if (id in this._cache.results) {
            this.#liberate(id);
        }
      
        this._cache.results[id] = { "raw": res };

        // making a copy to take ownership.
        if (copy) {
            selection = selection.slice();
        }
        this._selections[id] = selection;
        return;
    }

    removeSelection(id) {
        this.#liberate(id);
        delete this._cache.results[id];
        delete this._selections[id];
        return;
    }

    fetchResults(id) {
        return this._cache.results[id].raw;
    }

    fetchSelectionIndices(id, { copy = true } = {}) {
        let raw = this._selections[id];
        if (copy) {
            raw = raw.slice();
        }
        return raw;
    }

    fetchSelections({ copy = true, force = null } = {}) {
        let replacement = {};

        for (const [k, v] of Object.entries(this._selections)) {
            let store = v;
            let needs_copy = copy;

            if (force !== null) {
                if (force == "Array") {
                    if (!(v instanceof Array)) {
                        store = Array.from(v);
                        needs_copy = false;
                    }
                } else if (force == "Int32Array") {
                    if (!(v instanceof Int32Array)) {
                        store = new Int32Array(v);
                        needs_copy = false;
                    }
                }
            } 

            if (needs_copy) {
                store = store.slice();
            }
            replacement[k] = store;
        }
        return replacement;        
    }

    computeVersus(left, right, to_use, matfun, block, lfc_threshold, compute_auc) {
        if (!("versus" in this._cache)) {
            this._cache["versus"] = {};
        }
        let cache = this._cache.versus;

        let cache_info = markers.locateVersusCache(left, right, cache);
        let left_index = (cache_info.left_small ? 0 : 1);
        let right_index = (cache_info.left_small ? 1 : 0);

        if (cache_info.run) {
            // No need to free this afterwards; we don't own the normalized matrices anyway.
            let matrices = new scran.MultiMatrix;
            for (const modality of to_use) {
                matrices.add(modality, matfun(modality));
            }

            let selections = this._selections;
            if (!(left in selections && right in selections)) {
                throw new Error("invalid selection ID requested in versus mode");
            }

            let leftsel = selections[left];
            let rightsel = selections[right];
            if (leftsel.length == 0 || rightsel.length == 0) {
                throw new Error("non-zero entries should be present for both requested selections in versus mode");
            }

            let triplets = [];
            leftsel.forEach(x => {
                triplets.push({ "index": x, "cluster": left_index });
            });
            rightsel.forEach(x => {
                triplets.push({ "index": x, "cluster": right_index });
            });

            triplets.sort((a, b) => a.index - b.index);
            let keep = triplets.map(x => x.index);
            let new_clusters = triplets.map(x => x.cluster);
            markers.computeVersusResults(matrices, new_clusters, block, keep, cache_info.cached, lfc_threshold, compute_auc);
        }

        return { 
            results: cache_info.cached,
            left: left_index,
            right: right_index
        };
    }
}

/********************
 ****** State *******
 ********************/

/**
 * Applications can perform marker detection on custom selections of cells.
 * This allows users to dynamically select cells on a UI and quickly obtain a list of distinguishing markers for that selection.
 * This wraps the [`scoreMarkers`](https://kanaverse.github.io/scran.js/global.html#scoreMarkers) function 
 * from [**scran.js**](https://github.com/kanaverse/scran.js).
 *
 * Users should not construct these instances manually; instead, they are automatically assembled by {@linkcode createAnalysis}.
 * Similarly, users should not directly call the {@linkcode CustomSelectionsCore#compute compute} method, which is instead invoked by {@linkcode runAnalysis}.
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class CustomSelectionsState {
    #filter;
    #norm_states;
    #enrichment;

    #manager;
    #parameters;

    constructor(filter, norm_states, enrich, parameters = null, cache = null) {
        if (!(filter instanceof filter_module.CellFilteringState)) {
            throw new Error("'filter' should be a CellFilteringState object");
        }
        this.#filter = filter;

        if (!(norm_states.RNA instanceof rna_norm_module.RnaNormalizationState)) {
            throw new Error("'norm_states.RNA' should be an RnaNormalizationState object");
        }
        if (!(norm_states.ADT instanceof adt_norm_module.AdtNormalizationState)) {
            throw new Error("'norm_states.ADT' should be an AdtNormalizationState object");
        }
        if (!(norm_states.CRISPR instanceof crispr_norm_module.CrisprNormalizationState)) {
            throw new Error("'norm_states.CRISPR' should be an CrisprNormalizationState object");
        }
        this.#norm_states = norm_states;

        if (!(enrich instanceof enrich_module.FeatureSetEnrichmentState)) {
            throw new Error("'enrich' should be an FeatureSetEnrichmentState object");
        }
        this.#enrichment = enrich;

        let selections = null;
        if (parameters !== null && "selections" in parameters) {
            selections = parameters.selections;
        }

        this.#manager = new SelectionManager(selections, cache);
        this.#parameters = {};
        this.changed = false;
    }

    /**
     * Frees all resources associated with this instance.
     */
    free() {
        this.#manager.free();
        return;
    }

    /**
     * @return {object} Object containing the parameters.
     */
    fetchParameters() {
        return {
            lfc_threshold: this.#parameters.lfc_threshold,
            compute_auc: this.#parameters.compute_auc
        };
    }

    /**
     * Add a custom selection and compute its markers.
     * It is assumed that {@linkcode runAnalysis} was already run on this instance before calling this method.
     *
     * @param {string} id A unique identifier for the new custom selection.
     * @param {Array|TypedArray} selection The indices of the cells in the selection.
     * Indices should refer to positions of cells in the QC-filtered matrix, not the original matrix.
     * @param {object} [options] - Optional parameters.
     * @param {boolean} [options.copy=true] - Whether to make a copy of `selection` before storing it inside this object.
     * If `false`, it is assumed that the caller makes no further modifications to the passed `selection`.
     *
     * @return The custom selection is added to the state and calculation of its markers is performed.
     * Nothing is returned.
     */
    addSelection(id, selection, { copy = true } = {}) {
        this.#addSelection(id, selection, copy, this.#parameters.lfc_threshold, this.#parameters.compute_auc);
    }

    #addSelection(id, selection, copy, lfc_threshold, compute_auc) {
        let to_use = utils.findValidUpstreamStates(this.#norm_states);
        this.#manager.addSelection(
            id, 
            selection, 
            to_use, 
            modality => this.#norm_states[modality].fetchNormalizedMatrix(),
            this.#filter.fetchFilteredBlock(),
            copy,
            lfc_threshold,
            compute_auc
        );
        return;
    }

    /**
     * Remove a custom selection and its results from the state.
     *
     * @param {string} id - An identifier for the selection to be removed.
     *
     * @return The specified selection and its results are removed from the state.
     * Nothing is returned.
     */
    removeSelection(id) {
        this.#manager.removeSelection(id);
        return;
    }

    /**
     * @param {string} id - An identifier for the desired selection.
     *
     * @return {object} Object containing the markers for the desired selection.
     * Each key is a modality name while each value is a {@linkplain external:ScoreMarkersResults ScoreMarkersResults} object,
     * containing the marker detection results across all features of the corresponding modality.
     * The set of cells in the selection is denoted as group 1, while all cells outside of the selection are denoted as group 0.
     */
    fetchMarkers(id) {
        return this.#manager.fetchResults(id);
    }

    // Soft-deprecated.
    fetchResults(id) {
        return this.fetchMarkers(id);
    }

    /**
     * Compute enrichment among the set of upregulated markers for a custom selection.
     *
     * @param {string} id - An identifier for the desired selection.
     * @param {string} effect_size - Effect size to use for ranking.
     * This should be one of `"cohen"`, `"auc"`, `"lfc"` or `"delta_detected"`.
     *
     * @return {object} Object where each entry corresponds to a feature set collection and contains enrichment statistics for all sets in that collection.
     * See {@linkcode FeatureSetEnrichment#computeEnrichment FeatureSetEnrichment.computeEnrichment} for details.
     */
    computeEnrichment(id, effect_size) {
        let x = this.fetchMarkers(id);
        if (!("RNA" in x)) {
            return null;
        }
        return this.#enrichment.computeEnrichment(x["RNA"], 0, effect_size, "mean"); // two groups, so the effect size doesn't really matter.
    }

    /**
     * Retrieve the indices for a selection of interest.
     *
     * @param {string} id - The identifier for the selection.
     * @param {object} [options] - Optional parameters.
     * @param {boolean} [options.copy=true] - Whether to make a copy of `selection` before returning it.
     * If `false`, it is assumed that the caller does not modify the selection.
     *
     * @return {Array|TypedArray} Array of indices in the requested selection.
     * Note that indices are relative to the filtered matrix - 
     * use {@linkcode CellFilteringState#undoFiltering CellFilteringState.undoFiltering} to convert them to indices on the original dataset.
     */
    fetchSelectionIndices(id, { copy = true } = {}) {
        return this.#manager.fetchSelectionIndices(id, { copy });
    }

    /**
     * Retrieve indices for all selections.
     *
     * @param {object} [options] - Optional parameters.
     * @param {boolean} [options.copy=true] - Whether to make a copy of `selection` before returning it.
     * If `false`, it is assumed that the caller does not modify the selection.
     * @param {?string} [force=null] - Whether to force each `selection` to be an `"Array"` or "`Int32Array"`.
     * If `null`, the existing type of each selection is used.
     *
     * @return {object} Object where the keys are the selection names and the values are arrays of indices for each selection.
     * Each array is a copy and can be modified without affecting the CustomSelectionsState.
     * See {@linkcode CustomSelectionsState#fetchSelectionIndices fetchSelectionIndices} for more details on the interpretation of the indices.
     */
    fetchSelections({ copy = true, force = null } = {}) {
        return this.#manager.fetchSelections({ copy, force });
    }

    /**
     * @param {object} parameters - Parameter object, equivalent to the `custom_selections` property of the `parameters` of {@linkcode runAnalysis}.
     * @param {number} parameters.lfc_threshold - Log-fold change threshold to use when computing the Cohen's d and AUC for each pairwise comparison.
     * @param {boolean} parameters.compute_auc - Whether to compute the AUCs.
     * Setting this to `false` will skip AUC calculations and improve speed and memory efficiency.
     *
     * @return The state is updated by removing stale selections if the QC filter was altered.
     */
    compute(parameters) {
        let { lfc_threshold, compute_auc } = parameters;
        this.changed = false;

        /* If the QC filter was re-run, all of the selections are invalidated as
         * the identity of the indices may have changed.
         */
        if (this.#filter.changed) {
            this.#manager.free();
            this.#manager = new SelectionManager;
            this.changed = true;
        }

        /* If the parameters changed, we recompute all the per-selection markers.
         * Technically we would need to re-run detection on the existing selections
         * if the normalization changed but the QC was the same. In practice, this
         * never happens, so we'll deal with it later.
         */
        if (lfc_threshold !== this.#parameters.lfc_threshold || compute_auc != this.#parameters.compute_auc) {
            for (const [key, value] of Object.entries(this.#manager._selections)) {
                this.#addSelection(key, value, false, lfc_threshold, compute_auc);
            }
            this.#parameters.lfc_threshold = lfc_threshold;
            this.#parameters.compute_auc = compute_auc;
            this.changed = true;
        }

        /**
         * The status of #enrichment doesn't really affect any of the cached
         * resources because we don't cache the enrichment results anyway, so
         * there's no need to do anything if this.#enrichment.changed == true.
         */

        return;
    }

    /**
     * @return {object} Default parameters that may be modified and fed into {@linkcode MarkerDetectionCore#compute compute}.
     */
    static defaults() {
        return {
            lfc_threshold: 0,
            compute_auc: true
        };
    }

    /**
     * Extract markers for a pairwise comparison between two selections for more detailed examination of the differences between them.
     *
     * @param {string} left - Identifier of one selection in which to find upregulated markers.
     * @param {string} right - Identifier of another selection to be compared against `left`.
     *
     * @return {object} Object containing:
     *
     * - `results`: object containing the marker statistics for the comparison between two clusters.
     *    Each key is a modality name and each value is a {@linkplain external:ScoreMarkersResults ScoreMarkersResults} object.
     * - `left`: index of the group corresponding to the `left` selection in each ScoreMarkersResults object.
     *    e.g., Cohen's d for the RNA markers of the `left` selection are defined as `output.results.RNA.cohen(output.left)`.
     * - `right`: index of the group corresponding to the `right` selection in each ScoreMarkersResults object.
     *    e.g., Cohen's d for the RNA markers of the `right` selection are defined as `output.results.RNA.cohen(output.right)`.
     */
    computeVersusMarkers(left, right) {
        let to_use = utils.findValidUpstreamStates(this.#norm_states);
        return this.#manager.computeVersus(
            left, 
            right, 
            to_use,
            modality => this.#norm_states[modality].fetchNormalizedMatrix(),
            this.#filter.fetchFilteredBlock(),
            this.#parameters.lfc_threshold,
            this.#parameters.compute_auc
        );
    }

    // Soft-deprecated.
    computeVersus(left, right) {
        return this.computeVersusMarkers(left, right);
    }

    /**
     * Compute the feature set enrichment among the set of markers upregulated in one selection compared to another.
     *
     * @param {string} left - Identifier of one selection in which to find upregulated markers.
     * @param {string} right - Identifier of another selection to be compared against `left`.
     * @param {string} effect_size - Effect size to use for ranking.
     * This should be one of `"cohen"`, `"auc"`, `"lfc"` or `"delta_detected"`.
     *
     * @return {object} Object where each entry corresponds to a feature set collection and contains enrichment statistics for all sets in that collection.
     * See {@linkcode FeatureSetEnrichment#computeEnrichment FeatureSetEnrichment.computeEnrichment} for details.
     */
    computeVersusEnrichment(left, right, effect_size) {
        let x = this.computeVersusMarkers(left, right);
        if (!("RNA" in x.results)) {
            return null;
        }
        return this.#enrichment.computeEnrichment(x.results["RNA"], x.left, effect_size, "mean"); // two groups, so the effect size doesn't really matter.
    }
}

/*************************
 ****** Standalone *******
 *************************/

/**
 * Standalone version of {@linkplain CustomSelectionsState} that provides the same functionality outside of {@linkcode runAnalysis}.
 * Users can supply their own normalized matrices and blocking factor to compute the various marker statistics for each custom selection.
 * Users are also responsible for ensuring that the lifetime of the supplied objects exceeds that of the constructed CustomSelectionsStandalone instance,
 * i.e., the Wasm-related `free()` methods are not called while the MarkerDetectionStandalone instance is still in operation.
 */
export class CustomSelectionsStandalone {
    #normalized;
    #block;
    #block_levels;

    #manager;
    #parameters;

    #missing_map;

    /**
     * @param {external:MultiMatrix} normalized - A {@linkplain external:MultiMatrix MultiMatrix} of log-normalized values for multiple modalities.
     * @param {object} [options={}] - Optional parameters.
     * @param {?(Array|TypedArray)} [options.block=null] - Array of length equal to the number of columns in any value of `normalized`, containing the block assignments for each column.
     * If `null`, all columns are assigned to the same block.
     */
    constructor(normalized, { block = null } = {}) {
        let N = null;
        for (const k of normalized.available()) {
            let alt = normalized.get(k).numberOfColumns();
            if (N != null) {
                if (alt != N) {
                    throw new Error("all matrices in 'normalized' should have the same number of columns as the length of 'groups'");
                }
            } else {
                N = alt;
            }
        }

        this.#block = null;
        this.#block_levels = null;
        this.#missing_map = null;

        if (block !== null) {
            if (block.length != N) {
                throw new Error("'block' should have the same length as the number of columns in each entry of 'normalized'");
            }

            let dump = utils.subsetInvalidFactors([ block ]);
            if (dump.retain !== null) {
                let revmap = new Int32Array(N);
                revmap.fill(-1);
                dump.retain.forEach((y, i) => { revmap[y] = i; });
                this.#missing_map = { to: revmap, from: dump.retain };

                let new_matrices = new scran.MultiMatrix;
                let temp = scran.createInt32WasmArray(dump.retain.length);
                try {
                    temp.set(dump.retain);
                    for (const k of normalized.available()) {
                        new_matrices.add(k, scran.subsetColumns(normalized.get(k), temp))
                    }
                } catch (e) {
                    scran.free(new_matrices);
                    throw e;
                } finally {
                    scran.free(temp);
                }

                this.#normalized = new_matrices;
            } else {
                this.#normalized = normalized.clone();
            }

            this.#block = dump.arrays[0].ids;
            this.#block_levels = dump.arrays[0].levels;
        } else {
            this.#normalized = normalized.clone();
        }

        this.#manager = new SelectionManager;
        this.#parameters = CustomSelectionsState.defaults();
        this.changed = false;
    }

    /**
     * Frees all resources associated with this instance.
     */
    free() {
        this.#manager.free();
        return;
    }

    /**
     * @return {object} Object containing the parameters.
     */
    fetchParameters() {
        return {
            lfc_threshold: this.#parameters.lfc_threshold,
            compute_auc: this.#parameters.compute_auc
        };
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.copy=true] - Whether to copy the return value on output.
     * Set to `false` for greater efficiency in strictly read-only applications.
     *
     * @return {Array} Array of levels for the blocking factor.
     * Block indices in the {@linkplain external:ScoreMarkersResults ScoreMarkersResults} instances returned by {@linkcode fetchResults} can be cross-referenced to this array.
     */
    fetchBlockLevels({ copy = true } = {}) {
        let ret = this.#block_levels;
        return (copy ? ret.slice() : ret);
    }

    // Testing functions to check that the sanitization worked correctly.
    _peekMatrices() {
        return this.#normalized;
    }

    _peekBlock() {
        return this.#block;
    }

    /**
     * If this method is not called, the parameters default to those in {@linkcode CustomSelectionsState#defaults CustomSelectionsState.defaults}.
     *
     * @param {object} parameters - Parameter object, see the argument of the same name in {@linkcode CustomSelectionsState#compute CustomSelectionsState.compute} for more details.
     *
     * @return The state is updated with the new parameters.
     */
    setParameters(parameters) {
        if (this.#parameters.lfc_threshold !== parameters.lfc_threshold || this.#parameters.compute_auc !== parameters.compute_auc) {
            this.free();
        }
        this.#parameters = parameters;
        return;
    }

    /**
     * Add a custom selection and compute its markers.
     * Users should run {@linkcode CustomSelectionsCore#compute compute} at least once before calling this function.
     *
     * @param {string} id - A unique identifier for the new custom selection.
     * @param {Array|TypedArray} selection - The indices of the cells in the selection.
     * @param {object} [options] - Optional parameters.
     * @param {boolean} [options.copy=true] - Whether to make a copy of `selection` before storing it inside this object.
     * If `false`, it is assumed that the caller makes no further modifications to the passed `selection`.
     *
     * @return The custom selection is added to the state and calculation of its markers is performed.
     * Nothing is returned.
     */
    addSelection(id, selection, { copy = true } = {}) {
        let selection_internal = selection;

        // Removing the invalid observations.
        if (this.#missing_map !== null) { 
            let collected = [];
            let revmap = this.#missing_map.to;
            for (const i of selection) {
                let j = revmap[i];
                if (j >= 0) {
                    collected.push(j);
                }
            }
            selection_internal = collected;
        }

        this.#manager.addSelection(
            id, 
            selection_internal,
            this.#normalized.available(),
            modality => this.#normalized.get(modality),
            this.#block,
            copy,
            this.#parameters.lfc_threshold,
            this.#parameters.compute_auc
        );

        return;
    }

    /**
     * Remove a custom selection and its results from the state.
     *
     * @param {string} id - An identifier for the selection to be removed.
     *
     * @return The specified selection and its results are removed from the state.
     * Nothing is returned.
     */
    removeSelection(id) {
        this.#manager.removeSelection(id);
        return;
    }

    /**
     * @param {string} id - An identifier for the desired selection.
     *
     * @return {object} Object containing the markers for the desired selection.
     * Each key is a modality name while each value is a {@linkplain external:ScoreMarkersResults ScoreMarkersResults} object,
     * containing the marker detection results across all features of the corresponding modality.
     * The set of cells in the selection is denoted as group 1, while all cells outside of the selection are denoted as group 0.
     */
    fetchResults(id) {
        return this.#manager.fetchResults(id);
    }

    #unmap(ids) {
        // Restoring the indices after adjusting for the invalid observations,
        // so that users get back indices relative to the input matrices.
        if (this.#missing_map !== null) {
            ids.forEach((x, i) => {
                ids[i] = this.#missing_map.from[x];
            });
        }
        return;
    }

    /**
     * Retrieve the indices for a selection of interest.
     *
     * @param {string} id - The identifier for the selection.
     * @param {object} [options] - Optional parameters.
     * @param {boolean} [options.copy=true] - Whether to make a copy of `selection` before returning it.
     * If `false`, it is assumed that the caller does not modify the selection.
     *
     * @return {Array|TypedArray} Array of indices in the requested selection.
     */
    fetchSelectionIndices(id, { copy = true } = {}) {
        let output = this.#manager.fetchSelectionIndices(id, { copy });
        this.#unmap(output);
        return output;
    }

    /**
     * Retrieve indices for all selections.
     *
     * @param {object} [options] - Optional parameters.
     * @param {boolean} [options.copy=true] - Whether to make a copy of `selection` before returning it.
     * If `false`, it is assumed that the caller does not modify the selection.
     * @param {?string} [force=null] - Whether to force each `selection` to be an `"Array"` or "`Int32Array"`.
     * If `null`, the existing type of each selection is used.
     *
     * @return {object} Object where the keys are the selection names and the values are arrays of indices for each selection.
     * Each array is a copy and can be modified without affecting the CustomSelectionsState.
     */
    fetchSelections({ copy = true, force = null } = {}) {
        let output = this.#manager.fetchSelections({ copy, force });
        for (const [k, v] of Object.entries(output)) {
            this.#unmap(v);
        }
        return output;
    }

    /**
     * Extract markers for a pairwise comparison between two selections, 
     * for more detailed examination of the differences between them.
     *
     * @param {string} left - Identifier of one selection in which to find upregulated markers.
     * @param {string} right - Identifier of another selection to be compared against `left`.
     *
     * @return {object} Object containing:
     *
     * - `results`: object containing the marker statistics for the comparison between two clusters.
     *    Each key is a modality name and each value is a {@linkplain external:ScoreMarkersResults ScoreMarkersResults} object.
     * - `left`: index of the group corresponding to the `left` selection in each ScoreMarkersResults object.
     *    e.g., Cohen's d for the RNA markers of the `left` selection are defined as `output.results.RNA.cohen(output.left)`.
     * - `right`: index of the group corresponding to the `right` selection in each ScoreMarkersResults object.
     *    e.g., Cohen's d for the RNA markers of the `right` selection are defined as `output.results.RNA.cohen(output.right)`.
     */
    computeVersus(left, right) {
        return this.#manager.computeVersus(
            left, 
            right, 
            this.#normalized.available(),
            modality => this.#normalized.get(modality),
            this.#block,
            this.#parameters.lfc_threshold,
            this.#parameters.compute_auc
        );
    }
}

/**************************
 ******** Loading *********
 **************************/

function fill_results(stats) {
    let ngenes = stats.means.length;
    let object = scran.emptyScoreMarkersResults(ngenes, 
        /* number of groups */ 2, 
        /* number of blocks */ 1,
        { computeAuc: ("auc" in stats) }
    );

    object.means(1, { fillable: true }).set(stats.means);
    object.detected(1, { fillable: true }).set(stats.detected);

    for (const index of Object.values(markers.summaries2int)) {
        object.cohen(1, { summary: index, fillable: true }).set(stats.cohen);
        object.lfc(1, { summary: index, fillable: true }).set(stats.lfc);
        object.deltaDetected(1, { summary: index, fillable: true }).set(stats.delta_detected);
        if ("auc" in stats) {
            object.auc(1, { summary: index, fillable: true }).set(stats.auc);
        }
    }

    return object;
}

export function unserialize(handle, permuters, filter, norm_states, enrich) {
    let ghandle = handle.open("custom_selections");

    let parameters = CustomSelectionsState.defaults();
    {
        let phandle = ghandle.open("parameters");
        if ("lfc_threshold" in phandle.children) {
            parameters.lfc_threshold = phandle.open("lfc_threshold", { load: true }).values[0];
        }
        if ("compute_auc" in phandle.children) {
            parameters.compute_auc = phandle.open("compute_auc", { load: true }).values[0] > 0;
        }

        parameters.selections = {};
        let shandle = phandle.open("selections");
        for (const key of Object.keys(shandle.children)) {
            let vals = shandle.open(key, { load: true }).values;

            // v1 wasn't sorted, so we make sure to sort things.
            for (var i = 1; i < vals.length; i++) {
                if (vals[i] < vals[i-1]) {
                    vals.sort();
                    break;
                }
            }

            parameters.selections[key] = vals;
        }
    }

    let cache = { results: {} };
    {
        let rhandle = ghandle.open("results");

        if ("markers" in rhandle.children) {
            // before v2.0
            let mhandle = rhandle.open("markers");
            for (const sel of Object.keys(mhandle.children)) {
                let current = markers.unserializeGroupStats(mhandle.open(sel), permuters["RNA"], { no_summaries: true, compute_auc: parameters.compute_auc });
                cache.results[sel] = { raw: { RNA: fill_results(current) } };
            }
        } else {
            // after v2.0.
            let phandle = rhandle.open("per_selection");
            for (const sel of Object.keys(phandle.children)) {
                let shandle = phandle.open(sel);
                let collected = {};
                for (const feat of Object.keys(shandle.children)) {
                    let current = markers.unserializeGroupStats(shandle.open(feat), permuters[feat], { no_summaries: true, compute_auc: parameters.compute_auc });
                    collected[feat] = fill_results(current);
                }
                cache.results[sel] = { raw: collected };
            }
        }
    }

    return new CustomSelectionsState(filter, norm_states, enrich, parameters, cache);
}

