import * as scran from "scran.js";
import * as utils from "./utils/general.js";
import * as markers from "./utils/markers.js";
import * as filter_module from "./cell_filtering.js";
import * as rna_norm_module from "./rna_normalization.js";
import * as adt_norm_module from "./adt_normalization.js";
import * as crispr_norm_module from "./crispr_normalization.js";

export const step_name = "custom_selections";

/**
 * Abstract class for handling custom selections.
 * Users should construct {@linkplain CustomSelectionsState} or {@linkplain CustomSelectionsStandalone} instances instead.
 * @hideconstructor
 */
export class CustomSelectionsCore {
    #generator;
    #cache;
    #parameters;

    constructor(generator, parameters = null, cache = null) {
        this.#generator = generator;
        this.#cache = (cache === null ? { "results": {} } : cache); 
        this.#parameters = (parameters === null ? { "selections": {} } : parameters);
        this.changed = false;
    }

    #liberate(i) {
        for (const [k, v] of Object.entries(this.#cache.results[i].raw)) {
            v.free();                                                
        }
    }

    /**
     * Frees all resources associated with this instance.
     */
    free() {
        utils.freeCache(this.#cache.buffer);
        for (const k of Object.keys(this.#cache.results)) {
            this.#liberate(k);
        }
        markers.freeVersusResults(this.#cache.versus);
    }

    /***************************
     ******** Setters **********
     ***************************/

    /**
     * Add a custom selection and compute its markers.
     * Users should run {@linkcode CustomSelectionsCore#compute compute} at least once before calling this function.
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
        let to_use = this.#generator.validModalities();
        let mat = this.#generator.matrix(to_use[0]);
        let ncells = mat.numberOfColumns();
        utils.checkIndices(selection, ncells);

        // Assumes that we have at least one cell in and outside the selection!
        var buffer = utils.allocateCachedArray(ncells, "Int32Array", this.#cache);
        buffer.fill(0);
        var tmp = buffer.array();
        selection.forEach(element => { tmp[element] = 1; });

        let res = {};
        for (const k of to_use) {
            let mat = this.#generator.matrix(k);
            res[k] = scran.scoreMarkers(mat, buffer, { 
                block: this.#generator.block(),
                lfcThreshold: this.#parameters.lfc_threshold,
                computeAuc: this.#parameters.compute_auc
            }); 
        }
              
        // Removing previous results, if there were any.
        if (id in this.#cache.results) {
            this.#liberate(id);
        }
      
        this.#cache.results[id] = { "raw": res };

        // making a copy to take ownership.
        if (copy) {
            selection = selection.slice();
        }
        this.#parameters.selections[id] = selection;
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
        this.#liberate(id);
        delete this.#cache.results[id];
        delete this.#parameters.selections[id];
        return;
    }

    /***************************
     ******** Getters **********
     ***************************/

    /**
     * @param {string} id - An identifier for the desired selection.
     *
     * @return {object} Object containing the markers for the desired selection.
     * Each key is a modality name while each value is a {@linkplain external:ScoreMarkersResults ScoreMarkersResults} object,
     * containing the marker detection results across all features of the corresponding modality.
     * The set of cells in the selection is denoted as group 1, while all cells outside of the selection are denoted as group 0.
     */
    fetchResults(id, feat_type) {
        return this.#cache.results[id].raw;
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
        let raw = this.#parameters.selections[id];
        if (copy) {
            raw = raw.slice();
        }
        return raw;
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
        let replacement = {};

        for (const [k, v] of Object.entries(this.#parameters.selections)) {
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

    /**
     * @return {object} Object containing the parameters.
     */
    fetchParameters() {
        return {
            lfc_threshold: this.#parameters.lfc_threshold,
            compute_auc: this.#parameters.compute_auc
        };
    }

    /***************************
     ******** Compute **********
     ***************************/

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
        if (this.#generator.allInvalid()) {
            for (const key of Object.keys(this.#cache.results)) {
                this.#liberate(key);
            }
            this.#parameters.selections = {};
            this.#cache.results = {};
            markers.freeVersusResults(this.#cache.versus);
            this.changed = true;
        }

        /* If the parameters changed, we recompute all the per-selection markers.
         * Technically we would need to re-run detection on the existing selections
         * if the normalization changed but the QC was the same. In practice, this
         * never happens, so we'll deal with it later.
         */
        if (lfc_threshold !== this.#parameters.lfc_threshold || compute_auc != this.#parameters.compute_auc) {
            this.#parameters.lfc_threshold = lfc_threshold;
            this.#parameters.compute_auc = compute_auc;
            for (const [key, value] of Object.entries(this.#parameters.selections)) {
                this.addSelection(key, value, { copy: false });
            }
            this.changed = true;
        }

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

    /*******************************
     ******** Versus mode **********
     *******************************/

    /**
     * Extract markers for a pairwise comparison between two selections, 
     * for more detailed examination of the differences between them.
     * Users should run {@linkcode CustomSelectionsCore#compute compute} at least once before calling this function.
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
    computeVersus(left, right, rank_type, feat_type) {
        // No need to free this afterwards; we don't own the normalized matrices anyway.
        let matrices = new scran.MultiMatrix;
        for (const modality of this.#generator.validModalities()) {
            matrices.add(modality, this.#generator.matrix(modality));
        }

        if (!("versus" in this.#cache)) {
            this.#cache["versus"] = {};
        }
        let cache = this.#cache.versus;

        let cache_info = markers.locateVersusCache(left, right, cache);
        let left_index = (cache_info.left_small ? 0 : 1);
        let right_index = (cache_info.left_small ? 1 : 0);

        if (cache_info.run) {
            let selections = this.#parameters.selections;
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

            let block = this.#generator.block();
            markers.computeVersusResults(matrices, new_clusters, block, keep, cache_info.cached, this.#parameters.lfc_threshold, this.#parameters.compute_auc);
        }

        return { 
            results: cache_info.cached,
            left: left_index,
            right: right_index
        };
    }
}

/*************************
 ****** Standalone *******
 *************************/

class StandaloneGenerator {
    #normalized;
    #block;
        
    constructor(normalized, block) {
        let N = null;
        for (const [k, v] of Object.entries(normalized)) {
            if (N != null) {
                if (v.numberOfColumns() != N) {
                    throw new Error("all matrices in 'normalized' should have the same number of columns as the length of 'groups'");
                }
            } else {
                N = v.numberOfColumns();
            }
        }
        this.#normalized = normalized;

        if (block !== null) {
            if (block.length && typeof block[0] !== "number") {
                throw new Error("'block' should encode each block as a non-negative integer starting from zero");
            }
            if (block.length != N) {
                throw new Error("'block' should have the same length as the number of columns in each entry of 'normalized'");
            }
        }
        this.#block = block;
    }

    free() {
        return;
    }

    validModalities() {
        return Object.keys(this.#normalized);
    }

    needsUpdate(modality) {
        return false;
    }

    allInvalid() {
        return false;
    }

    matrix(modality) {
        return this.#normalized[modality];
    }

    block() {
        return this.#block;
    }
}

/**
 * Standalone version of {@linkplain CustomSelectionsState} that provides the same functionality outside of {@linkcode runAnalysis}.
 * Users can supply their own normalized matrices and blocking factor to compute the various marker statistics for each custom selection.
 * Users are also responsible for ensuring that the lifetime of the supplied objects exceeds that of the constructed CustomSelectionsStandalone instance,
 * i.e., the Wasm-related `free()` methods are not called while the MarkerDetectionStandalone instance is still in operation.
 *
 * @extends CustomSelectionsCore
 */
export class CustomSelectionsStandalone extends CustomSelectionsCore {
    /**
     * @param {object} normalized - Object where each key is a modality name and each value is a {@linkcode external:ScranMatrix ScranMatrix} of log-normalized values.
     * Each ScranMatrix should have the same number of columns.
     * @param {object} [options={}] - Optional parameters.
     * @param {?(Array|TypedArray)} [options.block=null] - Array of length equal to the number of columns in any value of `normalized`.
     * This should contain the block assignments for each column, encoded as non-negative integers starting from zero.
     * If `null`, all columns are assigned to the same block.
     */
    constructor(normalized, { block = null } = {}) {
        let generator = new StandaloneGenerator(normalized, block);
        super(generator, null, null);
    }

    /**
     * See {@linkcode CustomSelectionsCore.defaults} for details.
     */
    static defaults() {
        return CustomSelectionsCore.defaults();
    }
}

/********************
 ****** State *******
 ********************/

class StateGenerator {
    #filter;
    #norm_states;
        
    constructor(filter, norm_states) {
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
    }

    free() {
        return;
    }

    validModalities() {
        return utils.findValidUpstreamStates(this.#norm_states);
    }

    needsUpdate(modality) {
        return this.#norm_states[modality].changed;
    }

    allInvalid() {
        return this.#filter.changed;
    }

    matrix(modality) {
        return this.#norm_states[modality].fetchNormalizedMatrix();
    }

    block() {
        return this.#filter.fetchFilteredBlock();
    }
}

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
 * @extends CustomSelectionsCore
 */
export class CustomSelectionsState extends CustomSelectionsCore {
    constructor(filter, norm_states, parameters = null, cache = null) {
        let gen = new StateGenerator(filter, norm_states);
        super(gen, parameters, cache);
    }

    /**
     * See {@linkcode CustomSelectionsCore.defaults} for details.
     */
    static defaults() {
        return CustomSelectionsCore.defaults();
    }

    /**
     * Soft-deprecated, construct a {@linkplain CustomSelectionsStandalone} instead.
     *
     * @param {string} left - Identifier of one selection.
     * @param {string} right - Identifier of another selection to be compared against `left`.
     * @param {external:MultiMatrix} matrices - A {@linkplain external:MultiMatrix MultiMatrix} object containing log-normalized matrices for each modality.
     * @param {object} selections - Object containing selections of cells.
     * Each key should be a selection identifier while each value is an Array, TypedArray or WasmArray.
     * Each array should contain integer column indices on `matrices`.
     * @param {object} [options={}] - Optional parameters.
     * @param {object} [options.cache={}] - Cache for the results.
     * @param {?Int32WasmArray} [options.block=null] - Blocking factor of length equal to the number of cells in `mat`, 
     * see {@linkcode CellFilteringState#fetchFilteredBlock CellFilteringState.fetchFilteredBlock}.
     * @param {number} [options.lfc_threshold=0] - Log-fold change threshold to use when computing the Cohen's d and AUC for each pairwise comparison.
     * @param {boolean} [options.compute_auc=true] - Whether to compute the AUCs.
     * Setting this to `false` will skip AUC calculations and improve speed and memory efficiency.
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
    static computeVersusCustom(left, right, matrices, selections, { cache = {}, block = null, lfc_threshold = 0, compute_auc = true } = {}) {
        let cache_info = markers.locateVersusCache(left, right, cache);
        let left_index = (cache_info.left_small ? 0 : 1);
        let right_index = (cache_info.left_small ? 1 : 0);

        if (cache_info.run) {
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

export function unserialize(handle, permuters, filter, norm_states) {
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

    return new CustomSelectionsState(filter, norm_states, parameters, cache);
}

