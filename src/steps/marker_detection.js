import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as markers from "./utils/markers.js";
import * as filter_module from "./cell_filtering.js";
import * as choice_module from "./choose_clustering.js";
import * as rna_norm_module from "./rna_normalization.js";
import * as adt_norm_module from "./adt_normalization.js";
import * as crispr_norm_module from "./crispr_normalization.js";

export const step_name = "marker_detection";

/**
 * Results of marker detection,
 * see [here](https://www.kanaverse.github.io/scran.js/ScoreMarkersResults.html) for details.
 *
 * @external ScoreMarkersResults
 */

/*********************
 ***** Internals *****
 *********************/

function _free_results(cache) {
    for (const v of Object.values(cache.raw)) {
        utils.freeCache(v);
    }
    cache.raw = {};
    markers.freeVersusResults(cache.versus);
    delete cache.versus;
}

function _computeVersus(left, right, matrices, clusters, { cache = {}, block = null, lfc_threshold = 0, compute_auc = true } = {}) {
    let cache_info = markers.locateVersusCache(left, right, cache);
    let left_index = (cache_info.left_small ? 0 : 1);
    let right_index = (cache_info.left_small ? 1 : 0);

    if (cache_info.run) {
        let new_clusters = [];
        let keep = [];
        let leftfound = false, rightfound = false;
        clusters.forEach((x, i) => {
            if (x == left) {
                new_clusters.push(left_index);
                keep.push(i);
                leftfound = true;
            } else if (x == right) {
                new_clusters.push(right_index);
                keep.push(i);
                rightfound = true;
            }
        });

        if (!leftfound || !rightfound) {
            throw new Error("non-zero entries should be present for both requested clusters in versus mode");
        }

        markers.computeVersusResults(matrices, new_clusters, block, keep, cache_info.cached, lfc_threshold, compute_auc);
    }

    return { 
        results: cache_info.cached,
        left: left_index,
        right: right_index
    };
}

/*****************
 ***** State *****
 *****************/

/**
 * This step performs marker detection for each cluster of cells by performing pairwise comparisons to each other cluster.
 * This wraps the [`scoreMarkers`](https://kanaverse.github.io/scran.js/global.html#scoreMarkers) function 
 * from [**scran.js**](https://github.com/kanaverse/scran.js).
 * The clustering is obtained from the upstream {@linkplain ChooseClusteringState}.
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class MarkerDetectionState {
    #filter;
    #norm_states;
    #choice;
    #parameters;
    #cache;

    constructor(filter, norm_states, choice, parameters = null, cache = null) {
        if (!(filter instanceof filter_module.CellFilteringState)) {
            throw new Error("'filter' should be a State object from './cell_filtering.js'");
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

        if (!(choice instanceof choice_module.ChooseClusteringState)) {
            throw new Error("'choice' should be a State object from './choose_clustering.js'");
        }
        this.#choice = choice;

        this.#parameters = (parameters === null ? {} : parameters);
        this.#cache = (cache === null ? { "raw": {} } : cache);
        this.changed = false;
    }

    /**
     * Frees all resources associated with this instance.
     */
    free() {
        _free_results(this.#cache);
    }

    /**
     * Obtain the marker detection results for all modalities.
     * It is assumed that {@linkcode runAnalysis} was already run on this MarkerDetectionState instance before calling this method.
     *
     * @return {object} Object where each key is a modality name and each value is an {@linkplain external:ScoreMarkersResults ScoreMarkerResults} object,
     * containing marker detection statistics for all clusters.
     */
    fetchResults() {
        return this.#cache.raw;
    }

    /**
     * Parameters used to compute the results.
     * It is assumed that {@linkcode runAnalysis} was already run on this MarkerDetectionState instance before calling this method.
     *
     * @return {object} Object containing the parameters.
     */
    fetchParameters() {
        return { ...this.#parameters }; // avoid pass-by-reference links.
    }

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     *
     * @param {object} parameters - Parameter object, equivalent to the `marker_detection` property of the `parameters` of {@linkcode runAnalysis}.
     * @param {number} parameters.lfc_threshold - Log-fold change threshold to use when computing the Cohen's d and AUC for each pairwise comparison.
     * @param {boolean} parameters.compute_auc - Whether to compute the AUCs.
     * Setting this to `false` will skip AUC calculations and improve speed and memory efficiency.
     *
     * @return The state is updated with new results.
     */
    compute(parameters) {
        let { lfc_threshold, compute_auc } = parameters;
        this.changed = false;
        let changed_params = (lfc_threshold !== this.#parameters.lfc_threshold || compute_auc !== this.#parameters.compute_auc);
        
        for (const [k, v] of Object.entries(this.#norm_states)) {
            if (!v.valid()) {
                continue;
            }

            if (this.#choice.changed || v.changed || changed_params) {
                var mat = v.fetchNormalizedMatrix();
                var clusters = this.#choice.fetchClusters();
                var block = this.#filter.fetchFilteredBlock();
                
                utils.freeCache(this.#cache.raw[k]);
                this.#cache.raw[k] = scran.scoreMarkers(mat, clusters, { block: block, threshold: lfc_threshold, computeAuc: compute_auc });

                this.changed = true;
            }
        }

        this.#parameters.lfc_threshold = lfc_threshold;
        this.#parameters.compute_auc = compute_auc;
        if (this.changed) {
            markers.freeVersusResults(this.#cache.versus);
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

    // Soft-deprecated.
    static computeVersusCustom(left, right, matrices, clusters, { cache = {}, block = null, lfc_threshold = 0, compute_auc = true } = {}) {
        return computeVersus(left, right, matrices, clusters, { cache, block, lfc_threshold, compute_auc });
    }

    /**
     * Extract markers for a pairwise comparison between two clusters for more detailed examination of the differences between them.
     * It is assumed that {@linkcode runAnalysis} was already run on this MarkerDetectionState instance before calling this method.
     *
     * @param {number} left - Index of one cluster in which to find upregulated markers.
     * @param {number} right - Index of another cluster to be compared against `left`.
     *
     * @return {object} Object containing:
     *
     * - `results`: object containing the marker statistics for the comparison between two clusters.
     *    Each key is a modality name and each value is a {@linkplain external:ScoreMarkersResults ScoreMarkersResults} object.
     * - `left`: index of the group corresponding to the `left` cluster in each ScoreMarkersResults object,
     *    e.g., Cohen's d for the RNA markers of the `left` cluster are defined as `output.results.RNA.cohensD(output.left)`.
     * - `right`: index of the group corresponding to the `right` cluster in each ScoreMarkersResults object.
     *    e.g., Cohen's d for the RNA markers of the `left` cluster are defined as `output.results.RNA.cohensD(output.right)`.
     */
    computeVersus(left, right) {
        var clusters = this.#choice.fetchClusters();
        var block = this.#filter.fetchFilteredBlock();

        // No need to free this afterwards; we don't own the normalized matrices anyway.
        let matrices = new scran.MultiMatrix;
        for (const [modality, state] of Object.entries(this.#norm_states)) {
            if (!state.valid()) {
                continue;
            }
            matrices.add(modality, state.fetchNormalizedMatrix());
        }

        if (!("versus" in this.#cache)) {
            this.#cache["versus"] = {};
        }

        return _computeVersus(left, right, matrices, clusters, { 
            cache: this.#cache.versus, 
            block: block,
            lfc_threshold: this.#parameters.lfc_threshold,
            compute_auc: this.#parameters.compute_auc
        });
    }
}

/**********************
 ***** Standalone *****
 **********************/

/**
 * Standalone version of {@linkplain MarkerDetectionState} that provides the same functionality outside of {@linkcode runAnalysis}.
 * Users can supply their own normalized matrices, groups and blocking factor to compute the various marker statistics for each group.
 * Users are also responsible for ensuring that the lifetime of the supplied objects exceeds that of the constructed MarkerDetectionStandalone instance,
 * i.e., the Wasm-related `free()` methods are not called while the MarkerDetectionStandalone instance is still in operation.
 */
export class MarkerDetectionStandalone {
    #matrices;

    #groups;
    #group_levels;
    #block;
    #block_levels;

    #cache;
    #parameters;

    /**
     * @param {external:MultiMatrix} normalized - A {@linkplain external:MultiMatrix MultiMatrix} of log-normalized values for multiple modalities.
     * @param {Array|TypedArray} groups - Array of length equal to the number of columns in any value of `normalized`, containing the group assignments for each column. 
     * @param {object} [options={}] - Optional parameters.
     * @param {?(Array|TypedArray)} [options.block=null] - Array of length equal to the number of columns in any value of `normalized`, containing the block assignments for each column.
     * If `null`, all columns are assigned to the same block.
     */
    constructor(normalized, groups, { block = null } = {}) {
        // Checking dimensions.
        let N = groups.length;

        for (const k of normalized.available()) {
            let v = normalized.get(k);
            if (v.numberOfColumns() != N) {
                throw new Error("all matrices in 'normalized' should have the same number of columns as the length of 'groups'");
            }
        }

        if (block !== null) {
            if (block.length != N) {
                throw new Error("'block' should have the same length as 'groups' if not null");
            }
        }

        let arrays = [ groups ];
        if (block !== null) {
            arrays.push(block);
        }
        let dump = utils.subsetInvalidFactors(arrays);

        // Maybe taking a subset to eliminate invalid entries.
        let new_matrices;
        if (dump.retain !== null) {
            new_matrices = new scran.MultiMatrix;
            let temp = scran.createInt32WasmArray(dump.retain.length);
            try {
                temp.set(dump.retain);
                for (const k of normalized.available()) {
                    new_matrices.add(k, scran.subsetColumns(normalized.get(k), temp))
                }
            } catch (e) {
                new_matrices.free();
                throw e;
            } finally {
                scran.free(temp);
            }
        } else {
            new_matrices = normalized.clone();
        }

        this.#matrices = new_matrices;
        this.#groups = dump.arrays[0].ids;
        this.#group_levels = dump.arrays[0].levels;

        if (block !== null) {
            this.#block = dump.arrays[1].ids;
            this.#block_levels = dump.arrays[1].levels;
        } else {
            this.#block = null;
            this.#block_levels = null;
        }

        this.#cache = { raw: {}, init: true };
        this.#parameters = MarkerDetectionState.defaults();
        this.changed = false;
    }

    /**
     * Frees all resources associated with this instance.
     */
    free() {
        scran.free(this.#groups);
        scran.free(this.#block);
        scran.free(this.#matrices);
        _free_results(this.#cache);
    }

    /**
     * @return {object} Marker detection results for the all modalities.
     * Each key is a modality name and each value is an {@linkplain external:ScoreMarkersResults ScoreMarkerResults} object,
     * containing marker detection statistics for all groups.
     * This is available after running {@linkcode MarkerDetectionState#compute compute}.
     */
    fetchResults() {
        return this.#cache.raw;
    }

    /**
     * @return {object} Object containing the parameters.
     */
    fetchParameters() {
        return { ...this.#parameters }; // avoid pass-by-reference links.
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.copy=true] - Whether to copy the return value on output.
     * Set to `false` for greater efficiency in strictly read-only applications.
     *
     * @return {Array} Array of levels for the grouping factor.
     * Group indices in the {@linkplain external:ScoreMarkersResults ScoreMarkersResults} instances returned by {@linkcode fetchResults} can be cross-referenced to this array.
     */
    fetchGroupLevels({ copy = true } = {}) {
        let ret = this.#group_levels;
        return (copy ? ret.slice() : ret);
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
        return this.#matrices;
    }

    _peekGroups() {
        return this.#groups;
    }

    _peekBlock() {
        return this.#block;
    }

    /**
     * If this method is not called, the parameters default to those in {@linkcode MarkerDetectionState#defaults MarkerDetectionState.defaults}.
     *
     * @param {object} parameters - Parameter object, see the argument of the same name in {@linkcode MarkerDetectionState#compute MarkerDetectionState.compute} for more details.
     *
     * @return The state is updated with new parameters.
     */
    setParameters(parameters) {
        if (this.#parameters.lfc_threshold !== parameters.lfc_threshold || this.#parameters.compute_auc !== parameters.compute_auc) {
            // Removing existing results, as they are now invalid.
            _free_results(this.#cache);
        }
        this.#parameters = { ...parameters };
        return;
    }

    /**
     * Compute markers for all groups.
     *
     * @return The state is updated with new results.
     */
    computeAll() {
        let { lfc_threshold, compute_auc } = this.#parameters;
        for (const k of this.#matrices.available()) {
            var mat = this.#matrices.get(k);
            utils.freeCache(this.#cache.raw[k]);
            this.#cache.raw[k] = scran.scoreMarkers(mat, this.#groups, { block: this.#block, threshold: lfc_threshold, computeAuc: compute_auc });
        }
        return;
    }

    /**
     * Extract markers for a pairwise comparison between two groups, 
     * for more detailed examination of the differences between them.
     *
     * @param {number} left - Index of one group in which to find upregulated markers.
     * @param {number} right - Index of another group to be compared against `left`.
     *
     * @return {object} Object containing:
     *
     * - `results`: object containing the marker statistics for the comparison between two groups.
     *    Each key is a modality name and each value is a {@linkplain external:ScoreMarkersResults ScoreMarkersResults} object.
     * - `left`: index of the group corresponding to the `left` group in each ScoreMarkersResults object,
     *    e.g., Cohen's d for the RNA markers of the `left` group are defined as `output.results.RNA.cohensD(output.left)`.
     * - `right`: index of the group corresponding to the `right` group in each ScoreMarkersResults object.
     *    e.g., Cohen's d for the RNA markers of the `left` group are defined as `output.results.RNA.cohensD(output.right)`.
     */
    computeVersus(left, right) {
        if (!("versus" in this.#cache)) {
            this.#cache["versus"] = {};
        }
        return _computeVersus(left, right, this.#matrices, this.#groups, { 
            cache: this.#cache.versus, 
            block: this.#block,
            lfc_threshold: this.#parameters.lfc_threshold,
            compute_auc: this.#parameters.compute_auc
        });
    }
}
