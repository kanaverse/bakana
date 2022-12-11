import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as nutils from "./utils/normalization.js";
import * as markers from "./utils/markers.js";
import * as filter_module from "./cell_filtering.js";
import * as choice_module from "./choose_clustering.js";

export const step_name = "marker_detection";

/**
 * This step performs marker detection for each cluster of cells by performing pairwise comparisons to each other cluster.
 * This wraps the `scoreMarkers` function from [**scran.js**](https://github.com/jkanche/scran.js).
 * The clustering is obtained from the {@linkcode choose_clustering} step.
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

        for (const norm of Object.values(norm_states)) {
            if (!(norm instanceof nutils.NormalizationStateBase)) {
                throw new Error("'norm' should be a NormalizationStateBase object");
            }
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

    free() {
        for (const v of Object.values(this.#cache.raw)) {
            utils.freeCache(v);
        }
        markers.freeVersusResults(this.#cache.versus);
    }

    /***************************
     ******** Getters **********
     ***************************/

    /**
     * @return {object} Marker detection results for the all modalities.
     * Each key is a modality name and each value is a ScoreMarkersResults object,
     * containing marker detection statistics for all clusters.
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

    /***************************
     ******** Compute **********
     ***************************/

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     *
     * @return The state is updated with new results.
     */
    compute() {
        this.changed = false;

        for (const [k, v] of Object.entries(this.#norm_states)) {
            if (!v.valid()) {
                continue;
            }

            if (v.changed || this.#choice.changed) {
                var mat = v.fetchNormalizedMatrix();
                var clusters = this.#choice.fetchClusters();
                var block = this.#filter.fetchFilteredBlock();
                
                utils.freeCache(this.#cache.raw[k]);
                this.#cache.raw[k] = scran.scoreMarkers(mat, clusters, { block: block });

                // No parameters to set.

                this.changed = true;
            }
        }

        if (this.changed) {
            markers.freeVersusResults(this.#cache);
        }

        return;
    }

    static defaults() {
        return {};
    }

    /*******************************
     ******** Versus mode **********
     *******************************/

    /**
     * Compute markers between two clusters as described for {@linkcode MarkerDetectionState#computeVersus computeVersus}.
     * This allows applications to run versus-mode comparisons on custom inputs without creating a MarkerDetectionState object.
     *
     * @param {number} left - Index of one cluster. 
     * @param {number} right - Index of another cluster to be compared against `left`.
     * @param {MultiMatrix} matrices - A MultiMatrix object containing log-normalized matrices for each modality.
     * @param {Array|TypedArray|Int32WasmArray} clusters - Integer array of length equal to the number of cells in `matrices`,
     * containing the cluster assignments for all cells in `matrices`.
     * @param {object} [options={}] - Optional parameters.
     * @param {object} [options.cache={}] - Cache for the results.
     * @param {?Int32WasmArray} [options.block=null] - Blocking factor of length equal to the number of cells in `mat`, 
     * see {@linkcode CellFilteringState#fetchFilteredBlock CellFilteringState.fetchFilteredBlock}.
     *
     * @return {object} Object containing:
     *
     * - `results`: object containing the marker statistics for the comparison between two clusters.
     *    Each key is a modality name and each value is a ScoreMarkersResults object.
     * - `left`: index of the group corresponding to the `left` cluster in each ScoreMarkersResults object.
     *    e.g., Cohen's d for the RNA markers of the `left` cluster are defined as `output.results.RNA.cohen(output.left)`.
     * - `right`: index of the group corresponding to the `right` cluster in each ScoreMarkersResults object.
     *    e.g., Cohen's d for the RNA markers of the `right` cluster are defined as `output.results.RNA.cohen(output.right)`.
     */
    static computeVersusCustom(left, right, matrices, clusters, { cache = {}, block = null } = {}) {
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

            markers.computeVersusResults(matrices, new_clusters, block, keep, cache_info.cached);
        }

        return { 
            results: cache_info.cached,
            left: left_index,
            right: right_index
        };
    }

    /**
     * Extract markers for a pairwise comparison between two clusters, 
     * for more detailed examination of the differences between them.
     *
     * @param {number} left - Index of one cluster in which to find upregulated markers.
     * @param {number} right - Index of another cluster to be compared against `left`.
     *
     * @return {object} Object containing:
     *
     * - `results`: object containing the marker statistics for the comparison between two clusters.
     *    Each key is a modality name and each value is a ScoreMarkersResults object.
     * - `left`: index of the group corresponding to the `left` cluster in each ScoreMarkersResults object,
     *    e.g., Cohen's d for the RNA markers of the `left` cluster are defined as `output.results.RNA.cohen(output.left)`.
     * - `right`: index of the group corresponding to the `right` cluster in each ScoreMarkersResults object.
     *    e.g., Cohen's d for the RNA markers of the `left` cluster are defined as `output.results.RNA.cohen(output.right)`.
     */
    computeVersus(left, right) {
        var clusters = this.#choice.fetchClusters();
        var block = this.#filter.fetchFilteredBlock();

        // No need to free this afterwards; we don't own the normalized matrices anyway.
        let matrices = new scran.MultiMatrix;
        for (const [modality, state] of Object.entries(this.#norm_states)) {
            let curmat = state.fetchNormalizedMatrix();
            if (curmat !== null) {
                matrices.add(modality, curmat);
            }
        }

        if (!("versus" in this.#cache)) {
            this.#cache["versus"] = {};
        }

        return this.constructor.computeVersusCustom(left, right, matrices, clusters, { 
            cache: this.#cache.versus, 
            block: block
        });
    }

    /*************************
     ******** Saving *********
     *************************/

    serialize(handle) {
        let ghandle = handle.createGroup("marker_detection");
        ghandle.createGroup("parameters");

        {
            let rhandle = ghandle.createGroup("results");
            let chandle = rhandle.createGroup("per_cluster");

            for (const [k, v] of Object.entries(this.#cache.raw)) {
                let khandle = chandle.createGroup(k);
                var num = v.numberOfGroups();
                for (var i = 0; i < num; i++) {
                    let ihandle = khandle.createGroup(String(i));
                    markers.serializeGroupStats(ihandle, v, i);
                }
            }
        }
    }
}

/**************************
 ******** Loading *********
 **************************/

class ScoreMarkersMimic {
    constructor(clusters) {
        this.clusters = clusters;
    }

    effect_grabber(key, group, summary, copy) {
        let sidx = markers.int2summaries[summary];
        let chosen = this.clusters[group][key][sidx];
        return utils.mimicGetter(chosen, copy);
    }

    lfc(group, { summary = 1, copy = true } = {}) {
        return this.effect_grabber("lfc", group, summary, copy);
    }

    deltaDetected(group, { summary = 1, copy = true } = {}) {
        return this.effect_grabber("delta_detected", group, summary, copy);
    }

    cohen(group, { summary = 1, copy = true } = {}) {
        return this.effect_grabber("cohen", group, summary, copy);
    }

    auc(group, { summary = 1, copy = true } = {}) {
        return this.effect_grabber("auc", group, summary, copy);
    }

    stat_grabber(key, group, copy) {
        let chosen = this.clusters[group][key];
        return utils.mimicGetter(chosen, copy);
    }

    means(group, { copy = true } = {}) {
        return this.stat_grabber("means", group, copy);
    }

    detected(group, { copy = true } = {}) {
        return this.stat_grabber("detected", group, copy);
    }

    numberOfGroups() {
        return Object.keys(this.clusters).length;
    }

    free() {}
}

export function unserialize(handle, permuters, filter, norm_states, choice) {
    let ghandle = handle.open("marker_detection");

    // No parameters to unserialize.
    let parameters = {};

    let cache = {};
    {
        let rhandle = ghandle.open("results");

        if ("clusters" in rhandle.children) { 
            // below v2.0
            let chandle = rhandle.open("clusters");
            let clusters = {};
            for (const cl of Object.keys(chandle.children)) {
                clusters[Number(cl)] = markers.unserializeGroupStats(chandle.open(cl), permuters["RNA"]);
            }
            cache.raw = { RNA: new ScoreMarkersMimic(clusters) };
        } else {
            // after v2.0.
            let chandle = rhandle.open("per_cluster");
            cache.raw = {};
            for (const a of Object.keys(chandle.children)) {
                let clusters = {};
                let ahandle = chandle.open(a);
                for (const cl of Object.keys(ahandle.children)) {
                    clusters[Number(cl)] = markers.unserializeGroupStats(ahandle.open(cl), permuters[a]);
                }
                cache.raw[a] = new ScoreMarkersMimic(clusters);
            }
        }
    }

    return new MarkerDetectionState(filter, norm_states, choice, parameters, cache);
}

