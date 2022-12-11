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
        markers.freeVersusResults(this.#cache);
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
     * Compute markers between two clusters as described for {@linkcode MarkerDetectionState#computeVersus computeVersus}, but for a custom matrix, clusters and cache.
     * This allows applications to run versus-mode comparisons on custom inputs without creating a MarkerDetectionState object.
     *
     * @param {number} left - Index of one cluster in which to find upregulated markers.
     * @param {number} right - Index of another cluster to be compared against `left`.
     * @param {string} rank_type - Effect size to use for ranking markers.
     * @param {string} feat_type - Feature type of interest.
     * @param {ScranMatrix} mat - Matrix containing log-expression data for all cells.
     * @param {Array|TypedArray|Int32WasmArray} clusters - Integer array of length containing the cluster assignments for all cells in `mat`.
     * @param {object} [options={}] - Optional parameters.
     * @param {object} [options.cache={}] - Cache for the results.
     * @param {?Int32WasmArray} [options.block=null] - Blocking factor of length equal to the number of cells in `mat`, 
     * see {@linkcode CellFilteringState#fetchFilteredBlock CellFilteringState.fetchFilteredBlock}.
     *
     * @return An object containing the marker statistics, see {@linkcode MarkerDetectionState#computeVersus computeVersus} for more details.
     */
    static computeVersusCustom(left, right, mat, clusters, { cache = {}, block = null } = {}) {

        let generate = (smal, bigg) => {
            let new_clusters = [];
            let keep = [];
            let smalfound = false, biggfound = false;
            clusters.forEach((x, i) => {
                if (x == smal) {
                    new_clusters.push(0);
                    keep.push(i);
                    smalfound = true;
                } else if (x == bigg) {
                    new_clusters.push(1);
                    keep.push(i);
                    biggfound = true;
                }
            });

            if (!smalfound || !biggfound) {
                throw new Error("non-zero entries should be present for both requested clusters in versus mode");
            }

            let new_block = null;
            if (block !== null) {
                let block_arr = block.array();
                new_block = keep.map(i => block_arr[i]);
                markers.dropUnusedBlocks(new_block);
            }

            let output;
            let sub;
            try {
                sub = scran.subsetColumns(mat, keep);
                output = scran.scoreMarkers(sub, new_clusters, { block: new_block });
            } finally {
                utils.freeCache(sub);
            }

            return output;
        };

        return markers.generateVersusResults(left, right, rank_type, feat_type, cache, generate);
    }

    /**
     * Extract markers for a pairwise comparison between two clusters, for more detailed examination of the differences between them.
     *
     * @param {number} left - Index of one cluster in which to find upregulated markers.
     * @param {number} right - Index of another cluster to be compared against `left`.
     * Effect sizes are defined as `left - right`.
     * @param {string} rank_type - Effect size to use for ranking markers.
     * This should be one of `lfc`, `cohen`, `auc` or `delta_detected`.
     *
     * @param {string} feat_type - The feature type of interest, usually `"RNA"` or `"ADT"`.
     *
     * @return An object containing the marker statistics for this pairwise comparison, sorted by the specified effect from `rank_type`.
     * This has the same structure as described for {@linkcode MarkerDetectionState#fetchGroupResults fetchGroupResults}.
     */
    computeVersus(left, right, rank_type, feat_type) {
        var clusters = this.#choice.fetchClusters();
        var block = this.#filter.fetchFilteredBlock();

        for (const [feat_type, mat] of Object.entries(this.#norm_states)) {
            return this.constructor.computeVersusCustom(left, right, rank_type, feat_type, mat, clusters, { cache: this.#cache, block: block });
        }
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

