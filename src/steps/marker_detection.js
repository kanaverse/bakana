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
    }

    /***************************
     ******** Getters **********
     ***************************/

    /**
     * Fetch the marker statistics for a given cluster.
     *
     * @param {number} group - An integer specifying the cluster of interest.
     * This should be less than the value returned by {@linkcode MarkerDetectionState#numberOfGroups numberOfGroups}.
     * @param {string} rank_type - Summarized effect size to use for ranking markers.
     * This should follow the format of `<effect>-<summary>` where `<effect>` may be `lfc`, `cohen`, `auc` or `delta_detected`,
     * and `<summary>` may be `min`, `mean` or `min-rank`.
     * @param {string} feat_type - The feature type of interest, usually `"RNA"` or `"ADT"`.
     *
     * @return An object containing the marker statistics for the selection, sorted by the specified effect and summary size from `rank_type`.
     * This contains:
     *   - `means`: a `Float64Array` of length equal to the number of genes, containing the mean expression within the selection.
     *   - `detected`: a `Float64Array` of length equal to the number of genes, containing the proportion of cells with detected expression inside the selection.
     *   - `lfc`: a `Float64Array` of length equal to the number of genes, containing the log-fold changes for the comparison between cells inside and outside the selection.
     *   - `delta_detected`: a `Float64Array` of length equal to the number of genes, containing the difference in the detected proportions between cells inside and outside the selection.
     */
    fetchGroupResults(group, rank_type, feat_type) {
        return markers.fetchGroupResults(this.#cache.raw[feat_type], group, rank_type); 
    }

    /** 
     * @return The number of clusters for which markers were computed.
     */
    numberOfGroups() {
        let first = Object.values(this.#cache.raw)[0];
        return first.numberOfGroups();
    }

    fetchGroupMeans(group, feat_type, { copy = true }) {
        return this.#cache.raw[feat_type].means(group, { copy: copy });
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
                var clusters = this.#choice.fetchClustersAsWasmArray();
                var block = this.#filter.fetchFilteredBlock();
                
                utils.freeCache(this.#cache.raw[k]);
                this.#cache.raw[k] = scran.scoreMarkers(mat, clusters, { block: block });

                // No parameters to set.

                this.changed = true;
            }
        }

        return;
    }

    static defaults() {
        return {};
    }

    /***************************
     ******** Results **********
     ***************************/

    /**
     * Obtain a summary of the state, typically for display on a UI like **kana**.
     *
     * @return An empty object.
     * This is returned for consistency with the other steps.
     */
    summary() {
        return {};
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

    lfc(group, { summary, copy }) {
        return this.effect_grabber("lfc", group, summary, copy);
    }

    deltaDetected(group, { summary, copy }) {
        return this.effect_grabber("delta_detected", group, summary, copy);
    }

    cohen(group, { summary, copy }) {
        return this.effect_grabber("cohen", group, summary, copy);
    }

    auc(group, { summary, copy }) {
        return this.effect_grabber("auc", group, summary, copy);
    }

    stat_grabber(key, group, copy) {
        let chosen = this.clusters[group][key];
        return utils.mimicGetter(chosen, copy);
    }

    means(group, { copy }) {
        return this.stat_grabber("means", group, copy);
    }

    detected(group, { copy }) {
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

    return {
        state: new MarkerDetectionState(filter, norm_states, choice, parameters, cache),
        parameters: { ...parameters } // make a copy to avoid pass-by-reference links with state's internal parameters
    };
}

