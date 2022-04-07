import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as markers from "./utils/markers.js";
import * as qc_module from "./quality_control.js";
import * as norm_module from "./normalization.js";
import * as choice_module from "./choose_clustering.js";

/**
 * This step performs marker detection for each cluster of cells by performing pairwise comparisons to each other cluster.
 * This wraps the `scoreMarkers` function from [**scran.js**](https://github.com/jkanche/scran.js).
 * The clustering is obtained from the {@linkcode choose_clustering} step.
 *
 * The parameters in {@linkcode runAnalysis} should be an empty object.
 *
 * Calling the **`results()`** method for the relevant state instance will return an empty object.
 *
 * The state instance for this step has a **`numberOfGroups()`** method.
 * This returns the number of clusters for which markers were computed.
 *
 * The state instance for this step has a **`fetchGroupResults(rank_type, group)`** method.
 * 
 * - `group` should be an integer specifying the cluster of interest.
 *   This should be less than the value returned by `numberOfGroups()`. 
 * - `rank_type` should be a string specifying the effect size to use for ranking markers.
 *   This should follow the format of `<effect>-<summary>` where `<effect>` may be `lfc`, `cohen`, `auc` or `delta_detected`,
 *   and `<summary>` may be `min`, `mean` or `min-rank`.
 * - An object is returned containing the marker statistics for the selection, sorted by the specified effect and summary size from `rank_type`.
 *   This contains:
 *   - `means`: a `Float64Array` of length equal to the number of genes, containing the mean expression within the selection.
 *   - `detected`: a `Float64Array` of length equal to the number of genes, containing the proportion of cells with detected expression inside the selection.
 *   - `lfc`: a `Float64Array` of length equal to the number of genes, containing the log-fold changes for the comparison between cells inside and outside the selection.
 *   - `delta_detected`: a `Float64Array` of length equal to the number of genes, containing the difference in the detected proportions between cells inside and outside the selection.
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 *
 * @namespace marker_detection
 */

export class State {
    #qc;
    #norm;
    #choice;
    #parameters;
    #cache;

    constructor(qc, norm, choice, parameters = null, cache = null) {
        if (!(qc instanceof qc_module.State)) {
            throw new Error("'qc' should be a State object from './quality_control.js'");
        }
        this.#qc = qc;

        if (!(norm instanceof norm_module.State)) {
            throw new Error("'norm' should be a State object from './normalization.js'");
        }
        this.#norm = norm;

        if (!(choice instanceof choice_module.State)) {
            throw new Error("'choice' should be a State object from './choose_clustering.js'");
        }
        this.#choice = choice;

        this.#parameters = (parameters === null ? {} : parameters);
        this.#cache = (cache === null ? {} : cache);
        this.changed = false;
    }

    free() {
        utils.freeCache(this.#cache.raw);
    }

    /***************************
     ******** Getters **********
     ***************************/

    fetchGroupResults(group, rank_type) {
        return markers.fetchGroupResults(this.#cache.raw, group, rank_type); 
    }

    numberOfGroups() {
        return this.#cache.raw.numberOfGroups();
    }

    fetchGroupMeans(group, { copy = true }) {
        return this.#cache.raw.means(group, { copy: copy });
    }

    /***************************
     ******** Compute **********
     ***************************/

    compute() {
        this.changed = false;

        if (this.#norm.changed || this.#choice.changed) {
            var mat = this.#norm.fetchNormalizedMatrix();
            var clusters = this.#choice.fetchClustersAsWasmArray();
            var block = this.#qc.fetchFilteredBlock();
            
            utils.freeCache(this.#cache.raw);
            this.#cache.raw = scran.scoreMarkers(mat, clusters, { block: block });

            // No parameters to set.

            this.changed = true;
        }

        return;
    }

    /***************************
     ******** Results **********
     ***************************/

    results() {
        return {};
    }

    /*************************
     ******** Saving *********
     *************************/

    serialize(handle) {
        let ghandle = handle.createGroup("marker_detection");
        ghandle.createGroup("parameters");

        {
            let chandle = ghandle.createGroup("results");
            let rhandle = chandle.createGroup("clusters");

            var num = this.#cache.raw.numberOfGroups();
            for (var i = 0; i < num; i++) {
                let ihandle = rhandle.createGroup(String(i));
                markers.serializeGroupStats(ihandle, this.#cache.raw, i);
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

export function unserialize(handle, permuter, qc, norm, choice) {
    let ghandle = handle.open("marker_detection");

    // No parameters to unserialize.
    let parameters = {};

    let cache = {};
    {
        let chandle = ghandle.open("results");
        let rhandle = chandle.open("clusters");
        let clusters = {};
        for (const cl of Object.keys(rhandle.children)) {
            clusters[Number(cl)] = markers.unserializeGroupStats(rhandle.open(cl), permuter);
        }
        cache.raw = new ScoreMarkersMimic(clusters);
    }

    return {
        state: new State(qc, norm, choice, parameters, cache),
        parameters: { ...parameters }
    };
}

