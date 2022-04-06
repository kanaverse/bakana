import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as markers from "./utils/markers.js";
import * as qc_module from "./quality_control.js";
import * as norm_module from "./normalization.js";
import * as choice_module from "./choose_clustering.js";

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

    fetchGroupResults(rank_type, group) {
        return markers.fetchGroupResults(this.#cache.raw, rank_type, group); 
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
                markers.serializeGroupStats(rhandle, this.#cache.raw, i);
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

