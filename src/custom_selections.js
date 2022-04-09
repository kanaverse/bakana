import * as scran from "scran.js";
import * as utils from "./utils/general.js";
import * as markers from "./utils/markers.js";
import * as qc_module from "./quality_control.js";
import * as norm_module from "./normalization.js";

/**
 * Applications can perform marker detection on custom selections of cells.
 * This allows users to dynamically select cells on a UI and quickly obtain a list of distinguishing markers for that selection.
 * This wraps the `scoreMarkers` function from [**scran.js**](https://github.com/jkanche/scran.js).
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class CustomSelectionsState {
    #qc;
    #norm;
    #cache;
    #parameters;

    constructor(qc, norm, parameters = null, cache = null) {
        if (!(qc instanceof qc_module.QualityControlState)) {
            throw new Error("'qc' should be a State object from './quality_control.js'");
        }
        this.#qc = qc;

        if (!(norm instanceof norm_module.NormalizationState)) {
            throw new Error("'norm' should be a State object from './normalization.js'");
        }
        this.#norm = norm;

        this.#cache = (cache === null ? { "results": {} } : cache); 
        this.#parameters = (parameters === null ? { "selections": {} } : parameters);
        this.changed = false;
    }

    free() {
        utils.freeCache(this.#cache.buffer);
        for (const [k, v] of Object.entries(this.#cache.results)) {
            v.raw.free();
        }
    }

    /***************************
     ******** Setters **********
     ***************************/

    /**
     * Add a custom selection and compute its markers.
     *
     * @param {string} id A unique identifier for the new custom selection.
     * @param {Array|TypedArray} selection The indices of the cells in the selection.
     * Indices should refer to positions of cells in the QC-filtered matrix, not the original matrix.
     *
     * @return The custom selection is added to the state and calculation of its markers is performed.
     * Nothing is returned.
     */
    addSelection(id, selection) {
        var mat = this.#norm.fetchNormalizedMatrix();

        var buffer = utils.allocateCachedArray(mat.numberOfColumns(), "Int32Array", this.#cache);
        buffer.fill(0);
        var tmp = buffer.array();
        selection.forEach(element => { tmp[element] = 1; });

        // Assumes that we have at least one cell in and outside the selection!
        var res = scran.scoreMarkers(mat, buffer); 
      
        // Removing previous results, if there were any.
        if (id in this.#cache.results) {
            utils.freeCache(this.#cache.results[id].raw);
        }
      
        this.#cache.results[id] = { "raw": res };
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
        utils.freeCache(this.#cache.results[id].raw);
        delete this.#cache.results[id];
        delete this.#parameters.selections[id];
        return;
    }

    /***************************
     ******** Getters **********
     ***************************/

    /**
     * Fetch the marker results for a custom selection.
     * 
     * @param {string} id - An identifier for the desired selection.
     * @param {string} rank_type - Effect size to use for ranking markers.
     * This should be one of `lfc`, `cohen`, `auc` or `delta_detected`.
     *
     * @return An object containing the marker statistics for the selection, sorted by the specified effect and summary size from `rank_type`.
     * This contains:
     * - `means`: a `Float64Array` of length equal to the number of genes, containing the mean expression within the selection.
     * - `detected`: a `Float64Array` of length equal to the number of genes, containing the proportion of cells with detected expression inside the selection.
     * - `lfc`: a `Float64Array` of length equal to the number of genes, containing the log-fold changes for the comparison between cells inside and outside the selection.
     * - `delta_detected`: a `Float64Array` of length equal to the number of genes, containing the difference in the detected proportions between cells inside and outside the selection.
     */
    fetchResults(id, rank_type) {
        var current = this.#cache.results[id].raw;
        return markers.fetchGroupResults(current, 1, rank_type + "-mean"); 
    }

    /***************************
     ******** Compute **********
     ***************************/

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     *
     * @return The state is updated by removing stale selections if the QC filter was altered.
     */
    compute() {
        this.changed = false;

        /* If the QC filter was re-run, all of the selections are invalidated as
         * the identity of the indices may have changed.
         */
        if (this.#qc.changed) {
            for (const [key, val] of Object.entries(this.#cache.results)) {
                utils.freeCache(val.raw);                    
            }
            this.#parameters.selections = {};
            this.#cache.results = {};
            this.changed = true;
        }

        /*
         * Technically we would need to re-run detection on the existing selections
         * if the normalization changed but the QC was the same. In practice, this
         * never happens, so we'll deal with it later.
         */

        return;
    }

    /***************************
     ******** Results **********
     **************************/

    /**
     * Obtain a summary of the state, typically for display on a UI like **kana**.
     *
     * @return An empty object.
     * This is returned for consistency with the other steps.
     */
    summary() {
        return this.#parameters;
    }

    /*************************
     ******** Saving *********
     *************************/

    serialize(handle) {
        let ghandle = handle.createGroup("custom_selections");

        {
            let phandle = ghandle.createGroup("parameters");
            let rhandle = phandle.createGroup("selections");
            for (const [key, val] of Object.entries(this.#parameters.selections)) {
                rhandle.writeDataSet(String(key), "Int32", null, val);
            }
        }

        {
            let chandle = ghandle.createGroup("results");
            let rhandle = chandle.createGroup("markers");
            for (const [key, val] of Object.entries(this.#cache.results)) {
                let ihandle = rhandle.createGroup(key);
                markers.serializeGroupStats(ihandle, val.raw, 1, { no_summaries: true });
            }
        }
    }
}

/**************************
 ******** Loading *********
 **************************/

class CustomMarkersMimic {
    constructor(results) {
        this.results = results;
    }

    effect_grabber(key, group, summary, copy) {
        if (group != 1) {
            throw "only group 1 is supported for custom marker mimics";
        }
        if (summary != 1) {
            throw "only the mean effect size is supported for custom marker mimics";
        }
        let chosen = this.results[group][key];
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
        let chosen = this.results[group][key];
        return utils.mimicGetter(chosen, copy);
    }

    means(group, { copy }) {
        return this.stat_grabber("means", group, copy);
    }

    detected(group, { copy }) {
        return this.stat_grabber("detected", group, copy);
    }

    free() {}
}

export function unserialize(handle, permuter, qc, norm) {
    let ghandle = handle.open("custom_selections");

    let parameters = { selections: {} };
    {
        let phandle = ghandle.open("parameters");
        let rhandle = phandle.open("selections");

        for (const key of Object.keys(rhandle.children)) {
            parameters.selections[key] = rhandle.open(key, { load: true }).values;
        }
    }

    let cache = { results: {} };
    {
        let chandle = ghandle.open("results");
        let rhandle = chandle.open("markers");

        for (const sel of Object.keys(rhandle.children)) {
            let current = markers.unserializeGroupStats(rhandle.open(sel), permuter, { no_summaries: true });
            cache.results[sel] = { raw: new CustomMarkersMimic({ 1 : current }) };
        }
    }

    // Need to make a copy to avoid moving the buffers.
    let output = { selections: {} };
    for (const [k, v] of Object.entries(parameters.selections)) {
        output.selections[k] = v.slice();        
    }

    return {
        state: new CustomSelectionsState(qc, norm, parameters, cache),
        parameters: output
    };
}

