import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as markers from "./utils/markers.js";
import * as filter_module from "./cell_filtering.js";
import * as choice_module from "./choose_clustering.js";
import * as rna_norm_module from "./rna_normalization.js";
import * as adt_norm_module from "./adt_normalization.js";

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

        if (!(norm_states.RNA instanceof rna_norm_module.RnaNormalizationState)) {
            throw new Error("'norm_states.RNA' should be an RnaNormalizationState object");
        }
        if (!(norm_states.ADT instanceof adt_norm_module.AdtNormalizationState)) {
            throw new Error("'norm_states.ADT' should be an AdtNormalizationState object");
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
     * Each argument is taken from the property of the same name in the `marker_detection` property of the `parameters` of {@linkcode runAnalysis}.
     *
     * @param {number} lfc_threshold - Log-fold change threshold to use when computing the Cohen's d and AUC for each pairwise comparison.
     * @param {boolean} compute_auc - Whether to compute the AUCs.
     * Setting this to `false` will skip AUC calculations and improve speed and memory efficiency.
     *
     * @return The state is updated with new results.
     */
    compute(lfc_threshold, compute_auc) {
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
                this.#cache.raw[k] = scran.scoreMarkers(mat, clusters, { block: block, lfcThreshold: lfc_threshold, computeAuc: compute_auc });

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
     * @param {number} [options.lfc_threshold=0] - Log-fold change threshold to use when computing the Cohen's d and AUC for each pairwise comparison.
     * @param {boolean} [options.compute_auc=true] - Whether to compute the AUCs.
     * Setting this to `false` will skip AUC calculations and improve speed and memory efficiency.
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
    static computeVersusCustom(left, right, matrices, clusters, { cache = {}, block = null, lfc_threshold = 0, compute_auc = true } = {}) {
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
            if (!state.valid()) {
                continue;
            }
            matrices.add(modality, state.fetchNormalizedMatrix());
        }

        if (!("versus" in this.#cache)) {
            this.#cache["versus"] = {};
        }

        return this.constructor.computeVersusCustom(left, right, matrices, clusters, { 
            cache: this.#cache.versus, 
            block: block,
            lfc_threshold: this.#parameters.lfc_threshold,
            compute_auc: this.#parameters.compute_auc
        });
    }

    /*************************
     ******** Saving *********
     *************************/

    serialize(handle) {
        let ghandle = handle.createGroup("marker_detection");

        {
            let phandle = ghandle.createGroup("parameters");
            phandle.writeDataSet("lfc_threshold", "Float64", [], this.#parameters.lfc_threshold)
            phandle.writeDataSet("compute_auc", "Uint8", [], (this.#parameters.compute_auc ? 1 : 0));
        }

        {
            let rhandle = ghandle.createGroup("results");
            let chandle = rhandle.createGroup("per_cluster");

            for (const [k, v] of Object.entries(this.#cache.raw)) {
                let khandle = chandle.createGroup(k);
                var num = v.numberOfGroups();
                for (var i = 0; i < num; i++) {
                    let ihandle = khandle.createGroup(String(i));
                    markers.serializeGroupStats(ihandle, v, i, { compute_auc: this.#parameters.compute_auc });
                }
            }
        }
    }
}

/**************************
 ******** Loading *********
 **************************/

function fill_results(stats, num_blocks) {
    let keys = Object.keys(stats);
    let first = stats[keys[0]];
    let ngenes = first.means.length;
    let object = scran.emptyScoreMarkersResults(ngenes, keys.length, num_blocks, { computeAuc: ("auc" in first) });

    for (const k of keys) {
        let i = Number(k);
        let vals = stats[k];
        object.means(i, { fillable: true }).set(vals.means);
        object.detected(i, { fillable: true }).set(vals.detected);

        for (const [s, v] of Object.entries(vals.cohen)) {
            object.cohen(i, { summary: markers.summaries2int[s], fillable: true }).set(v);
        }

        for (const [s, v] of Object.entries(vals.lfc)) {
            object.lfc(i, { summary: markers.summaries2int[s], fillable: true }).set(v);
        }

        for (const [s, v] of Object.entries(vals.delta_detected)) {
            object.deltaDetected(i, { summary: markers.summaries2int[s], fillable: true }).set(v);
        }

        if ("auc" in vals) {
            for (const [s, v] of Object.entries(vals.auc)) {
                object.auc(i, { summary: markers.summaries2int[s], fillable: true }).set(v);
            }
        }
    }

    return object;
}

export function unserialize(handle, permuters, filter, norm_states, choice) {
    let ghandle = handle.open("marker_detection");

    let parameters = MarkerDetectionState.defaults();
    {
        let phandle = ghandle.open("parameters");
        if ("lfc_threshold" in phandle.children) {
            parameters.lfc_threshold = phandle.open("lfc_threshold", { load: true }).values[0];
        }
        if ("compute_auc" in phandle.children) {
            parameters.compute_auc = phandle.open("compute_auc", { load: true }).values[0] > 0;
        }
    }

    // Figure out the number of blocks.
    let num_blocks = 1;
    {
        let filtered = filter.fetchFilteredBlock();
        if (filtered != null) {
            filtered.forEach(x => {
                if (x + 1 > num_blocks) {
                    num_blocks = x + 1;
                }
            });
        }
    }

    // Set up the marker detection statistics.
    let cache = {};
    {
        let rhandle = ghandle.open("results");
        cache.raw = {};

        if ("clusters" in rhandle.children) { 
            // below v2.0
            let chandle = rhandle.open("clusters");
            let clusters = {};
            for (const cl of Object.keys(chandle.children)) {
                clusters[Number(cl)] = markers.unserializeGroupStats(chandle.open(cl), permuters["RNA"], { compute_auc: parameters.compute_auc });
            }
            cache.raw.RNA = fill_results(clusters, num_blocks);
        } else {
            // after v2.0.
            let chandle = rhandle.open("per_cluster");
            for (const a of Object.keys(chandle.children)) {
                let clusters = {};
                let ahandle = chandle.open(a);
                for (const cl of Object.keys(ahandle.children)) {
                    clusters[Number(cl)] = markers.unserializeGroupStats(ahandle.open(cl), permuters[a], { compute_auc: parameters.compute_auc });
                }
                cache.raw[a] = fill_results(clusters, num_blocks);
            }
        }
 
    }

    return new MarkerDetectionState(filter, norm_states, choice, parameters, cache);
}

