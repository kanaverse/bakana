import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as neighbor_module from "./neighbor_index.js";

export const step_name = "snn_graph_cluster";

/**
 * This step does SNN graph clustering based on the neighbor search index built by {@linkplain NeighborIndexState}.
 * This wraps [`clusterSNNGraph`](https://jkanche.com/scran.js/global.html#clusterSNNGraph) 
 * and related functions from [**scran.js**](https://github.com/jkanche/scran.js).
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class SnnGraphClusterState {
    #index;
    #parameters;
    #cache;

    constructor(index, parameters = null, cache = null) {
        if (!(index instanceof neighbor_module.NeighborIndexState)) {
            throw new Error("'index' should be a State object from './neighbor_index.js'");
        }
        this.#index = index;

        this.#parameters = (parameters === null ? {} : parameters);
        this.#cache = (cache === null ? {} : cache);
        this.changed = false;
    }

    free() {
        utils.freeCache(this.#cache.neighbors);
        utils.freeCache(this.#cache.graph);
        utils.freeCache(this.#cache.clusters);
    }

    /***************************
     ******** Getters **********
     ***************************/

    /**
     * @return {Int32WasmArray} Array of cluster assignments for each cell in the (filtered) dataset,
     * available after running {@linkcode SnnGraphClusterState#compute compute}.
     */
    fetchClusters() {
        if (!this.#valid()) {
            throw "cannot fetch SNN clusters from an invalid state";
        } else {
            return this.#cache.clusters.membership({ copy: "view" });
        }
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

    #valid() {
        return "clusters" in this.#cache;
    }

    #compute_neighbors(k) {
        this.#cache.neighbors = scran.findNearestNeighbors(this.#index.fetchIndex(), k);
        return;
    }

    #compute_graph(scheme) {
        if (!("neighbors" in this.#cache)) { // need to check as reloaded state will not populate the internals.
            this.#compute_neighbors(this.#parameters.k);
        }
        this.#cache.graph = scran.buildSNNGraph(this.#cache.neighbors, { scheme: scheme });
        return;
    }

    #compute_clusters(algorithm, multilevel_resolution, leiden_resolution, walktrap_steps) {
        if (!("graph" in this.#cache)) {
            this.#compute_graph(this.#parameters.scheme);
        }
        this.#cache.clusters = scran.clusterSNNGraph(this.#cache.graph, {
            method: algorithm,
            multiLevelResolution: multilevel_resolution,
            leidenResolution: leiden_resolution,
            leidenModularityObjective: true, // avoid problems with unstable interpretation of leidenResolution.
            walktrapSteps: walktrap_steps
        });
        return;
    }

    static defaults() {
        return { 
            k: 10,
            scheme: "rank",
            algorithm: "multilevel",
            multilevel_resolution: 1,
            leiden_resolution: 1,
            walktrap_steps: 4
        };
    }

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     * Each argument is taken from the property of the same name in the `snn_graph_cluster` property of the `parameters` of {@linkcode runAnalysis}.
     * The exception is `run_me`, which is computed internally and does not correspond to any parameter in `snn_graph_cluster`.
     *
     * @param {boolean} run_me - Whether or not to run this step, depending on the clustering method chosen by the user (see {@linkplain ChooseClusteringState}).
     * @param {number} k - Number of nearest neighbors used to construct the graph.
     * @param {string} scheme - Weighting scheme for graph construction.
     * This can be one of `"rank"`, `"number"` or `"jaccard"`.
     * @param {string} algorithm - Algorithm to use for community detection.
     * This can be one of `"multilevel"`, `"walktrap"` or `"leiden"`.
     * @param {number} multilevel_resolution - Resolution of the multi-level community detection.
     * @param {number} leiden_resolution - Resolution of the Leiden community detection.
     * @param {number} walktrap_steps - Number of merge steps for the Walktrap algorithm.
     *
     * @return The object is updated with the new results.
     */
    compute(run_me, k, scheme, algorithm, multilevel_resolution, leiden_resolution, walktrap_steps) {
        this.changed = false;

        if (this.#index.changed || k !== this.#parameters.k) {
            utils.freeCache(this.#cache.neighbors);
            if (run_me) {
                this.#compute_neighbors(k);
            } else {
                delete this.#cache.neighbors; // ensuring that this is re-run on future calls to compute() with run_me = true.
            }
            this.#parameters.k = k;
            this.changed = true;
        }

        if (this.changed || scheme !== this.#parameters.scheme) {
            utils.freeCache(this.#cache.graph);
            if (run_me) {
                this.#compute_graph(scheme);
            } else {
                delete this.#cache.graph;
            }
            this.#parameters.scheme = scheme;
            this.changed = true 
        }

        if (this.changed 
            || algorithm !== this.#parameters.algorithm 
            || multilevel_resolution !== this.#parameters.multilevel_resolution 
            || leiden_resolution !== this.#parameters.leiden_resolution 
            || walktrap_steps !== this.#parameters.walktrap_steps 
            || (!this.#valid() && run_me))
        {
            utils.freeCache(this.#cache.clusters);
            if (run_me) {
                this.#compute_clusters(algorithm, multilevel_resolution, leiden_resolution, walktrap_steps);
            } else {
                delete this.#cache.clusters;
            }

            this.#parameters.algorithm = algorithm;
            this.#parameters.multilevel_resolution = multilevel_resolution;
            this.#parameters.leiden_resolution = leiden_resolution;
            this.#parameters.walktrap_steps = walktrap_steps;
            this.changed = true;
        }

        return;
    }
}

/**************************
 ******** Loading *********
 **************************/

export function unserialize(handle, index) {
    let ghandle = handle.open("snn_graph_cluster");

    let parameters = SnnGraphClusterState.defaults();
    {
        let phandle = ghandle.open("parameters");
        parameters.k = phandle.open("k", { load: true }).values[0];

        parameters.scheme = phandle.open("scheme", { load: true }).values[0];
        if (typeof parameters.scheme !== "string") { // because I stuffed up and tried to save a string as an int in v1.0, oops.
            parameters.scheme = "rank";
        }

        if ("algorithm" in phandle.children) {
            // v3.0
            parameters.algorithm = phandle.open("algorithm", { load: true }).values[0];
            parameters.multilevel_resolution = phandle.open("multilevel_resolution", { load: true }).values[0];
            parameters.leiden_resolution = phandle.open("leiden_resolution", { load: true }).values[0];
            parameters.walktrap_steps = phandle.open("walktrap_steps", { load: true }).values[0];
        } else {
            // v2.0
            parameters.multilevel_resolution = phandle.open("resolution", { load: true }).values[0];
        }
    }

    let cache = {};
    {
        let rhandle = ghandle.open("results");
        if ("clusters" in rhandle.children) {
            let clusters = rhandle.open("clusters", { load: true }).values;
            cache.clusters = scran.emptyClusterSNNGraphResults(clusters.length, 1);
            cache.clusters.setBest(0); // whatever.
            cache.clusters.membership({ fillable: true }).set(clusters);
        }
    }

    return new SnnGraphClusterState(index, parameters, cache);
}


