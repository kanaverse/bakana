import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as neighbor_module from "./neighbor_index.js";

export const step_name = "snn_graph_cluster";

/**
 * This step does SNN graph clustering based on the neighbor search index built by {@linkplain NeighborIndexState}.
 * This wraps [`clusterGraph`](https://kanaverse.github.io/scran.js/global.html#clusterGraph) 
 * and related functions from [**scran.js**](https://github.com/kanaverse/scran.js).
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
        this.#cache.graph = scran.buildSnnGraph(this.#cache.neighbors, { scheme: scheme });
        return;
    }

    #compute_clusters(algorithm, multilevel_resolution, leiden_resolution, walktrap_steps) {
        if (!("graph" in this.#cache)) {
            this.#compute_graph(this.#parameters.scheme);
        }
        this.#cache.clusters = scran.clusterGraph(this.#cache.graph, {
            method: algorithm,
            multiLevelResolution: multilevel_resolution,
            leidenResolution: leiden_resolution,
            leidenModularityObjective: true, // avoid problems with unstable interpretation of leidenResolution.
            walktrapSteps: walktrap_steps
        });
        return;
    }

    /**
     * @return {object} Object containing default parameters,
     * see the `parameters` argument in {@linkcode SnnGraphClusterState#compute compute} for details.
     */
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
     *
     * @param {boolean} run_me - Whether or not to run this step, depending on the clustering method chosen by the user (see {@linkplain ChooseClusteringState}).
     * @param {object} parameters - Parameter object, equivalent to the `snn_graph_cluster` property of the `parameters` of {@linkcode runAnalysis}.
     * @param {number} [parameters.k] - Number of nearest neighbors used to construct the graph.
     * @param {string} [parameters.scheme] - Weighting scheme for graph construction.
     * This can be one of `"rank"`, `"number"` or `"jaccard"`.
     * @param {string} [parameters.algorithm] - Algorithm to use for community detection.
     * This can be one of `"multilevel"`, `"walktrap"` or `"leiden"`.
     * @param {number} [parameters.multilevel_resolution] - Resolution of the multi-level community detection.
     * @param {number} [parameters.leiden_resolution] - Resolution of the Leiden community detection.
     * @param {number} [parameters.walktrap_steps] - Number of merge steps for the Walktrap algorithm.
     *
     * @return The object is updated with the new results.
     */
    compute(run_me, parameters) {
        parameters = utils.defaultizeParameters(parameters, SnnGraphClusterState.defaults());
        this.changed = false;

        if (this.#index.changed || parameters.k !== this.#parameters.k) {
            utils.freeCache(this.#cache.neighbors);
            if (run_me) {
                this.#compute_neighbors(parameters.k);
            } else {
                delete this.#cache.neighbors; // ensuring that this is re-run on future calls to compute() with run_me = true.
            }
            this.changed = true;
        }

        if (this.changed || parameters.scheme !== this.#parameters.scheme) {
            utils.freeCache(this.#cache.graph);
            if (run_me) {
                this.#compute_graph(parameters.scheme);
            } else {
                delete this.#cache.graph;
            }
            this.changed = true 
        }

        if (this.changed 
            || parameters.algorithm !== this.#parameters.algorithm 
            || parameters.multilevel_resolution !== this.#parameters.multilevel_resolution 
            || parameters.leiden_resolution !== this.#parameters.leiden_resolution 
            || parameters.walktrap_steps !== this.#parameters.walktrap_steps 
            || (!this.#valid() && run_me))
        {
            utils.freeCache(this.#cache.clusters);
            if (run_me) {
                this.#compute_clusters(parameters.algorithm, parameters.multilevel_resolution, parameters.leiden_resolution, parameters.walktrap_steps);
            } else {
                delete this.#cache.clusters;
            }
            this.changed = true;
        }

        this.#parameters = parameters;
        return;
    }
}
