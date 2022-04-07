import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as neighbor_module from "./neighbor_index.js";

/**
 * This step does SNN graph clustering based on the neighbor search index built by {@linkcode neighbor_index}.
 * This wraps `clusterSNNGraph` and related functions from [**scran.js**](https://github.com/jkanche/scran.js).
 *
 * The parameters in {@linkcode runAnalysis} should be an object containing:
 *
 * - `k`: the number of nearest neighbors used to construct the graph.
 * - `scheme`: string specifying the weighting scheme for graph construction.
 *   This can be one of `"rank"`, `"number"` or `"jaccard"`.
 * - `resolution`: number containing the resolution of the community detection.
 *
 * Calling the **`results()`** method for the relevant state instance will return an empty object.
 * See the {@linkcode choose_clustering} step to actually obtain the cluster assignments.
 * 
 * Methods not documented here are not part of the stable API and should not be used by applications.
 *
 * @namespace snn_graph_cluster
 */

export class State {
    #index;
    #parameters;
    #cache;

    constructor(index, parameters = null, cache = null) {
        if (!(index instanceof neighbor_module.State)) {
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

    fetchClustersAsWasmArray() {
        if (!this.#valid()) {
            throw "cannot fetch SNN clusters from an invalid state";
        } else {
            return this.#cache.clusters.membership({ copy: "view" });
        }
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

    #compute_clusters(resolution) {
        if (!("graph" in this.#cache)) {
            this.#compute_graph(this.#parameters.scheme);
        }
        this.#cache.clusters = scran.clusterSNNGraph(this.#cache.graph, { resolution: resolution });
        return;
    }

    compute(run_me, k, scheme, resolution) {
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

        if (this.changed || resolution !== this.#parameters.resolution || (!this.#valid() && run_me)) {
            utils.freeCache(this.#cache.clusters);
            if (run_me) {
                this.#compute_clusters(resolution);
            } else {
                delete this.#cache.clusters;
            }
            this.#parameters.resolution = resolution;
            this.changed = true;
        }

        return;
    }

    /***************************
     ******** Results **********
     ***************************/

    results() {
        // Cluster IDs will be passed to main thread in 
        // choose_clustering, so no need to do it here.
        return {};
    }

    /*************************
     ******** Saving *********
     *************************/

    serialize(handle) {
        let ghandle = handle.createGroup("snn_graph_cluster");

        {
            let phandle = ghandle.createGroup("parameters");
            phandle.writeDataSet("k", "Int32", [], this.#parameters.k);
            phandle.writeDataSet("scheme", "String", [], this.#parameters.scheme);
            phandle.writeDataSet("resolution", "Float64", [], this.#parameters.resolution);
        }

        {
            let rhandle = ghandle.createGroup("results");
            if (this.#valid()) {
                let clusters = this.fetchClustersAsWasmArray();
                rhandle.writeDataSet("clusters", "Int32", null, clusters);
            }
        }

        return;
    }
}

/**************************
 ******** Loading *********
 **************************/

class SNNClusterMimic {
    constructor(clusters) {
        this.buffer = scran.createInt32WasmArray(clusters.length);
        this.buffer.set(clusters);
    }

    membership({ copy }) {
        return utils.mimicGetter(this.buffer, copy);
    }

    free() {
        this.buffer.free();
    }
}

export function unserialize(handle, index) {
    let ghandle = handle.open("snn_graph_cluster");

    let parameters = {};
    {
        let phandle = ghandle.open("parameters");
        parameters.k = phandle.open("k", { load: true }).values[0];
        parameters.scheme = phandle.open("scheme", { load: true }).values[0];
        parameters.resolution = phandle.open("resolution", { load: true }).values[0];
    }

    let cache = {};
    {
        let rhandle = ghandle.open("results");
        if ("clusters" in rhandle.children) {
            let clusters = rhandle.open("clusters", { load: true }).values;
            cache.clusters = new SNNClusterMimic(clusters);
        }
    }

    return {
        state: new State(index, parameters, cache),
        parameters: { ...parameters }
    };
}


