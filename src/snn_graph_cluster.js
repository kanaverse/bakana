import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as index from "./neighbor_index.js";

var cache = {};
var parameters = {};

export var changed = false;

/***************************
 ******** Compute **********
 ***************************/

function valid() {
    return "clusters" in cache;
}

export function compute_neighbors(k) {
    cache.neighbors = scran.findNearestNeighbors(index.fetchIndex(), k);
    return;
}

export function compute_graph(scheme) {
    if (!("neighbors" in cache)) { // need to check as reloaded state will not populate the internals.
        compute_neighbors(parameters.k);
    }
    cache.graph = scran.buildSNNGraph(cache.neighbors, { scheme: scheme });
    return;
}

export function compute_clusters(resolution) {
    if (!("graph" in cache)) {
        compute_graph(parameters.scheme);
    }
    cache.clusters = scran.clusterSNNGraph(cache.graph, { resolution: resolution });
    return;
}

export function compute(run_me, k, scheme, resolution) {
    changed = false;

    if (index.changed || k !== parameters.k) {
        utils.freeCache(cache.neighbors);
        if (run_me) {
            compute_neighbors(k);
        } else {
            delete cache.neighbors; // ensuring that this is re-run on future calls to compute() with run_me = true.
        }
        parameters.k = k;
        changed = true;
    }

    if (changed || scheme !== parameters.scheme) {
        utils.freeCache(cache.graph);
        if (run_me) {
            compute_graph(scheme);
        } else {
            delete cache.graph;
        }
        parameters.scheme = scheme;
        changed = true 
    }

    if (changed || resolution !== parameters.resolution || (!valid() && run_me)) {
        utils.freeCache(cache.clusters);
        if (run_me) {
            compute_clusters(resolution);
        } else {
            delete cache.clusters;
        }
        parameters.resolution = resolution;
        changed = true;
    }

    return;
}

/***************************
 ******** Results **********
 ***************************/

export function results() {
    // Cluster IDs will be passed to main thread in 
    // choose_clustering, so no need to do it here.
    return {};
}

/*************************
 ******** Saving *********
 *************************/

export function serialize(handle) {
    let ghandle = handle.createGroup("snn_graph_cluster");

    {
        let phandle = ghandle.createGroup("parameters");
        phandle.writeDataSet("k", "Int32", [], parameters.k);
        phandle.writeDataSet("scheme", "String", [], ["rank", "number", "jaccard"][parameters.scheme]); // TODO: parameters.scheme should just directly be the string.
        phandle.writeDataSet("resolution", "Float64", [], parameters.resolution);
    }

    {
        let rhandle = ghandle.createGroup("results");
        if (valid()) {
            let clusters = fetchClustersAsWasmArray();
            rhandle.writeDataSet("clusters", "Int32", null, clusters);
        }
    }

    return;
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

export function unserialize(handle) {
    let ghandle = handle.open("snn_graph_cluster");

    {
        let phandle = ghandle.open("parameters");
        parameters = {
            k: phandle.open("k", { load: true }).values[0],
            scheme: phandle.open("scheme", { load: true }).values[0],
            resolution: phandle.open("resolution", { load: true }).values[0]
        };
        parameters.scheme = { "rank": 0, "number": 1, "jaccard": 2 }[parameters.scheme];
    }

    utils.freeCache(cache.neighbors);
    utils.freeCache(cache.graph);
    utils.freeCache(cache.clusters);
    cache = {};
    changed = false;

    {
        let rhandle = ghandle.open("results");
        if ("clusters" in rhandle.children) {
            let clusters = rhandle.open("clusters", { load: true }).values;
            cache.clusters = new SNNClusterMimic(clusters);
        }
    }

    return { ...parameters };
}

/***************************
 ******** Getters **********
 ***************************/

export function fetchClustersAsWasmArray() {
    if (!valid()) {
        throw "cannot fetch SNN clusters from an invalid state";
    } else {
        return cache.clusters.membership({ copy: "view" });
    }
}


