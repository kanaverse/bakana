import * as scran from "scran.js";
import * as utils from "./utils/general.js";
import * as snn_module from "./snn_graph_cluster.js";
import * as kmeans_module from "./kmeans_cluster.js";

export class State {
    #snn_cluster;
    #kmeans_cluster;
    #parameters;
    #cache;

    constructor(snn, kmeans, parameters = null, cache = null) {
        if (!(snn instanceof snn_module.State)) {
            throw new Error("'snn' should be a State object from './snn_graph_cluster.js'");
        }
        this.#snn_cluster = snn;

        if (!(kmeans instanceof kmeans_module.State)) {
            throw new Error("'kmeans' should be a State object from './kmeans_cluster.js'");
        }
        this.#kmeans_cluster = kmeans;

        this.#parameters = (parameters === null ? {} : parameters);
        this.#cache = (cache === null ? {} : cache);
        this.changed = false;
    }

    free() {}

    /***************************
     ******** Getters **********
     ***************************/

    fetchClustersAsWasmArray() {
        if (this.#parameters.method == "snn_graph") {
            return this.#snn_cluster.fetchClustersAsWasmArray();
        } else if (this.#parameters.method == "kmeans") {
            return this.#kmeans_cluster.fetchClustersAsWasmArray();
        }
    }

    /***************************
     ******** Compute **********
     ***************************/

    compute(method) {
        this.changed = true;
        
        if (method == this.#parameters.method) {
            if (method == "snn_graph") {
                if (!this.#snn_cluster.changed) {
                    this.changed = false;
                }
            } else if (method == "kmeans") {
                if (!this.#kmeans_cluster.changed) {
                    this.changed = false;
                }
            }
        }

        this.#parameters.method = method;
        return;
    }

    /***************************
     ******** Results **********
     ***************************/

    results() {
        var clusters = this.fetchClustersAsWasmArray();
        return { "clusters": clusters.slice() };
    }

    /*************************
     ******** Saving *********
     *************************/

    serialize(handle) {
        let ghandle = handle.createGroup("choose_clustering");

        {
            let phandle = ghandle.createGroup("parameters");
            phandle.writeDataSet("method", "String", [], this.#parameters.method);
        }

        // No need to serialize the cluster IDs as this is done for each step.
        ghandle.createGroup("results");
        return;
    }
}

/**************************
 ******** Loading *********
 **************************/

export function unserialize(handle, snn, kmeans) {
    let ghandle = handle.open("choose_clustering");

    let parameters;
    {
        let phandle = ghandle.open("parameters");
        parameters = {
            method: phandle.open("method", { load: true }).values[0]
        };
    }

    let cache = {};

    return {
        state: new State(snn, kmeans, parameters, cache),
        parameters: { ...parameters }
    };
}
