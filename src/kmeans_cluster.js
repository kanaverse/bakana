import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as pca_module from "./pca.js";

export class State {
    #pca;
    #parameters;
    #cache;

    constructor(pca, parameters = null, cache = null) {
        if (!(pca instanceof pca_module.State)) {
            throw new Error("'pca' should be a State object from './pca.js'");
        }
        this.#pca = pca;

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

    fetchClustersAsWasmArray() {
        if (!this.#valid()) {
            throw new Error("cannot fetch k-means clusters from an invalid state");
        } else {
            return this.#cache.raw.clusters({ copy: "view" });
        }
    }

    /***************************
     ******** Compute **********
     ***************************/

    #valid() {
        return "raw" in this.#cache;
    }

    compute(run_me, k) {
        this.changed = false;

        if (this.#pca.changed || k != this.#parameters.k || (!this.#valid() && run_me)) {
            utils.freeCache(this.#cache.raw);

            if (run_me) {
                var pcs = this.#pca.fetchPCs();
                this.#cache.raw = scran.clusterKmeans(pcs.pcs, k, { numberOfDims: pcs.num_pcs, numberOfCells: pcs.num_obs, initMethod: "pca-part" });
            } else {
                delete this.#cache.raw; // ensure this step gets re-run later when run_me = true. 
            }

            this.#parameters.k = k;
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
        let ghandle = handle.createGroup("kmeans_cluster");

        {
            let phandle = ghandle.createGroup("parameters");
            phandle.writeDataSet("k", "Int32", [], this.#parameters.k);
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

class KmeansMimic {
    constructor(clusters) {
        this.buffer = scran.createInt32WasmArray(clusters.length);
        this.buffer.set(clusters);
    }

    clusters({ copy }) {
        return utils.mimicGetter(this.buffer, copy);
    }

    free() {
        this.buffer.free();
    }
}

export function unserialize(handle, pca) {
    let parameters = {
        k: 10
    };
    let cache = {};

    // Protect against old analysis states that don't have kmeans_cluster.
    if ("kmeans_cluster" in handle.children) {
        let ghandle = handle.open("kmeans_cluster");

        {
            let phandle = ghandle.open("parameters");
            parameters.k = phandle.open("k", { load: true }).values[0];
        }

        {
            let rhandle = ghandle.open("results");
            if ("clusters" in rhandle.children) {
                let clusters = rhandle.open("clusters", { load: true }).values;
                cache.raw = new KmeansMimic(clusters);
            }
        }
    }

    return {
        state: new State(pca, parameters, cache),
        parameters: { ...parameters }
    };
}
