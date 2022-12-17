import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as correct_module from "./batch_correction.js";

/**
 * This step performs k-means clustering on the PCs, wrapping the `clusterKmeans` function from [**scran.js**](https://github.com/jkanche/scran.js).
 *
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class KmeansClusterState {
    #correct;
    #parameters;
    #cache;

    constructor(correct, parameters = null, cache = null) {
        if (!(correct instanceof correct_module.BatchCorrectionState)) {
            throw new Error("'correct' should be a BatchCorrectionState object");
        }
        this.#correct = correct;

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

    /**
     * @return {Int32WasmArray} Array of cluster assignments for each cell in the (filtered) dataset,
     * available after running {@linkcode KmeansClusterState#compute compute}.
     */
    fetchClusters() {
        if (!this.#valid()) {
            throw new Error("cannot fetch k-means clusters from an invalid state");
        } else {
            return this.#cache.raw.clusters({ copy: "view" });
        }
    }

    /**
     * @return {object} Object containing the parameters.
     */
    fetchParameters() {
        return { ...this.#parameters };
    };

    /***************************
     ******** Compute **********
     ***************************/

    #valid() {
        return "raw" in this.#cache;
    }

    /** 
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     * Each argument is taken from the property of the same name in the `kmeans_cluster` property of the `parameters` of {@linkcode runAnalysis}.
     * The exception is `run_me`, which is computed internally and does not correspond to any parameter in `kmeans_cluster`.
     *
     * @param {boolean} run_me - Whether or not to run this step, depending on the clustering method chosen by the user (see {@linkplain ChooseClusteringState}).
     * @param {number} k - Number of clusters to create.
     *
     * @return The object is updated with the new results.
     */
    compute(run_me, k) {
        this.changed = false;

        if (this.#correct.changed || k != this.#parameters.k || (!this.#valid() && run_me)) {
            utils.freeCache(this.#cache.raw);

            if (run_me) {
                var pcs = this.#correct.fetchCorrected();
                this.#cache.raw = scran.clusterKmeans(pcs, k, { 
                    numberOfDims: this.#correct.fetchNumberOfDimensions(),
                    numberOfCells: this.#correct.fetchNumberOfCells(),
                    initMethod: "pca-part" 
                });
            } else {
                delete this.#cache.raw; // ensure this step gets re-run later when run_me = true. 
            }

            this.#parameters.k = k;
            this.changed = true;
        }

        return;
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
                let clusters = this.fetchClusters();
                rhandle.writeDataSet("clusters", "Int32", null, clusters);
             }
        }

        return;
    }
}

/**************************
 ******** Loading *********
 **************************/

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
                cache.raw = scran.emptyClusterKmeansResults(clusters.length, parameters.k, pca.fetchNumberOfDimensions());
                cache.raw.clusters({ copy: false }).set(clusters);
            }
        }
    }

    return new KmeansClusterState(pca, parameters, cache);
}
