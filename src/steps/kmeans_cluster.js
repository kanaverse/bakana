import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as correct_module from "./batch_correction.js";

export const step_name = "kmeans_cluster";

/**
 * This step performs k-means clustering on the PCs, 
 * wrapping the [`clusterKmeans`](https://kanaverse.github.io/scran.js/global.html#clusterKmeans) function 
 * from [**scran.js**](https://github.com/kanaverse/scran.js).
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
     * @return {object} Object containing default parameters,
     * see the `parameters` argument in {@linkcode KmeansClusterState#compute compute} for details.
     */
    static defaults() {
        return { k: 10 };
    }

    /** 
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     *
     * @param {boolean} run_me - Whether or not to run this step, depending on the clustering method chosen by the user (see {@linkplain ChooseClusteringState}).
     * @param {object} parameters - Parameter object, equivalent to the `choose_clustering` property of the `parameters` of {@linkcode runAnalysis}.
     * @param {number} [parameters.k] - Number of clusters to create.
     *
     * @return The object is updated with the new results.
     */
    compute(run_me, parameters) {
        parameters = utils.defaultizeParameters(parameters, KmeansClusterState.defaults());
        this.changed = false;

        if (this.#correct.changed || parameters.k != this.#parameters.k || (!this.#valid() && run_me)) {
            utils.freeCache(this.#cache.raw);

            if (run_me) {
                var pcs = this.#correct.fetchCorrected();
                this.#cache.raw = scran.clusterKmeans(
                    pcs,
                    parameters.k,
                    { 
                        numberOfDims: this.#correct.fetchNumberOfDimensions(),
                        numberOfCells: this.#correct.fetchNumberOfCells(),
                        initMethod: "var-part" 
                    }
                );
            } else {
                delete this.#cache.raw; // ensure this step gets re-run later when run_me = true. 
            }

            this.changed = true;
        }

        this.#parameters = parameters;
        return;
    }
}
