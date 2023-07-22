import * as scran from "scran.js";
import * as vizutils from "./utils/viz_parent.js";
import * as utils from "./utils/general.js";
import * as neighbor_module from "./neighbor_index.js";
import * as aworkers from "./abstract/worker_parent.js";

/**
 * This creates a UMAP embedding based on the neighbor index constructed at {@linkplain NeighborIndexState}.
 * This wraps [`runUMAP`](https://kanaverse.github.io/scran.js/global.html#runUMAP)
 * and related functions from [**scran.js**](https://github.com/kanaverse/scran.js).
 * 
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class UmapState {
    #index;
    #parameters;
    #reloaded;

    #worker_id;

    #ready;
    #run;

    constructor(index, parameters = null, reloaded = null) {
        if (!(index instanceof neighbor_module.NeighborIndexState)) {
            throw new Error("'index' should be a State object from './neighbor_index.js'");
        }
        this.#index = index;

        this.#parameters = (parameters === null ? {} : parameters);
        this.#reloaded = reloaded;
        this.changed = false;

        let worker = aworkers.createUmapWorker();
        let { worker_id, ready } = vizutils.initializeWorker(worker, vizutils.scranOptions);
        this.#worker_id = worker_id;
        this.#ready = ready;

        this.#run = null;
    }

    ready() {
        // It is assumed that the caller will await the ready()
        // status before calling any other methods of this instance.
        return this.#ready;
    }

    free() {
        return vizutils.killWorker(this.#worker_id);
    }

    /***************************
     ******** Getters **********
     ***************************/

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.copy=true] - Whether to create a copy of the coordinates,
     * if the caller might mutate them.
     *
     * @return {object} Object containing:
     *
     * - `x`: a Float64Array containing the x-coordinate for each cell.
     * - `y`: a Float64Array containing the y-coordinate for each cell.
     * - `iterations`: the number of iterations processed.
     *
     * @async
     */
    async fetchResults({ copy = true } = {}) {
        if (this.#reloaded !== null) {
            let output = {
                x: this.#reloaded.x,
                y: this.#reloaded.y
            };

            if (copy) {
                output.x = output.x.slice();
                output.y = output.y.slice();
            }

            output.iterations = this.#parameters.num_epochs;
            return output;
        } else {
            // Vectors that we get from the worker are inherently
            // copied, so no need to do anything extra here.
            await this.#run;
            return vizutils.sendTask(this.#worker_id, { "cmd": "FETCH" });
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

    #core(num_neighbors, num_epochs, min_dist, animate, reneighbor) {
        var nn_out = null;
        if (reneighbor) {
            nn_out = vizutils.computeNeighbors(this.#index, num_neighbors);
        }

        let args = {
            "num_neighbors": num_neighbors,
            "num_epochs": num_epochs,
            "min_dist": min_dist,
            "animate": animate
        };

        // This returns a promise but the message itself is sent synchronously,
        // which is important to ensure that the UMAP runs in its worker in
        // parallel with other analysis steps. Do NOT put the runWithNeighbors
        // call in a .then() as this may defer the message sending until 
        // the current thread is completely done processing.
        this.#run = vizutils.runWithNeighbors(this.#worker_id, args, nn_out);
        return;
    }

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     *
     * @param {object} parameters - Parameter object, equivalent to the `umap` property of the `parameters` of {@linkcode runAnalysis}.
     * @param {number} parameters.num_neighbors - Number of neighbors to use to construct the simplicial sets.
     * @param {number} parameters.num_epochs - Number of epochs to run the algorithm.
     * @param {number} parameters.min_dist - Number specifying the minimum distance between points.
     * @param {boolean} parameters.animate - Whether to process animation iterations, see {@linkcode setVisualizationAnimate} for details.
     *
     * @return UMAP coordinates are computed in parallel on a separate worker thread.
     * A promise that resolves when the calculations are complete.
     */
    compute(parameters) {
        let { num_neighbors, num_epochs, min_dist, animate } = parameters;

        let same_neighbors = (!this.#index.changed && this.#parameters.num_neighbors === num_neighbors);
        if (same_neighbors && num_epochs === this.#parameters.num_epochs && min_dist === this.#parameters.min_dist) {
            this.changed = false;
            return new Promise(resolve => resolve(null));
        }

        // In the reloaded state, we must send the neighbor
        // information, because it hasn't ever been sent before.
        if (this.#reloaded !== null) {
            same_neighbors = false;
            this.#reloaded = null;
        }

        this.#core(num_neighbors, num_epochs, min_dist, animate, !same_neighbors);

        this.#parameters.num_neighbors = num_neighbors;
        this.#parameters.num_epochs = num_epochs;
        this.#parameters.min_dist = min_dist;
        this.#parameters.animate = animate;

        this.changed = true;
        return this.#run;
    }

    /***************************
     ******** Getters **********
     ***************************/

    /**
     * Repeat the animation iterations.
     * It is assumed that {@linkcode setVisualizationAnimate} has been set appropriately to process each iteration.
     *
     * @return A promise that resolves on successful completion of all iterations.
     */
    animate() {
        if (this.#reloaded !== null) {
            this.#reloaded = null;

            // We need to reneighbor because we haven't sent the neighbors across yet.
            this.#core(this.#parameters.num_neighbors, this.#parameters.num_epochs, this.#parameters.min_dist, true, true);
      
            // Mimicking the response from the re-run.
            return this.#run
                .then(contents => { 
                    return {
                        "type": "umap_rerun",
                        "data": { "status": "SUCCESS" }
                    };
                });
        } else {
            return vizutils.sendTask(this.#worker_id, { "cmd": "RERUN" });
        }
    }
}
