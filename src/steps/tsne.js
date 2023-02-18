import * as scran from "scran.js";
import * as vizutils from "./utils/viz_parent.js";
import * as utils from "./utils/general.js";
import * as neighbor_module from "./neighbor_index.js";
import * as aworkers from "./abstract/worker_parent.js";

/**
 * This creates a t-SNE embedding based on the neighbor index constructed by {@linkplain NeighborIndexState}.
 * This wraps [`runTSNE`](https://kanaverse.github.io/scran.js/global.html#runTSNE)
 * and related functions from [**scran.js**](https://github.com/kanaverse/scran.js).
 * 
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class TsneState {
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

        let worker = aworkers.createTsneWorker();
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
     * @return {object} Object containing the parameters.
     */
    fetchParameters() {
        return { ...this.#parameters }; // avoid pass-by-reference links.
    }

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
        
            output.iterations = this.#parameters.iterations;
            return output;
        } else {
            // Vectors that we get from the worker are inherently
            // copied, so no need to do anything extra here.
            await this.#run;
            return vizutils.sendTask(this.#worker_id, { "cmd": "FETCH" });
        }
    }


    /***************************
     ******** Compute **********
     ***************************/

    #core(perplexity, iterations, animate, reneighbor) {
        var nn_out = null;
        if (reneighbor) {
            var k = scran.perplexityToNeighbors(perplexity);
            nn_out = vizutils.computeNeighbors(this.#index, k);
        }

        let args = {
            "perplexity": perplexity,
            "iterations": iterations,
            "animate": animate
        };

        // This returns a promise but the message itself is sent synchronously,
        // which is important to ensure that the t-SNE runs in its worker in
        // parallel with other analysis steps. Do NOT put the runWithNeighbors
        // call in a .then() as this may defer the message sending until 
        // the current thread is completely done processing.
        this.#run = vizutils.runWithNeighbors(this.#worker_id, args, nn_out);
        return;
    }

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     *
     * @param {object} parameters - Parameter object, equivalent to the `tsne` property of the `parameters` of {@linkcode runAnalysis}.
     * @param {number} parameters.perplexity - Number specifying the perplexity for the probability calculations.
     * @param {number} parameters.iterations - Number of iterations to run the algorithm.
     * @param {boolean} parameters.animate - Whether to process animation iterations, see {@linkcode setVisualizationAnimate} for details.
     *
     * @return t-SNE coordinates are computed in parallel on a separate worker thread.
     * A promise is returned that resolves when those calculations are complete.
     */
    compute(parameters) {
        let { perplexity, iterations, animate } = parameters;

        let same_neighbors = (!this.#index.changed && perplexity === this.#parameters.perplexity);
        if (same_neighbors && iterations == this.#parameters.iterations) {
            this.changed = false;
            return new Promise(resolve => resolve(null));
        }

        // In the reloaded state, we must send the neighbor
        // information, because it hasn't ever been sent before.
        if (this.#reloaded !== null) {
            same_neighbors = false;
            this.#reloaded = null;
        }

        this.#core(perplexity, iterations, animate, !same_neighbors);

        this.#parameters.perplexity = perplexity;
        this.#parameters.iterations = iterations;
        this.#parameters.animate = animate;

        this.changed = true;
        return this.#run;
    }

    /***************************
     ******* Animators *********
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
            this.#core(this.#parameters.perplexity, this.#parameters.iterations, true, true);

            // Mimicking the response from the re-run.
            return this.#run
                .then(contents => {
                    return {
                        "type": "tsne_rerun",
                        "data": { "status": "SUCCESS" }
                    };
                });
        } else {
            return vizutils.sendTask(this.#worker_id, { "cmd": "RERUN" });
        }
    }
}

/**************************
 ******** Loading *********
 **************************/

export async function unserialize(handle, index) {
    let ghandle = handle.open("tsne");

    let parameters;
    {
        let phandle = ghandle.open("parameters");
        parameters = {
            perplexity: phandle.open("perplexity", { load: true }).values[0],
            iterations: phandle.open("iterations", { load: true }).values[0],
            animate: phandle.open("animate", { load: true }).values[0] > 0
        };
    }

    let reloaded;
    {
        let rhandle = ghandle.open("results");
        reloaded = {
            x: rhandle.open("x", { load: true }).values,
            y: rhandle.open("y", { load: true }).values
        };
    }

    let output = new TsneState(index, parameters, reloaded);
    await output.ready();
    return output;
}
