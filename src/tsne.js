import * as scran from "scran.js";
import * as vizutils from "./utils/viz_parent.js";
import * as utils from "./utils/general.js";
import * as neighbor_module from "./neighbor_index.js";
import * as aworkers from "./abstract/worker_parent.js";

/**
 * This creates a t-SNE embedding based on the neighbor index constructed by {@linkplain NeighborIndexState}.
 * This wraps `runTSNE` and related functions from [**scran.js**](https://github.com/jkanche/scran.js).
 * 
 * Methods not documented here are not part of the stable API and should not be used by applications.
 * @hideconstructor
 */
export class TsneState {
    #index;
    #parameters;
    #cache;
    #reloaded;

    #worker;
    #worker_id;

    #ready;
    #run;

    constructor(index, parameters = null, reloaded = null) {
        if (!(index instanceof neighbor_module.NeighborIndexState)) {
            throw new Error("'index' should be a State object from './neighbor_index.js'");
        }
        this.#index = index;

        this.#parameters = (parameters === null ? {} : parameters);
        this.#cache = { "counter": 0, "promises": {} };
        this.#reloaded = reloaded;
        this.changed = false;

        let worker = aworkers.createTsneWorker();
        let { worker_id, ready } = vizutils.initializeWorker(worker, this.#cache, vizutils.scranOptions);
        this.#worker = worker;
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
        this.#run = vizutils.runWithNeighbors(this.#worker, args, nn_out, this.#cache);
        return;
    }

    /**
     * This method should not be called directly by users, but is instead invoked by {@linkcode runAnalysis}.
     * Each argument is taken from the property of the same name in the `tsne` property of the `parameters` of {@linkcode runAnalysis}.
     *
     * @param {number} perplexity - Number specifying the perplexity for the probability calculations.
     * @param {number} iterations - Number of iterations to run the algorithm.
     * @param {boolean} animate - Whether o process animation iterations, see {@linkcode setVisualizationAnimate} for details.
     *
     * @return t-SNE coordinates are computed in parallel on a separate worker thread.
     * A promise is returned that resolves when those calculations are complete.
     */
    compute(perplexity, iterations, animate) {
        let same_neighbors = (!this.#index.changed && perplexity === this.#parameters.perplexity);
        if (same_neighbors && iterations == this.#parameters.iterations) {
            this.changed = false;
            return;
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
     ******** Results **********
     ***************************/

    async #fetch_results(copy) {
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
            return vizutils.sendTask(this.#worker, { "cmd": "FETCH" }, this.#cache);
        }
    }

    /**
     * Obtain a summary of the state, typically for display on a UI like **kana**.
     *
     * @return A promise that resolves to an object containing:
     *
     * - `x`: a Float64Array containing the x-coordinate for each cell.
     * - `y`: a Float64Array containing the y-coordinate for each cell.
     * - `iterations`: the number of iterations processed.
     */
    summary() {
        return this.#fetch_results(true);
    }

    /*************************
     ******** Saving *********
     *************************/

    async serialize(handle) {
        let ghandle = handle.createGroup("tsne");

        {
            let phandle = ghandle.createGroup("parameters");
            phandle.writeDataSet("perplexity", "Float64", [], this.#parameters.perplexity);
            phandle.writeDataSet("iterations", "Int32", [], this.#parameters.iterations);
            phandle.writeDataSet("animate", "Uint8", [], Number(this.#parameters.animate));
        }

        {
            let res = await this.#fetch_results(false);
            let rhandle = ghandle.createGroup("results");
            rhandle.writeDataSet("x", "Float64", null, res.x);
            rhandle.writeDataSet("y", "Float64", null, res.y);
        }

        return;
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
            return vizutils.sendTask(this.#worker, { "cmd": "RERUN" }, this.#cache);
        }
    }
}

/**************************
 ******** Loading *********
 **************************/

export function unserialize(handle, index) {
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

    return {
        state: new TsneState(index, parameters, reloaded),
        parameters: { ...parameters }
    };
}
