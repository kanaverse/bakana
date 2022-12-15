import * as scran from "scran.js";
import * as utils from "./utils/general.js";
import * as vizutils from "./utils/viz_child.js";
import * as aworkers from "./abstract/worker_child.js";

var cache = {};
var init_changed = false;
var init_parameters = {};
var run_parameters = {};

function rerun(animate, iterations) {
    var num_obs = cache.init.numberOfCells(); 
    var delay = vizutils.chooseDelay(animate);
    var current_status = cache.init.clone();

    try {
        for (; current_status.iterations() < iterations; ) {
            scran.runTSNE(current_status, { runTime: delay, maxIterations: iterations }); 
  
            if (animate) {
                let xy = current_status.extractCoordinates();
                aworkers.sendMessage({
                    "type": "tsne_iter",
                    "x": xy.x,
                    "y": xy.y,
                    "iteration": current_status.iterations()
                }, [xy.x.buffer, xy.y.buffer]);
            }
        }
        cache.final = current_status.extractCoordinates();

    } finally {
        current_status.free();
    }
}

var loaded;
aworkers.registerCallback(msg => {
    var id = msg.data.id;
  
    if (msg.data.cmd == "INIT") {
        loaded = scran.initialize(msg.data.scranOptions); 
        loaded
            .then(x => {
                aworkers.sendMessage({
                    "id": id,
                    "type": "init_worker",
                    "data": { "status": "SUCCESS" }
                });
            })
            .catch(error => {
                aworkers.sendMessage({ 
                    "id": id,
                    "type": "error",
                    "error": error
                });
            });
  
    } else if (msg.data.cmd == "RUN") {
        loaded
            .then(x => {
                var new_neighbors;
                if ("neighbors" in msg.data) {
                    utils.freeCache(cache.neighbors);
                    cache.neighbors = vizutils.recreateNeighbors(msg.data.neighbors);
                    new_neighbors = true;
                } else {
                    new_neighbors = false;
                }

                var init_args = { "perplexity": msg.data.params.perplexity };
                if (!new_neighbors && !utils.changedParameters(init_args, init_parameters)) {
                    init_changed = false;
                } else {
                    utils.freeCache(cache.init);
                    cache.init = scran.initializeTSNE(cache.neighbors, { perplexity: init_args.perplexity });
                    init_parameters = init_args;
                    init_changed = true;
                }

                // Nothing downstream depends on the run results, so we don't set any changed flag.
                var run_args = { "iterations": msg.data.params.iterations };
                if (init_changed || utils.changedParameters(run_args, run_parameters)) {
                    rerun(msg.data.params.animate, run_args.iterations);
                    run_parameters = run_args;
                }
          
                aworkers.sendMessage({
                    "id": id,
                    "type": "tsne_run",
                    "data": { "status": "SUCCESS" }
                });
            })
            .catch(error => {
                aworkers.sendMessage({ 
                    "id": id,
                    "type": "error",
                    "error": error
                });
            });
  
    } else if (msg.data.cmd == "RERUN") {
        loaded
            .then(x => {
                rerun(true, run_parameters.iterations);
                aworkers.sendMessage({
                    "id": id,
                    "type": "tsne_rerun",
                    "data": { "status": "SUCCESS" }
                });
            })
            .catch(error => {
                aworkers.sendMessage({ 
                    "id": id,
                    "type": "error",
                    "error": error
                });
            });

    } else if (msg.data.cmd == "FETCH") {
        loaded
            .then(x => {
                var info = {
                    "x": cache.final.x.slice(),
                    "y": cache.final.y.slice(),
                    "iterations": run_parameters.iterations
                };
                var transfer = [info.x.buffer, info.y.buffer];
                aworkers.sendMessage({
                    "id": id,
                    "type": "tsne_fetch",
                    "data": info
                }, transfer);
            })
            .catch(error => {
                aworkers.sendMessage({ 
                    "id": id,
                    "type": "error",
                    "error": error
                });
            });

    } else if (msg.data.cmd == "KILL") {
        loaded
            .then(x => {
                scran.terminate();
                aworkers.sendMessage({
                    "id": id,
                    "type": "tsne_killed",
                    "data": null
                });
            })
            .catch(error => {
                aworkers.sendMessage({ 
                    "id": id,
                    "type": "error",
                    "error": error
                });
            });

    } else {
        aworkers.sendMessage({
            "id": id,
            "type": "error",
            "error": "unknown message type '" + JSON.stringify(msg.data) + "'"
        });
    }
});
