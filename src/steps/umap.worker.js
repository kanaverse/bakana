import * as scran from "scran.js";
import * as utils from "./utils/general.js";
import * as vizutils from "./utils/viz_child.js";
import * as aworkers from "../abstract/worker_child.js";

var cache = {};
var init_changed = false;
var init_parameters = {};
var run_parameters = {};

function rerun(animate) {
    var delay = vizutils.chooseDelay(animate);
    var current_status = cache.init.clone();

    try {
        cache.total = current_status.totalEpochs();
        for (; current_status.currentEpoch() < cache.total; ) {
            scran.runUMAP(current_status, { runTime: delay });

            if (animate) {
                var xy = current_status.extractCoordinates();
                aworkers.sendMessage({
                    "type": "umap_iter",
                    "x": xy.x,
                    "y": xy.y,
                    "iteration": current_status.currentEpoch()
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
        
                var init_args = { "min_dist": msg.data.params.min_dist, "num_epochs": msg.data.params.num_epochs };
                if (!new_neighbors && !utils.changedParameters(init_args, init_parameters)) {
                    init_changed = false;
                } else {
                    utils.freeCache(cache.init);
                    cache.init = scran.initializeUMAP(cache.neighbors, { epochs: init_args.num_epochs, minDist: init_args.min_dist });
                    init_parameters = init_args;
                    init_changed = true;
                }
        
                // Nothing downstream depends on the run results, so we don't set any changed flag.
                var run_args = {};
                if (init_changed || utils.changedParameters(run_args, run_parameters)) {
                    rerun(msg.data.params.animate);
                    run_parameters = run_args;
                }
        
                aworkers.sendMessage({
                    "id": id,
                    "type": "umap_run",
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
                rerun(true);
                aworkers.sendMessage({
                    "id": id,
                    "type": "umap_rerun",
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
                    "iterations": cache.total
                };
                var transfer = [info.x.buffer, info.y.buffer];
                aworkers.sendMessage({
                    "id": id,
                    "type": "umap_fetch",
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
    }
});
