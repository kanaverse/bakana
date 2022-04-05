import * as scran from "scran.js";
import * as utils from "./utils/general.js";
import * as pca from "./pca.js";

var cache = {};
var parameters = {};

export var changed = false;

/***************************
 ******** Compute **********
 ***************************/

export function raw_compute(approximate) {
    var pcs = pca.fetchPCs();
    cache.raw = scran.buildNeighborSearchIndex(pcs.pcs, { approximate: approximate, numberOfDims: pcs.num_pcs, numberOfCells: pcs.num_obs });
    return;
}

export function compute(approximate) {
    changed = false;

    if (pca.changed || approximate != parameters.approximate) {
        utils.freeCache(cache.raw);
        raw_compute(approximate);
        parameters.approximate = approximate;
        changed = true;
    }

    return;
}

/***************************
 ******** Results **********
 ***************************/

export function results() {
    return {};
}

/*************************
 ******** Saving *********
 *************************/

export function serialize(handle) {
    let ghandle = handle.createGroup("neighbor_index");

    {
        let phandle = ghandle.createGroup("parameters");
        phandle.writeDataSet("approximate", "Uint8", [], Number(parameters.approximate));
    }

    ghandle.createGroup("results");
    return;
}

/**************************
 ******** Loading *********
 **************************/

export function unserialize(handle) {
    let ghandle = handle.open("neighbor_index");

    {
        let phandle = ghandle.open("parameters");
        parameters = {
            approximate: phandle.open("approximate", { load: true }).value > 0
        };
    }

    utils.freeCache(cache.raw);
    cache = {};

    return { ...parameters };
}

/***************************
 ******** Getters **********
 ***************************/

export function fetchIndex() {
    if (!("raw" in cache)) {
        raw_compute(parameters.approximate);
    }
    return cache.raw;
}
