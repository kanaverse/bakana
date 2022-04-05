import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as qc from "./quality_control.js";

var cache = {};
var parameters = {};

export var changed = false;

/***************************
 ******** Compute **********
 ***************************/

function raw_compute() {
    var mat = qc.fetchFilteredMatrix();
    var buffer = utils.allocateCachedArray(mat.numberOfColumns(), "Float64Array", cache);

    var discards = qc.fetchDiscards();
    var sums = qc.fetchSums({ unsafe: true }); // Better not have any more allocations in between now and filling of size_factors!

    // Reusing the totals computed earlier.
    var size_factors = buffer.array();
    var j = 0;
    discards.array().forEach((x, i) => {
        if (!x) {
            size_factors[j] = sums[i];
            j++;
        }
    });

    if (j != mat.numberOfColumns()) {
        throw "normalization and filtering are not in sync";
    }

    var block = qc.fetchFilteredBlock();

    utils.freeCache(cache.matrix);
    cache.matrix = scran.logNormCounts(mat, { sizeFactors: buffer, block: block });
    return;
}

export function compute() {
    changed = false;
    if (qc.changed) {
        changed = true;
    } 

    if (changed) {
        raw_compute();
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
    // Token effort.
    let ghandle = handle.createGroup("normalization");
    ghandle.createGroup("parameters"); 
    ghandle.createGroup("results"); 
}

/**************************
 ******** Loading *********
 **************************/

export function unserialize(path) {
    utils.freeCache(cache.matrix);
    delete cache.matrix;
    changed = false;
    return;
}

/***************************
 ******** Getters **********
 ***************************/

export function fetchNormalizedMatrix() {
    if (!("matrix" in cache)) {
        raw_compute();
    }
    return cache.matrix;
}

export function fetchExpression(index) {
    var mat = fetchNormalizedMatrix();
    var buffer = utils.allocateCachedArray(mat.numberOfColumns(), "Float64Array", cache); // re-using the buffer.
    mat.row(index, { buffer: buffer });
    return buffer.slice();
}
