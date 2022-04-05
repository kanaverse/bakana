import * as scran from "scran.js"; 
import * as utils from "./utils/general.js";
import * as inputs from "./inputs.js";
import { mito } from "./mito.js";

var cache = {};
var parameters = {};

export var changed = false;

/***************************
 ******** Compute **********
 ***************************/

function apply_filters() {
    var mat = inputs.fetchCountMatrix();
    var disc = fetchDiscards();

    utils.freeCache(cache.matrix);
    cache.matrix = scran.filterCells(mat, disc);

    let block = inputs.fetchBlock();
    cache.blocked = (block !== null);

    if (cache.blocked) {
        let bcache = utils.allocateCachedArray(cache.matrix.numberOfColumns(), "Int32Array", cache, "block_buffer");

        let bcache_arr = bcache.array();
        let block_arr = block.array();
        let disc_arr = disc.array();
        let j = 0;
        for (let i = 0; i < block_arr.length; i++) {
            if (disc_arr[i] == 0) {
                bcache_arr[j] = block_arr[i];
                j++;
            }
        }
    }

    return;
}

export function compute(use_mito_default, mito_prefix, nmads) {
    changed = false;

    if (inputs.changed || use_mito_default !== parameters.use_mito_default || mito_prefix !== parameters.mito_prefix) {
        var mat = inputs.fetchCountMatrix();

        // TODO: add more choices.
        var nsubsets = 1;
        var subsets = utils.allocateCachedArray(mat.numberOfRows() * nsubsets, "Uint8Array", cache, "metrics_buffer");
        subsets.fill(0);

        // Finding the prefix.
        // TODO: use the guessed features to narrow the Ensembl/symbol search.
        var gene_info = inputs.fetchGenes();
        var sub_arr = subsets.array();
        for (const [key, val] of Object.entries(gene_info)) {
            if (use_mito_default) {
                val.forEach((x, i) => {
                    if (mito.symbol.has(x) || mito.ensembl.has(x)) {
                        sub_arr[i] = 1;
                    }
                });
            } else {
                var lower_mito = mito_prefix.toLowerCase();
                val.forEach((x, i) => {
                    if(x.toLowerCase().startsWith(lower_mito)) {
                        sub_arr[i] = 1;
                    }
                });
            }
        }

        utils.freeCache(cache.metrics);
        cache.metrics = scran.computePerCellQCMetrics(mat, subsets);

        parameters.use_mito_default = use_mito_default;
        parameters.mito_prefix = mito_prefix;
        changed = true;
    }

    if (changed || nmads !== parameters.nmads) {
        let block = inputs.fetchBlock();
        utils.freeCache(cache.filters);
        cache.filters = scran.computePerCellQCFilters(cache.metrics, { numberOfMADs: nmads, block: block });

        parameters.nmads = nmads;
        changed = true;
    }

    if (changed) {
        apply_filters();
        changed = true;
    }

    return;
}

/***************************
 ******** Results **********
 ***************************/

function format_metrics({ copy = true } = {}) {
    copy = (copy ? true : "view");
    return {
        sums: cache.metrics.sums({ copy: copy }),
        detected: cache.metrics.detected({ copy: copy }),
        proportion: cache.metrics.subsetProportions(0, { copy: copy })
    };
}

function format_thresholds({ copy = true } = {}) {
    copy = (copy ? true : "view");
    return {
        sums: cache.filters.thresholdsSums({ copy: copy }),
        detected: cache.filters.thresholdsDetected({ copy: copy }),
        proportion: cache.filters.thresholdsSubsetProportions(0, { copy: copy })
    }
}

export function results() {
    // This function should not do any Wasm allocations, so copy: false is safe.

    var data = {};
    var blocks = inputs.fetchBlockLevels();
    if (blocks === null) {
        blocks = [ "default" ];
        data["default"] = format_metrics();
    } else {
        let metrics = format_metrics({ copy: false });
        let bids = inputs.fetchBlock();
        let barray = bids.array();

        for (var b = 0; b < blocks.length; b++) {
            let current = {};
            for (const [key, val] of Object.entries(metrics)) {
                current[key] = val.array().filter((x, i) => barray[i] == b);
            }
            data[blocks[b]] = current;
        }
    }

    var thresholds = {};
    let listed = format_thresholds({ copy: false });
    for (var b = 0; b < blocks.length; b++) {
        let current = {};
        for (const [key, val] of Object.entries(listed)) {
            current[key] = val[b];
        }
        thresholds[blocks[b]] = current;
    }

    var ranges = {};
    for (var b = 0; b < blocks.length; b++) {
        let curranges = {};
        let curdata = data[blocks[b]];

        for (const [key, val] of Object.entries(curdata)) {
            var max = -Infinity, min = Infinity;
            val.forEach(function (x) {
                if (max < x) {
                    max = x;
                }
                if (min > x) {
                    min = x;
                }
            });
            curranges[key] = [min, max];
        }

        ranges[blocks[b]] = curranges;
    }

    let remaining = 0;
    if ("matrix" in cache) {
        remaining = cache.matrix.numberOfColumns();
    } else {
        fetchDiscards().array().forEach(x => {
            if (x == 0) {
                remaining++;
            }
        });
    }

    return { 
        "data": data, 
        "ranges": ranges,
        "thresholds": thresholds,
        "retained": remaining
    };
}

/*************************
 ******** Saving *********
 *************************/

export function serialize(handle) {
    let ghandle = handle.createGroup("quality_control");

    {
        let phandle = ghandle.createGroup("parameters"); 
        phandle.writeDataSet("use_mito_default", "Uint8", [], Number(parameters.use_mito_default));
        phandle.writeDataSet("mito_prefix", "String", [], parameters.mito_prefix);
        phandle.writeDataSet("nmads", "Float64", [], parameters.nmads);
    }

    {
        let rhandle = ghandle.createGroup("results"); 

        {
            let mhandle = rhandle.createGroup("metrics");
            let data = format_metrics({ copy: false });
            mhandle.writeDataSet("sums", "Float64", null, data.sums)
            mhandle.writeDataSet("detected", "Int32", null, data.detected);
            mhandle.writeDataSet("proportion", "Float64", null, data.proportion);
        }

        {
            let thandle = rhandle.createGroup("thresholds");
            let thresholds = format_thresholds({ copy: false });
            for (const x of [ "sums", "detected", "proportion" ]) {
                let current = thresholds[x];
                thandle.writeDataSet(x, "Float64", null, current);
            }
        }

        let disc = fetchDiscards();
        rhandle.writeDataSet("discards", "Uint8", null, disc);
    }
}

/**************************
 ******** Loading *********
 **************************/

class QCFiltersMimic {
    constructor(sums, detected, proportion, discards) {
        this.sums_ = sums;
        this.detected_ = detected;
        this.proportion_ = proportion;
        this.discards = scran.createUint8WasmArray(discards.length);
        this.discards.set(discards);
    }

    thresholdsSums({ copy }) {
        return utils.mimicGetter(this.sums_, copy);
    }

    thresholdsDetected({ copy }) {
        return utils.mimicGetter(this.detected_, copy);
    }

    thresholdsSubsetProportions(index, { copy }) {
        if (index != 0) {
            throw "only 'index = 0' is supported for mimics";
        }
        return utils.mimicGetter(this.proportion_, copy);
    }

    discardOverall({ copy }) {
        return utils.mimicGetter(this.discards, copy);
    }

    free() {
        this.discards.free();
    }
}

export function unserialize(handle) {
    let ghandle = handle.open("quality_control");

    {
        let phandle = ghandle.open("parameters"); 
        parameters = {
            use_mito_default: phandle.open("use_mito_default", { load: true }).values[0] > 0,
            mito_prefix: phandle.open("mito_prefix", { load: true }).values[0],
            nmads: phandle.open("nmads", { load: true }).values[0]
        }
    }

    utils.freeCache(cache.metrics);
    utils.freeCache(cache.filters);
    utils.freeCache(cache.block_buffer);
    cache = {};

    {
        let rhandle = ghandle.open("results");

        let mhandle = rhandle.open("metrics");
        let sums = mhandle.open("sums", { load: true }).values;

        cache.metrics = scran.emptyPerCellQCMetricsResults(sums.length, 1);
        cache.metrics.sums({ copy: false }).set(sums);

        let detected = mhandle.open("detected", { load: true }).values;
        cache.metrics.detected({ copy: false }).set(detected);
        let proportions = mhandle.open("proportion", { load: true }).values;
        cache.metrics.subsetProportions(0, { copy: false }).set(proportions);

        let thandle = rhandle.open("thresholds");
        let thresholds_sums = thandle.open("sums", { load: true }).values;
        let thresholds_detected = thandle.open("detected", { load: true }).values;
        let thresholds_proportion = thandle.open("proportion", { load: true }).values;

        let discards = rhandle.open("discards", { load: true }).values; 
        cache.filters = new QCFiltersMimic(
            thresholds_sums, 
            thresholds_detected,
            thresholds_proportion,
            discards
        );
    }

    return { ...parameters };
}

/***************************
 ******** Getters **********
 ***************************/

export function fetchSums({ unsafe = false } = {}) {
    // Unsafe, because we're returning a raw view into the Wasm heap,
    // which might be invalidated upon further allocations.
    return cache.metrics.sums({ copy: !unsafe });
}

export function fetchDiscards() {
    return cache.filters.discardOverall({ copy: "view" });
}

export function fetchFilteredMatrix() {
    if (!("matrix" in cache)) {
        apply_filters();
    }
    return cache.matrix;
}

export function fetchFilteredBlock() {
    if (!("blocked" in cache)) {
        apply_filters();
    }
    if (cache.blocked) {
        return cache.block_buffer;
    } else {
        return null;
    }
}
