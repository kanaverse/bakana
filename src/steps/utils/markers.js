import * as scran from "scran.js";
import * as wa from "wasmarrays.js";

// SOFT-DEPRECATED BEGIN.
export const summaries2int = { "min": 0, "mean": 1, "min_rank": 4 };

export function unserializeGroupStats(handle, permuter, { no_summaries = false, compute_auc = true } = {}) {
    let output = {};
    for (const x of [ "means", "detected" ]) {
        output[x] = permuter(handle.open(x, { load: true }).values);
    }

    for (const i of [ "lfc", "delta_detected", "auc", "cohen" ]) {
        if (i == "auc" && !compute_auc) {
            continue;
        }

        if (no_summaries) {
            output[i] = handle.open(i, { load: true }).values;
        } else {
            let rhandle = handle.open(i);
            let current = {};
            for (const j of Object.keys(rhandle.children)) {
                current[j] = permuter(rhandle.open(j, { load: true }).values);
            }
            output[i] = current;
        }
    }

    return output;
}

export function fillGroupStats(object, i, vals) {
    object.means(i, { copy: false }).set(vals.means);
    object.detected(i, { copy: false }).set(vals.detected);

    for (const [s, v] of Object.entries(vals.cohen)) {
        object.cohensD(i, { summary: summaries2int[s], copy: false }).set(v);
    }

    for (const [s, v] of Object.entries(vals.lfc)) {
        object.deltaMean(i, { summary: summaries2int[s], copy: false }).set(v);
    }

    for (const [s, v] of Object.entries(vals.delta_detected)) {
        object.deltaDetected(i, { summary: summaries2int[s], copy: false }).set(v);
    }

    if ("auc" in vals) {
        for (const [s, v] of Object.entries(vals.auc)) {
            object.auc(i, { summary: summaries2int[s], copy: false }).set(v);
        }
    }
}

export function formatMarkerResults(results, group, rankEffect) {
    if (!rankEffect || rankEffect === undefined) {
        rankEffect = "cohen-min-rank";
    }

    var ordering;
    {
        // Choosing the ranking statistic. Do NOT do any Wasm allocations
        // until 'ranking' is fully consumed!
        let increasing = false;
        let summary = "mean";
        if (rankEffect.match(/-min$/)) {
            summary = "min";
        } else if (rankEffect.match(/-min-rank$/)) {
            summary = "min-rank";
            increasing = true;
        }

        let ranking;
        if (rankEffect.match(/^cohen-/)) {
            ranking = results.cohensD(group, { summary: summary, copy: false });
        } else if (rankEffect.match(/^auc-/)) {
            ranking = results.auc(group, { summary: summary, copy: false });
        } else if (rankEffect.match(/^lfc-/)) {
            ranking = results.deltaMean(group, { summary: summary, copy: false });
        } else if (rankEffect.match(/^delta-d-/)) {
            ranking = results.deltaDetected(group, { summary: summary, copy: false });
        } else {
            throw "unknown rank type '" + rankEffect + "'";
        }
  
        // Computing the ordering based on the ranking statistic.
        ordering = new Int32Array(ranking.length);
        for (var i = 0; i < ordering.length; i++) {
            ordering[i] = i;
        }
        if (increasing) {
            ordering.sort((f, s) => (ranking[f] - ranking[s]));
        } else {
            ordering.sort((f, s) => (ranking[s] - ranking[f]));
        }
    }
  
    // Apply that ordering to each statistic of interest.
    var reorder = function(stats) {
        var thing = new Float64Array(stats.length);
        for (var i = 0; i < ordering.length; i++) {
            thing[i] = stats[ordering[i]];
        }
        return thing;
    };
  
    var stat_detected = reorder(results.detected(group, { copy: false }));
    var stat_mean = reorder(results.mean(group, { copy: false }));
    var stat_lfc = reorder(results.deltaMean(group, { summary: 1, copy: false }));
    var stat_delta_d = reorder(results.deltaDetected(group, { summary: 1, copy: false }));

    return {
        "ordering": ordering,
        "means": stat_mean,
        "detected": stat_detected,
        "lfc": stat_lfc,
        "delta_detected": stat_delta_d
    };
}
// SOFT-DEPRECATED END.

export function locateVersusCache(left, right, cache) {
    let left_small = left < right;

    let bigg = (left_small ? right : left);
    if (!(bigg in cache)) {
        cache[bigg] = {};
    }
    let biggversus = cache[bigg];

    let smal = (left_small ? left : right); 
    let rerun = !(smal in biggversus);
    if (rerun) {
        biggversus[smal] = {};
    }

    return { 
        cached: biggversus[smal],
        run: rerun,
        left_small: left_small
    };
}

export function freeVersusResults(cache) {
    if (cache) {
        for (const v of Object.values(cache)) {
            for (const v2 of Object.values(v)) {
                for (const m of Object.values(v2)) {
                    scran.free(m);
                }
            }
        }
        for (const k of Object.keys(cache)) {
            delete cache[k];
        }
    }
}

export function computeVersusResults(matrices, clusters, block, keep, cache, lfc_threshold, compute_auc) {
    let new_block = null;
    if (block !== null) {
        new_block = wa.subsetWasmArray(block, keep);
        scran.dropUnusedLevels(new_block);
    }

    for (const modality of matrices.available()) {
        let modmat = matrices.get(modality);
        let sub;
        try {
            sub = scran.subsetColumns(modmat, keep);
            cache[modality] = scran.scoreMarkers(sub, clusters, { block: new_block, threshold: lfc_threshold, computeAuc: compute_auc });
        } finally {
            scran.free(sub);
        }
    }
}
