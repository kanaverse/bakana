import * as scran from "scran.js";

export const summaries2int = { "min": 0, "mean": 1, "min_rank": 4 };
export const int2summaries = { 0: "min", 1: "mean", 4: "min_rank" };

export function serializeGroupStats(ihandle, obj, group, { no_summaries = false } = {}) {
    for (const x of [ "means", "detected" ]) {
        let y= obj[x](group, { copy: "hdf5" });
        ihandle.writeDataSet(x, "Float64", null, y);
    }

    for (const i of [ "lfc", "delta_detected", "auc", "cohen" ]) {
        let i0 = i;
        if (i == "delta_detected") {
            i0 = "deltaDetected";
        }

        let extractor = (index) => obj[i0](group, { summary: index, copy: "hdf5" });
        if (no_summaries) {
            let y = extractor(summaries2int["mean"]);
            ihandle.writeDataSet(i, "Float64", null, y);
        } else {
            let curhandle = ihandle.createGroup(i);
            for (const [j, k] of Object.entries(summaries2int)) {
                let y = extractor(k);
                curhandle.writeDataSet(j, "Float64", null, y);
            }
        }
    }
}

export function unserializeGroupStats(handle, permuter, { no_summaries = false } = {}) {
    let output = {};
    for (const x of [ "means", "detected" ]) {
        output[x] = handle.open(x, { load: true }).values;
        permuter(output[x]);
    }

    for (const i of [ "lfc", "delta_detected", "auc", "cohen" ]) {
        if (no_summaries) {
            output[i] = handle.open(i, { load: true }).values;
        } else {
            let rhandle = handle.open(i);
            let current = {};
            for (const j of Object.keys(summaries2int)) {
                current[j] = rhandle.open(j, { load: true }).values;
                permuter(current[j]);
            }
            output[i] = current;
        }
    }

    return output;
}

/**
 * Report marker results for a given group or cluster, ordered so that the strongest candidate markers appear first.
 *
 * @param {ScoreMarkersResults} results - The marker results object generated by the `scoreMarkers` function in **scran.js**.
 * @param {number} group - Integer specifying the group or cluster of interest.
 * Any number can be used if it was part of the `groups` passed to `scoreMarkers`.
 * @param {string} rankEffect - Summarized effect size to use for ranking markers.
 * This should follow the format of `<effect>-<summary>` where `<effect>` may be `lfc`, `cohen`, `auc` or `delta_detected`,
 * and `<summary>` may be `min`, `mean` or `min-rank`.
 *
 * @return An object containing the marker statistics for the selection, sorted by the specified effect and summary size from `rankEffect`.
 * This contains:
 *   - `means`: a Float64Array of length equal to the number of genes, containing the mean expression within the selection.
 *   - `detected`: a Float64Array of length equal to the number of genes, containing the proportion of cells with detected expression inside the selection.
 *   - `lfc`: a Float64Array of length equal to the number of genes, containing the log-fold changes for the comparison between cells inside and outside the selection.
 *   - `delta_detected`: a Float64Array of length equal to the number of genes, containing the difference in the detected proportions between cells inside and outside the selection.
 */
export function formatMarkerResults(results, group, rankEffect) {
    if (!rankEffect || rankEffect === undefined) {
        rankEffect = "cohen-min-rank";
    }

    var ordering;
    {
        // Choosing the ranking statistic. Do NOT do any Wasm allocations
        // until 'ranking' is fully consumed!
        let ranking;
        let increasing = false;
      
        let index = 1;
        if (rankEffect.match(/-min$/)) {
            index = 0;
        } else if (rankEffect.match(/-min-rank$/)) {
            increasing = true;
            index = 4;
        }

        if (rankEffect.match(/^cohen-/)) {
            ranking = results.cohen(group, { summary: index, copy: false });
        } else if (rankEffect.match(/^auc-/)) {
            ranking = results.auc(group, { summary: index, copy: false });
        } else if (rankEffect.match(/^lfc-/)) {
            ranking = results.lfc(group, { summary: index, copy: false });
        } else if (rankEffect.match(/^delta-d-/)) {
            ranking = results.deltaDetected(group, { summary: index, copy: false });
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
    var stat_mean = reorder(results.means(group, { copy: false }));
    var stat_lfc = reorder(results.lfc(group, { summary: 1, copy: false }));
    var stat_delta_d = reorder(results.deltaDetected(group, { summary: 1, copy: false }));

    return {
        "ordering": ordering,
        "means": stat_mean,
        "detected": stat_detected,
        "lfc": stat_lfc,
        "delta_detected": stat_delta_d
    };
}

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
        delete cache.versus;
    }
}

export function computeVersusResults(matrices, clusters, block, keep, cache) {
    let new_block = null;
    if (block !== null) {
        new_block = scran.subsetBlock(block, keep);
        scran.dropUnusedBlocks(new_block);
    }

    for (const modality of matrices.available()) {
        let modmat = matrices.get(modality);
        let sub;
        try {
            sub = scran.subsetColumns(modmat, keep);
            cache[modality] = scran.scoreMarkers(sub, clusters, { block: new_block });
        } finally {
            scran.free(sub);
        }
    }
}
