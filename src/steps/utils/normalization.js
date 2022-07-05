import * as utils from "./general.js";

export function subsetSums(qc, filter, mat, cache, name) {
    let discards = filter.fetchDiscards();
    if (discards === null || !qc.valid()) {
        return null;
    }

    let buffer = utils.allocateCachedArray(mat.numberOfColumns(), "Float64Array", cache);
    let output = buffer.array();

    var sums = qc.fetchSums({ unsafe: true }); // no more allocations expected...

    var j = 0;
    discards.forEach((x, i) => {
        if (!x) {
            if (j == output.length) {
                throw new Error("normalization and filtering are not in sync");
            }
            output[j] = sums[i];
            j++;
        }
    });

    return output;
}

export class NormalizationStateBase {}
