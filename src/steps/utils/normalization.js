import * as utils from "./general.js";

export function subsetSums(qc, filter, mat, cache, name) {
    let output = utils.allocateCachedArray(mat.numberOfColumns(), "Float64Array", cache, name);
    let discards = filter.fetchDiscards();

    // unsafe, so no more Wasm allocations past this point. 
    let sums = qc.fetchMetrics().sums({ copy: false }); 

    if (discards == null) {
        output.set(sums);
    } else {
        let oarr = output.array();
        var j = 0;
        discards.forEach((x, i) => {
            if (!x) {
                if (j == output.length) {
                    throw new Error("normalization and filtering are not in sync");
                }
                oarr[j] = sums[i];
                j++;
            }
        });
        if (j !== output.length) {
            throw new Error("normalization and filtering are not in sync");
        }
    }

    return output;
}
