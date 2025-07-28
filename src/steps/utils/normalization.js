export function subsetSums(qc, filter, mat) {
    let output = new Float64Array(mat.numberOfColumns());
    let keep = filter.fetchKeep();

    // unsafe, so no more Wasm allocations past this point. 
    let sums = qc.fetchMetrics().sum({ copy: false }); 

    if (keep === null) {
        output.set(sums);
    } else {
        var j = 0;
        keep.forEach((x, i) => {
            if (x) {
                if (j == output.length) {
                    throw new Error("normalization and filtering are not in sync");
                }
                output[j] = sums[i];
                j++;
            }
        });
        if (j !== output.length) {
            throw new Error("normalization and filtering are not in sync");
        }
    }

    return output;
}
