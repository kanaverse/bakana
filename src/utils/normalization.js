export function subsetSums(qc, filter, output) {
    var discards = filter.fetchDiscards();
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
}
