export function splitMetricsByBlock(metrics, blockLevels, blockIds) {
    var output = {};
    var blocks = blockIds.array();
    for (var b = 0; b < blockLevels.length; b++) {
        let current = {};
        for (const [key, val] of Object.entries(metrics)) {
            current[key] = val.array().filter((x, i) => blocks[i] == b);
        }
        output[blockLevels[b]] = current;
    }
    return output;
}

export function splitThresholdsByBlock(thresholds, blockLevels) {
    var output = {};
    for (const x of blockLevels) {
        output[x] = {};
    }

    for (const [key, val] of Object.entries(thresholds)) {
        let arr = val.array();
        for (var b = 0; b < blockLevels.length; b++) {
            output[blockLevels[b]][key] = arr[b];
        }
    }

    return output;
}
