import * as path from "path";
import * as fs from "fs";
import * as bakana from "../src/index.js";
import * as valkana from "valkana";

export async function initializeAll() {
    await bakana.initialize({ localFile: true });
    await valkana.initialize({ localFile: true });
}

export function validateState(path, embedded = true) {
    valkana.validateState(path, embedded, bakana.kanaFormatVersion);
}

export function baseParams() {
    let output = bakana.analysisDefaults();

    // Cut down on the work.
    output.pca.num_pcs = 10;

    // Avoid getting held up by pointless iterations.
    output.tsne.iterations = 10;
    output.umap.num_epochs = 10;

    // Actually do something.
    output.cell_labelling = {
        mouse_references: [ "ImmGen" ],
        human_references: [ "BlueprintEncode" ]
    };
    return output;
}

bakana.setCellLabellingDownload(url => {
    let fpath = path.basename(url);
    return fs.readFileSync("files/references/" + fpath).slice(); // Mimic a buffer from fetch().
});

export function mockOffsets(paths) {
    let offsets = {};
    let sofar = 0;
    for (const p of paths) {
        offsets[sofar] = p;
        sofar += fs.statSync(p).size;
    }
    return offsets;
}

function is_same(left, right) {
    if (left.length != right.length) {
        throw new Error("mismatch in array lengths");
    }
    for (var i = 0; i < left.length; i++) {
        if (left[i] != right[i]) {
            return false;
        }
    }
    return true;
}

export function checkReorganization(matrix, ids, names, loadedMatrix, loadedIds, loadedNames, { mustDiffer = true, referenceSubset = false } = {}) {
    if (!referenceSubset && ids !== null) { 
        throw new Error("reference matrix should not be reorganized");
    } else if (referenceSubset && ids === null) {
        throw new Error("subsetted reference matrix should be reorganized");
    }
    if (loadedIds === null) {
        throw new Error("loaded matrix should be reorganized");
    }

    let NR = matrix.numberOfRows();
    let NC = matrix.numberOfColumns();
    if (loadedMatrix.numberOfRows() != NR || NC != loadedMatrix.numberOfColumns()) {
        throw new Error("loaded and reference matrix have different dimensions");
    }

    if (mustDiffer) {
        let same = true;
        for (var i = 0; i < loadedIds.length; i++) {
            if (loadedIds[i] != i) {
                same = false;
                break;
            }
        }
        if (same) {
            throw new Error("identities should not be trivial after reorganization");
        }
    }

    let sorted_ids = loadedIds.slice().sort((a, b) => a - b);
    let permuter;
    if (referenceSubset) {
        // Assume reference IDs are already sorted when referenceSubset = true.
        if (!is_same(ids, sorted_ids)) {
            throw new Error("reference and loaded identities should have identical elements");
        }

        // Creating a mapping to permute the reference to the loaded order.
        let mapping = {};
        ids.forEach((x, i) => {
            mapping[x] = i;
        });
        permuter = Array.from(loadedIds).map(x => mapping[x]);
    } else {
        // Should contain everything at least once.
        for (var i = 0; i < sorted_ids.length; i++) {
            if (sorted_ids[i] != i) {
                throw new error("loaded identities should contain all consecutive integers");
            }
        }
        permuter = Array.from(loadedIds);
    }

    // Checking that the reorganization matches up with the reference for every 100th column.
    let at_least_one_difference = false;

    for (var c = 0; c < NC; c += 100) {
        let loaded_first = loadedMatrix.column(c);
        let reference_first = matrix.column(c);
        if (mustDiffer && !is_same(loaded_first, reference_first)) {
            at_least_one_difference = true;
        }

        let converted = new reference_first.constructor(NR);
        permuter.forEach((x, i) => {
            converted[i] = reference_first[x];
        });
        if (!is_same(loaded_first, converted)) {
            throw new Error("loaded matrix's reorganization is not consistent"); 
        }
    }

    if (mustDiffer && !at_least_one_difference) {
        throw new Error("loaded and reference columns should have some kind of difference");
    }

    // Now checking the identifiers.
    if (mustDiffer && is_same(names, loadedNames)) {
        throw new Error("loaded and reference names should not have been the same");
    }

    let converted = new names.constructor(NR);
    permuter.forEach((x, i) => {
        converted[i] = names[x];
    });

    if (!is_same(converted, loadedNames)) {
        throw new Error("loaded matrix's name reorganization is not consistent");
    }
}

