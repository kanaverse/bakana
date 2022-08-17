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
    let fpath = path.basename(decodeURIComponent(url));
    let obj = fs.readFileSync("files/references/" + fpath);
    return obj.buffer.slice(obj.byteOffset, obj.byteOffset + obj.byteLength);
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

    if (referenceSubset) {
        if (!is_same(ids.slice().sort(), loadedIds.slice().sort())) {
            throw new Error("reference and loaded identities should have identical elements");
        }
        let mapping = {};
        ids.forEach((x, i) => {
            mapping[x] = i;
        });
        loadedIds = loadedIds.map(x => mapping[x]);
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
        loadedIds.forEach((x, i) => {
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
    loadedIds.forEach((x, i) => {
        converted[i] = names[x];
    });

    if (!is_same(converted, loadedNames)) {
        throw new Error("loaded matrix's name reorganization is not consistent");
    }
}

