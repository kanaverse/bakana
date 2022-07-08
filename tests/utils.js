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

