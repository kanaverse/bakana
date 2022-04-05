import * as path from "path";
import * as fs from "fs";
import * as bakana from "../src/index.js";

export var baseParams = bakana.analysisDefaults();
baseParams.cell_labelling = {
    mouse_references: [ "ImmGen" ],
    human_references: [ "BlueprintEncode" ]
};

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

