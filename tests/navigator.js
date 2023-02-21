import * as fs from "fs";
import * as utils from "./utils.js";
import * as bakana from "../src/index.js";

export const baseDirectory = "miscellaneous/from-tests";

export async function obtainLocalTestPath() {
    let h5path = baseDirectory + "/H5AD";
    if (fs.existsSync(h5path)) {
        return "H5AD";
    }

    let fpath = "files/datasets/zeisel-brain.h5ad";
    let files = { 
        default: new bakana.H5adDataset(fpath)
    };

    let altpath = "ArtifactDB-testing";
    let fullalt = baseDirectory + "/" + altpath;
    if (fs.existsSync(fullalt)) {
        return altpath;
    }

    let state = await bakana.createAnalysis();
    try {
        let params = utils.baseParams();
        await bakana.runAnalysis(state, files, params);
        await bakana.saveSingleCellExperiment(state, altpath, { directory: baseDirectory });
    } finally {
        await bakana.freeAnalysis(state);
    }

    return altpath;
}

export class LocalProjectDirectoryNavigator {
    #directory;

    constructor(d) {
        this.#directory = d;
    }

    serialize() {
        let enc = new TextEncoder;
        let contents = enc.encode(this.#directory);
        let f = new bakana.SimpleFile(contents, { name: this.#directory });
        return [{ type: "directory", file: f }];
    }

    metadata(p) {
        let fullpath = this.#directory + "/" + p;

        while (1) {
            if (!fullpath.endsWith(".json")) { 
                fullpath += ".json";
            }
            let contents = fs.readFileSync(fullpath);

            let dec = new TextDecoder;
            let json = dec.decode(contents);
            let values = JSON.parse(json);

            if (values["$schema"].startsWith("redirection/")){
                fullpath = this.#directory + "/" + values.redirection.targets[0].location;
            } else {
                return values;
            }
        }
    }

    file(p) {
        return this.#directory + "/" + p;
    }
}
