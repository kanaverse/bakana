import * as fs from "fs";
import * as utils from "./utils.js";
import * as bakana from "../src/index.js";
import JSZip from "jszip";

export const baseDirectory = "miscellaneous/from-tests";

export function pathExists(chosen) {
    let path = baseDirectory + "/" + chosen;
    if (fs.existsSync(path)) {
        return chosen;
    } else {
        return null;
    }
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

function list_json(directory, prefix, files) {
    const target = (prefix == null ? directory : directory + "/" + prefix);
    return fs.readdirSync(target).forEach(f => {
        const host = (prefix == null ? f : prefix + "/" + f);
        const full_path = directory + "/" + host;
        if (fs.statSync(full_path).isDirectory()) {
            return list_json(directory, host, files);
        } else if (f.endsWith(".json")) {
            files.push(host);
        }
    });
}

export function zipDirectory(directory, children) {
    // Finding all JSON documents.
    let all_files = [];
    for (const x of children) {
        if (fs.statSync(directory + "/" + x).isDirectory()) {
            list_json(directory, x, all_files);
        } else {
            all_files.push(x);
        }
    }

    // Mimicking the output of saveSingleCellExperiment.
    let prep_files = [];
    const dec = new TextDecoder;
    for (const f of all_files) {
        let buffer = fs.readFileSync(directory + "/" + f);
        let str = dec.decode(buffer);
        let parsed = JSON.parse(str);

        if (parsed.path == f || parsed["$schema"].startsWith("redirection/")) {
            prep_files.push({ metadata: parsed });
        } else {
            // Sometimes we report a Uint8Array, sometimes a file path,
            // to make sure we get test coverage of both possibilities in zipFiles.
            let obj = { metadata: parsed, contents: directory + "/" + parsed.path };
            if (Math.random() > 0.5) {
                obj.contents = fs.readFileSync(obj.contents);
            }
            prep_files.push(obj);
        }
    }

    const enc = new TextEncoder;
    let ex = {
        path: "README.md",
        contents: enc.encode("Hi! Read me please!")
    };

    return bakana.zipFiles(prep_files, { extras: [ex] });
}
