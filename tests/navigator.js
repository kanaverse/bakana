import * as fs from "fs";
import * as utils from "./utils.js";
import * as bakana from "../src/index.js";

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
