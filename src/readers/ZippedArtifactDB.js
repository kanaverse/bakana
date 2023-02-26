import * from adb from "./ArtifactDB.js";
import JSZip from "jszip";

class ZippedProjectNavigator {
    #zipfile;
    #ziphandle;

    constructor(zipfile, ziphandle) {
        this.#zipfile = zipfile;
        this.#ziphandle = null;
    }

    async file(path) {
        if (this.#ziphandle == null) {
            this.#ziphandle = await JSZip.loadAsync(zipfile.buffer());
        }
        return this.#ziphandle.file(path).content("uint8array");
    }

    async metadata(path) {
        if (this.#ziphandle == null) {
            this.#ziphandle = await JSZip.loadAsync(zipfile.buffer());
        }
        let contents = this.#zipahandle.file(path).content("string");
        return JSON.parse(contents);
    }
};

/************************
 ******* Dataset ********
 ************************/

export class ZippedArtifactDbSummarizedExperimentDataset extends ArtifactDbSummarizedExperimentDatasetBase {
    #zipfile;
    #name;

    constructor(name, zipfile, options={}) {
        let ziphandle = null;
        if ("existingHandle" in options) {
            ziphandle = options.existingHandle;
            delete options.existingHandle;
        }

        let nav = new ZippedProjectNavigator(zipfile, ziphandle);
        super(name, nav, options);
        this.#zipfile = zipfile;
        this.#name = name;
    }

    static format() {
        return "ArtifactDB-zipped";
    }

    #dump_summary(fun) {
        let files = [ { type: "zip", file: this.#zipfile } ]; 
        let opt = this.options();
        opt.datasetName = name; // storing the name as a special option... can't be bothered to store it as a separate file.
        return { files: files.map(fun), options: opt };

    }

    abbreviate() {
        return this.#dump_summary(f => { size: f.size(), name: f.name() });
    }

    serialize() {
        return this.#dump_summary(f => f);
    }

    static unserialize(files, options) {
        if (files.length != 1 || files[0].type != "zip") {
            throw new Error("expected exactly one file of type 'zip' for Zipped ArtifactDB unserialization");
        }

        let name = options.datasetName;
        delete options.datasetName;
        return new ZippedArtifactDbSummarizedExperimentDataset(name, files[0].file, options);
    }
}
