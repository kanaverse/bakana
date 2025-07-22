import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";
import * as bioc from "bioconductor";
import * as fs from "fs";
import jszip from "jszip";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

export class LocalDirectoryNavigator {
    #host;

    constructor(d) {
        this.#host = d;
    }

    #full(path) {
        if (!path.startsWith("/")) {
            path = "/" + path;
        }
        return this.#host + path;
    }

    get(path, asBuffer) {
        let full = this.#full(path);
        if (asBuffer) {
            let contents = fs.readFileSync(full);
            return new Uint8Array(contents);
        } else {
            return full;
        }
    }

    exists(path) {
        return fs.existsSync(this.#full(path));
    }

    clean(getPath) {}
}

class LocalAlabasterDataset extends bakana.AbstractAlabasterDataset {
    #dir;

    constructor(dir, options={}) {
        super(new LocalDirectoryNavigator(dir), options);
        this.#dir = dir;
    }

    static format() {
        return "alabaster-local-directory";
    }

    static unserialize(files, options) {
        let dec = new TextDecoder;
        let payload = files[0].file.buffer();
        let info = JSON.stringify(dec.decode(payload));
        return new LocalAlabasterDataset(info.directory, options);
    }

    abbreviate() {
        return {
            directory: this.#dir,
            options: this.options()
        };
    }

    serialize() {
        let enc = new TextEncoder;
        let config = { directory: this.#dir };
        let buffer = enc.encode(JSON.stringify(config));
        return {
            files: [ { file: new bakana.SimpleFile(buffer, { name: "config.json" }), type: "config" } ],
            options: this.options()
        }
    }
}

/***********************************************/

test("AlabasterAbstractDataset for simple datasets", async () => {
    let stripped_ds = new LocalAlabasterDataset("files/alabaster/zeisel-brain-stripped");

    let summ = await stripped_ds.summary({ cache: true });
    expect(summ.modality_features[""] instanceof bioc.DataFrame).toBe(true);
    expect(summ.modality_features[""].numberOfRows()).toBeGreaterThan(0);
    expect(summ.cells instanceof bioc.DataFrame).toBe(true);

    let preview = await stripped_ds.previewPrimaryIds({ cache: true });
    expect("RNA" in preview).toBe(true);
    expect(preview.RNA.length).toBeGreaterThan(0);

    let loaded = await stripped_ds.load({ cache: true });
    expect(loaded.matrix.numberOfColumns()).toEqual(loaded.cells.numberOfRows());
    expect(loaded.matrix.available()).toEqual(["RNA"]);
    expect(loaded.matrix.get("RNA").numberOfRows()).toEqual(loaded.features["RNA"].numberOfRows());
    expect(loaded.primary_ids["RNA"]).toEqual(loaded.features["RNA"].rowNames());

    stripped_ds.clear();
})

//TEst("alabaster summary and loading works with multiple modalities", async () => {
//    let files = { super: new LocalAlabasterDataset(target_adt, nav.baseDirectory) };
//
//    let summ = await files.super.summary();
//    expect(Object.keys(summ.modality_features).sort()).toEqual(["", "ADT"]);
//    expect(summ.modality_assay_names).toEqual({ "": ["counts", "logcounts"], "ADT": [ "counts", "logcounts" ] });
//
//    // Trying with a name for the experiment.
//    files.super.setOptions({ adtExperiment: "ADT" });
//    {
//        let everything = await files.super.load();
//
//        let nRna = everything.features["RNA"].numberOfRows();
//        expect(nRna).toBeGreaterThan(0);
//        expect(everything.matrix.get("RNA").numberOfRows()).toEqual(nRna);
//
//        let nAdt = everything.features["ADT"].numberOfRows();
//        expect(nAdt).toBeGreaterThan(0);
//        expect(nAdt).toBeLessThan(nRna);
//        expect(everything.matrix.get("ADT").numberOfRows()).toEqual(nAdt);
//
//        expect(everything.matrix.numberOfColumns()).toEqual(everything.cells.numberOfRows());
//        everything.matrix.free();
//    }
//
//    // Trying with a numeric index for the ADT experiment.
//    files.super.setOptions({ adtExperiment: 0 });
//    {
//        let everything = await files.super.load();
//
//        let nRna = everything.features["RNA"].numberOfRows();
//        expect(nRna).toBeGreaterThan(0);
//        expect(everything.matrix.get("RNA").numberOfRows()).toEqual(nRna);
//
//        let nAdt = everything.features["ADT"].numberOfRows();
//        expect(nAdt).toBeGreaterThan(0);
//        expect(nAdt).toBeLessThan(nRna);
//        expect(everything.matrix.get("ADT").numberOfRows()).toEqual(nAdt);
//
//        expect(everything.matrix.numberOfColumns()).toEqual(everything.cells.numberOfRows());
//        everything.matrix.free();
//    }
//})
