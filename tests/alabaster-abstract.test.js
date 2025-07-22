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
    expect(Object.keys(summ.modality_features)).toEqual([""]);
    expect(summ.modality_features[""] instanceof bioc.DataFrame).toBe(true);
    expect(summ.modality_features[""].numberOfRows()).toBeGreaterThan(0);
    expect(summ.cells instanceof bioc.DataFrame).toBe(true);

    let preview = await stripped_ds.previewPrimaryIds({ cache: true });
    expect(Object.keys(preview)).toEqual(["RNA"]);
    expect(preview.RNA.length).toBeGreaterThan(0);

    let loaded = await stripped_ds.load({ cache: true });
    expect(loaded.matrix.numberOfColumns()).toEqual(loaded.cells.numberOfRows());
    expect(loaded.matrix.available()).toEqual(["RNA"]);
    const mat = loaded.matrix.get("RNA");
    const feat = loaded.features["RNA"];
    expect(mat.numberOfRows()).toEqual(feat.numberOfRows());
    expect(mat.isSparse()).toBe(true);
    expect(loaded.primary_ids["RNA"]).toEqual(feat.rowNames());

    stripped_ds.clear();
})

test("AlabasterAbstractDataset for complex datasets", async () => {
    let full_ds = new LocalAlabasterDataset("files/alabaster/zeisel-brain");

    let summ = await full_ds.summary({ cache: true });
    expect(Object.keys(summ.modality_features)).toEqual(["endogenous", "ERCC", "repeat"]);
    for (const mod of Object.values(summ.modality_features)) {
        expect(mod instanceof bioc.DataFrame).toBe(true);
        expect(mod.numberOfRows()).toBeGreaterThan(0);
    }
    expect(summ.cells instanceof bioc.DataFrame).toBe(true);

    full_ds.setOptions({
        rnaExperiment: "endogenous",
        adtExperiment: "ERCC",
        crisprExperiment: "repeat",
        rnaCountAssay: "counts",
        adtCountAssay: "counts",
        crisprCountAssay: "counts"
    });
    let preview = await full_ds.previewPrimaryIds({ cache: true });
    expect(Object.keys(preview)).toEqual(["RNA", "ADT", "CRISPR"]);
    for (const mod of Object.values(preview)) {
        expect(mod.length).toBeGreaterThan(0);
    }

    let loaded = await full_ds.load({ cache: true });
    expect(loaded.matrix.numberOfColumns()).toEqual(loaded.cells.numberOfRows());
    expect(loaded.matrix.available()).toEqual(["RNA", "ADT", "CRISPR"]);
    for (const mod of loaded.matrix.available()) {
        const mat = loaded.matrix.get(mod);
        const feat = loaded.features[mod];
        expect(mat.numberOfRows()).toEqual(feat.numberOfRows());
        expect(mat.isSparse()).toEqual(true);
        expect(loaded.primary_ids[mod]).toEqual(feat.rowNames());
    }

    full_ds.clear();
})

test("AlabasterAbstractDataset for complex matrices", async () => {
    let dense_ds = new LocalAlabasterDataset("files/alabaster/zeisel-brain-dense");
    dense_ds.setOptions({ rnaExperiment: "endogenous", adtExperiment: null, crisprExperiment: null });

    let loaded = await dense_ds.load({ cache: true });
    expect(loaded.matrix.available()).toEqual(["RNA"]);
    const mat = loaded.matrix.get("RNA");
    const feat = loaded.features["RNA"];
    expect(mat.numberOfRows()).toEqual(feat.numberOfRows());
    expect(mat.isSparse()).toBe(false);
    expect(loaded.primary_ids["RNA"]).toEqual(feat.rowNames());
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
