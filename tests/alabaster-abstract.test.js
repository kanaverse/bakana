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
    expect(summ.cells.numberOfRows()).toBeGreaterThan(0);

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

test("AlabasterAbstractDataset for dense matrices", async () => {
    let dense_ds = new LocalAlabasterDataset("files/alabaster/zeisel-brain-dense");
    dense_ds.setOptions({ rnaExperiment: "endogenous", adtExperiment: null, crisprExperiment: null });

    let loaded = await dense_ds.load({ cache: true });
    expect(loaded.matrix.available()).toEqual(["RNA"]);
    const mat = loaded.matrix.get("RNA");
    const feat = loaded.features["RNA"];
    expect(mat.numberOfRows()).toEqual(feat.numberOfRows());
    expect(mat.isSparse()).toBe(true); // we force it to be sparse, because that's how it is.
    expect(loaded.primary_ids["RNA"]).toEqual(feat.rowNames());
})

/***********************************************/

class LocalAlabasterResult extends bakana.AbstractAlabasterResult {
    #dir;
    constructor(dir, options={}) {
        super(new LocalDirectoryNavigator(dir), options);
        this.#dir = dir;
    }
}

test("AlabasterAbstractResult behaves correctly ", async () => {
    let res = new LocalAlabasterResult("files/alabaster/zeisel-brain-with-results");

    let summ = await res.summary({ cache: true });
    expect(Object.keys(summ.modality_features)).toEqual(["endogenous", "ERCC", "repeat"]);
    for (const mod of Object.values(summ.modality_features)) {
        expect(mod instanceof bioc.DataFrame).toBe(true);
        expect(mod.numberOfRows()).toBeGreaterThan(0);
    }

    expect(summ.cells instanceof bioc.DataFrame).toBe(true);
    expect(summ.cells.numberOfRows()).toBeGreaterThan(0);
    expect(summ.reduced_dimension_names).toEqual(["PCA", "TSNE", "UMAP"]);

    expect(summ.modality_assay_names["endogenous"]).toEqual(["counts", "logcounts"]);
    expect(summ.modality_assay_names["ERCC"]).toEqual(["counts"]);
    expect(summ.modality_assay_names["repeat"]).toEqual(["counts"]);

    // Checking that we can load the log counts.
    res.setOptions({
        primaryAssay: { endogenous: "logcounts", "ERCC": 0 },
        reducedDimensionNames: [ "TSNE", "UMAP" ]
    });

    let loaded = await res.load({ cache: true });
    expect(loaded.matrix.available()).toEqual(["endogenous", "ERCC"]);
    expect(loaded.matrix.numberOfColumns()).toEqual(loaded.cells.numberOfRows());
    expect(loaded.matrix.get("endogenous").numberOfRows()).toEqual(loaded.features["endogenous"].numberOfRows());
    expect(loaded.matrix.get("ERCC").numberOfRows()).toEqual(loaded.features["ERCC"].numberOfRows());

    expect(Object.keys(loaded.reduced_dimensions)).toEqual(["TSNE", "UMAP"]);
    expect(loaded.reduced_dimensions["TSNE"].length).toEqual(2);
    expect(loaded.reduced_dimensions["TSNE"][0].length).toEqual(loaded.matrix.numberOfColumns());
    expect(loaded.reduced_dimensions["UMAP"].length).toEqual(2);
    expect(loaded.reduced_dimensions["UMAP"][1].length).toEqual(loaded.matrix.numberOfColumns());
})
