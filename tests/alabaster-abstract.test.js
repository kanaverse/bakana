import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";
import * as bioc from "bioconductor";
import * as fs from "fs";

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
        let info = JSON.parse(dec.decode(payload));
        let output = new LocalAlabasterDataset(info.directory);
        output.setOptions(options);
        return output;
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
    let stripped_ds = new LocalAlabasterDataset("files/datasets/alabaster/zeisel-brain-stripped");
    await utils.checkDatasetGeneral(stripped_ds);

    let summ = await utils.checkDatasetSummary(stripped_ds);
    expect(Object.keys(summ.modality_features)).toEqual([""]);

    let loaded = await utils.checkDatasetLoad(stripped_ds);
    expect(loaded.matrix.available()).toEqual(["RNA"]);

    let copy = await utils.checkDatasetSerialize(stripped_ds);
    utils.sameDatasetSummary(summ, await copy.summary());
    utils.sameDatasetLoad(loaded, await copy.load());

    stripped_ds.clear();
})

test("AlabasterAbstractDataset for multimodal datasets", async () => {
    for (const mode of ["sparse", "dense"]) {
        let full_ds = new LocalAlabasterDataset("files/datasets/alabaster/zeisel-brain-" + mode);
        await utils.checkDatasetGeneral(full_ds);

        let summ = await utils.checkDatasetSummary(full_ds);
        expect(Object.keys(summ.modality_features)).toEqual(["gene", "repeat", "ERCC"]);

        full_ds.setOptions({
            rnaExperiment: "gene",
            adtExperiment: "ERCC",
            crisprExperiment: "repeat",
            rnaCountAssay: "counts",
            adtCountAssay: "counts",
            crisprCountAssay: "counts"
        });

        let loaded = await utils.checkDatasetLoad(full_ds);
        expect(loaded.matrix.available()).toEqual(["RNA", "ADT", "CRISPR"]);

        let copy = await utils.checkDatasetSerialize(full_ds);
        utils.sameDatasetSummary(summ, await copy.summary());
        utils.sameDatasetLoad(loaded, await copy.load());

        full_ds.clear();
    }
})

/***********************************************/

class LocalAlabasterResult extends bakana.AbstractAlabasterResult {
    #dir;
    constructor(dir, options={}) {
        super(new LocalDirectoryNavigator(dir), options);
        this.#dir = dir;
    }
}

test("AlabasterAbstractResult behaves with simple results", async () => {
    for (const extra of [ "", "-delayed", "-delayed-external" ]) {
        let res = new LocalAlabasterResult("files/datasets/alabaster/zeisel-brain-sparse-results" + extra);

        let summ = await utils.checkResultSummary(res);
        expect(Object.keys(summ.modality_features)).toEqual(["rna"]);
        expect(summ.reduced_dimension_names).toEqual(["pca", "tsne", "umap"]);
        expect(summ.modality_assay_names["rna"]).toEqual(["filtered", "normalized"]);

        // Checking that we can load the log counts.
        res.setOptions({
            primaryAssay: "normalized",
            reducedDimensionNames: [ "tsne", "umap" ]
        });

        let loaded = await utils.checkResultLoad(res);
        expect(loaded.matrix.available()).toEqual(["rna"]);

        expect(Object.keys(loaded.reduced_dimensions)).toEqual(["tsne", "umap"]);
        for (const [key, val] of Object.entries(loaded.reduced_dimensions)) {
            expect(val.length).toEqual(2);
        }
    }
})

test("AlabasterAbstractResult performs normalization", async () => {
    let res = new LocalAlabasterResult("files/datasets/alabaster/zeisel-brain-sparse-results");
    res.setOptions({
        primaryAssay: "filtered",
        isPrimaryNormalized: false
    });

    let loaded = await utils.checkResultLoad(res);
    expect(loaded.matrix.available()).toEqual(["rna"]);
})

test("AlabasterAbstractResult behaves with multimodal results", async () => {
    for (const extra of [ "", "-delayed", "-delayed-external" ]) {
        let res = new LocalAlabasterResult("files/datasets/alabaster/zeisel-brain-dense-multimodal-results" + extra);

        let summ = await utils.checkResultSummary(res);
        expect(Object.keys(summ.modality_features)).toEqual(["rna", "adt"]);
        expect(summ.reduced_dimension_names).toEqual(["pca", "combined.pca", "tsne", "umap"]);
        expect(summ.modality_assay_names["rna"]).toEqual(["filtered", "normalized"]);
        expect(summ.modality_assay_names["adt"]).toEqual(["filtered", "normalized"]);

        // Checking that we can load the log counts.
        res.setOptions({
            primaryAssay: { rna: "normalized", adt: 0 },
            reducedDimensionNames: [ "tsne", "umap" ]
        });

        let loaded = await utils.checkResultLoad(res);
        expect(loaded.matrix.available()).toEqual(["rna", "adt"]);

        expect(Object.keys(loaded.reduced_dimensions)).toEqual(["tsne", "umap"]);
        for (const [key, val] of Object.entries(loaded.reduced_dimensions)) {
            expect(val.length).toEqual(2);
        }
    }
})
