import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";
import * as bioc from "bioconductor";
import * as fs from "fs";
import * as os from "os";
import * as path from "path";
import JSZip from "jszip";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

function compressDirectoryInternal(dir, handle) {
    const listing = fs.readdirSync(dir);
    for (const f of listing) {
        const full = path.join(dir, f);
        let is_dir = fs.statSync(full).isDirectory();
        if (is_dir) {
            var subhandle = handle.folder(f);
            compressDirectoryInternal(full, subhandle);
        } else {
            handle.file(f, fs.readFileSync(full)); 
        }
    }
}

function compressDirectory(dir, prefix = ".") {
    let zip = new JSZip;
    let handle = zip;
    if (prefix !== ".") {
        handle = handle.folder(prefix);
    }
    compressDirectoryInternal(dir, handle);
    return zip.generateAsync({ type: "uint8array", compression: "DEFLATE" })
}

test("searchZippedAlabaster works correctly", async () => {
    const prefix = path.join(os.tmpdir(), "tmp");

    // Top-level OBJECT.
    {
        let dir = fs.mkdtempSync(prefix);
        fs.writeFileSync(path.join(dir, "OBJECT"), JSON.stringify({ type: "summarized_experiment", summarized_experiment: { version: "1.0", dimensions: [10, 20] }}));
        fs.mkdirSync(path.join(dir, "foo"));
        fs.writeFileSync(path.join(dir, "foo", "OBJECT"), JSON.stringify({ type: "summarized_experiment", summarized_experiment: { version: "1.0", dimensions: [20, 40] }}));
        fs.mkdirSync(path.join(dir, "bar"));
        fs.writeFileSync(path.join(dir, "bar", "OBJECT"), JSON.stringify({ type: "summarized_experiment", summarized_experiment: { version: "1.0", dimensions: [30, 60] }}));
        let payload = await compressDirectory(dir);
        let handle = await JSZip.loadAsync(payload);
        let found = await bakana.searchZippedAlabaster(handle);
        expect(found.size).toEqual(1);
        expect(found.get(".")).toEqual([10,20]);
    }

    // No top-level OBJECT.
    {
        let dir = fs.mkdtempSync(prefix);
        fs.mkdirSync(path.join(dir, "foo"));
        fs.writeFileSync(path.join(dir, "foo", "OBJECT"), JSON.stringify({ type: "summarized_experiment", summarized_experiment: { version: "1.0", dimensions: [20, 40] }}));
        fs.mkdirSync(path.join(dir, "bar"));
        fs.writeFileSync(path.join(dir, "bar", "OBJECT"), JSON.stringify({ type: "summarized_experiment", summarized_experiment: { version: "1.0", dimensions: [30, 60] }}));
        let payload = await compressDirectory(dir);
        let handle = await JSZip.loadAsync(payload);
        let found = await bakana.searchZippedAlabaster(handle);
        expect(found.size).toEqual(2);
        expect(found.get("foo")).toEqual([20,40]);
        expect(found.get("bar")).toEqual([30,60]);
    }

    // Nested OBJECT.
    {
        let dir = fs.mkdtempSync(prefix);
        fs.mkdirSync(path.join(dir, "foo"));
        fs.writeFileSync(path.join(dir, "foo", "OBJECT"), JSON.stringify({ type: "summarized_experiment", summarized_experiment: { version: "1.0", dimensions: [20, 40] }}));
        fs.mkdirSync(path.join(dir, "foo", "bar"));
        fs.writeFileSync(path.join(dir, "foo", "bar", "OBJECT"), JSON.stringify({ type: "summarized_experiment", summarized_experiment: { version: "1.0", dimensions: [30, 60] }}));
        let payload = await compressDirectory(dir);
        let handle = await JSZip.loadAsync(payload);
        let found = await bakana.searchZippedAlabaster(handle);
        expect(found.size).toEqual(1);
        expect(found.get("foo")).toEqual([20,40]);
    }
})

test("ZippedAlabasterDataset works correctly", async () => {
    for (const loc of [ ".", "foo/bar" ]) {
        let contents = await compressDirectory("files/datasets/alabaster/zeisel-brain-stripped", loc);
        let stripped_ds = new bakana.ZippedAlabasterDataset(loc, new bakana.SimpleFile(contents, { name: "my.zip" }));

        let summ = await stripped_ds.summary({ cache: true });
        expect(Object.keys(summ.modality_features)).toEqual([""]);
        expect(summ.modality_features[""] instanceof bioc.DataFrame).toBe(true);
        expect(summ.modality_features[""].numberOfRows()).toBeGreaterThan(0);
        expect(summ.cells instanceof bioc.DataFrame).toBe(true);
        expect(summ.cells.numberOfRows()).toBeGreaterThan(0);

        let abbrev = stripped_ds.abbreviate();
        expect(abbrev.files[0].type).toEqual("zip");

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
    }
})

test("ZippedAlabasterResult works correctly", async () => {
    for (const loc of [ ".", "foo/bar" ]) {
        let contents = await compressDirectory("files/datasets/alabaster/zeisel-brain-sparse-results-delayed", loc);
        let res = new bakana.ZippedAlabasterResult(loc, new bakana.SimpleFile(contents, { name: "my.zip" }));

        let summ = await res.summary({ cache: true });
        expect(Object.keys(summ.modality_features)).toEqual(["rna"]);
        const mod = summ.modality_features["rna"];
        expect(mod instanceof bioc.DataFrame).toBe(true);
        expect(mod.numberOfRows()).toBeGreaterThan(0);

        expect(summ.cells instanceof bioc.DataFrame).toBe(true);
        expect(summ.cells.numberOfRows()).toBeGreaterThan(0);
        expect(summ.reduced_dimension_names).toEqual(["pca", "tsne", "umap"]);

        expect(summ.modality_assay_names["rna"]).toEqual(["filtered", "normalized"]);

        // Checking that we can load the log counts.
        res.setOptions({
            primaryAssay: "normalized",
            reducedDimensionNames: [ "tsne", "umap" ]
        });

        let loaded = await res.load({ cache: true });
        expect(loaded.matrix.available()).toEqual(["rna"]);
        expect(loaded.matrix.numberOfColumns()).toEqual(loaded.cells.numberOfRows());
        expect(loaded.matrix.get("rna").numberOfRows()).toEqual(loaded.features["rna"].numberOfRows());

        expect(Object.keys(loaded.reduced_dimensions)).toEqual(["tsne", "umap"]);
        for (const [key, val] of Object.entries(loaded.reduced_dimensions)) {
            expect(val.length).toEqual(2);
            expect(val[0].length).toEqual(loaded.matrix.numberOfColumns());
        }
    }
})
