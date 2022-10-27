import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";

beforeAll(async () => await bakana.initialize({ localFile: true }));
afterAll(async () => await bakana.terminate());

test("annotation preflight works correctly (one file)", async () => {
    let res = await bakana.validateAnnotations(
        {
            default: new bakana.TenxMatrixMarketDataset(
                    "files/datasets/pbmc3k-matrix.mtx.gz",
                    "files/datasets/pbmc3k-features.tsv.gz",
                    "files/datasets/pbmc3k-barcodes.tsv.gz"
                )
        }
    );

    expect(res.features.RNA.common).toBeGreaterThan(0);
    console.log(res.features.RNA.fields);
    expect(typeof res.features.RNA.fields.default).toBe("string");

    let default_anno = res.annotations.default;
    let default_keys = Object.keys(default_anno);
    expect(default_keys.length).toBeGreaterThan(0);
    
    let first = default_anno[default_keys[0]];
    expect(first.values.length).toBeGreaterThan(0);
    expect(first.truncated).toBe(true); // as these are just the cell names, so they'll be cut off.

    // still works without any annotations.
    let res2 = await bakana.validateAnnotations(
        {
            default: new bakana.TenxMatrixMarketDataset("files/datasets/pbmc3k-matrix.mtx.gz", null, null)
        }
    );

    expect(res2.annotations.default).toEqual({});
    expect(res2.features.RNA.common).toBeGreaterThan(0);
})

test("annotation preflight works correctly (two files)", async () => {
    let res = await bakana.validateAnnotations(
        {
            "3k": new bakana.TenxMatrixMarketDataset(
                    "files/datasets/pbmc3k-matrix.mtx.gz",
                    "files/datasets/pbmc3k-features.tsv.gz",
                    "files/datasets/pbmc3k-barcodes.tsv.gz"
                ),
            "4k": new bakana.TenxHdf5Dataset("files/datasets/pbmc4k-tenx.h5")
        }
    );

    expect(res.features.RNA.common).toBeGreaterThan(30000);
    expect(res.features.RNA.fields["3k"]).toBe("id");
    expect(res.features.RNA.fields["4k"]).toBe("id");

    // Checking the annotations while we're here.
    expect(Object.keys(res.annotations["3k"]).length).toBeGreaterThan(0);
    expect(res.annotations["4k"]).toEqual({});
})

test("annotation preflight works correctly for H5ADs", async () => {
    let res = await bakana.validateAnnotations(
        {
            "brain": new bakana.H5adDataset("files/datasets/zeisel-brain.h5ad")
        }
    );

    let brain_anno = res.annotations.brain;
    let brain_keys = Object.keys(brain_anno);
    expect(brain_keys.length).toBeGreaterThan(0);
    
    let cells = brain_anno["cell_id"];
    expect(cells.values.length).toBeGreaterThan(0);
    expect(typeof cells.values[0]).toBe("string");
    expect(cells.truncated).toBe(true); 

    let sex = brain_anno["sex"];
    expect(sex.min).toBe(1);
    expect(sex.max).toBe(3);
})

test("annotation preflight works correctly for SummarizedExperiments", async () => {
    let res = await bakana.validateAnnotations(
        {
            "brain": new bakana.SummarizedExperimentDataset("files/datasets/zeisel-brain.rds")
        }
    );

    expect(res.features.RNA.common).toBeGreaterThan(0);

    let brain_anno = res.annotations.brain;
    let brain_keys = Object.keys(brain_anno);
    expect(brain_keys.length).toBeGreaterThan(0);

    console.log(brain_anno);
    let cells = brain_anno["cell_id"];
    console.log(cells);
    expect(cells.values.length).toBeGreaterThan(0);
    expect(typeof cells.values[0]).toBe("string");
    expect(cells.truncated).toBe(true); 

    let sex = brain_anno["sex"];
    expect(sex.min).toBe(1);
    expect(sex.max).toBe(3);
})

test("annotation preflight works correctly for SingleCellExperiments with altExps", async () => {
    let res = await bakana.validateAnnotations(
        {
            "immune": new bakana.SummarizedExperimentDataset("files/datasets/immune_3.0.0-tenx.rds")
        }
    );

    expect(res.features.RNA.common).toBeGreaterThan(0);
    expect(res.features.ADT.common).toBeGreaterThan(0);

    let immune_anno = res.annotations.immune;
    let immune_keys = Object.keys(immune_anno);

    let samples = immune_anno['Sample'];
    expect(samples.values.length).toBeGreaterThan(0);
    expect(typeof samples.values[0]).toBe("string");
})

test("annotation preflight fails correctly (two files, wrong species)", async () => {
    let res;
    let err;
    try {
        res = await bakana.validateAnnotations(
            {
                "brain": new bakana.H5adDataset("files/datasets/zeisel-brain.h5ad"),
                "4k": new bakana.TenxHdf5Dataset("files/datasets/pbmc4k-tenx.h5")
            }
        );
    } catch (e) {
        err = e.toString();
    }

    expect(res).toBeUndefined();
    expect(err).toMatch("common feature type");
})

test("annotation preflight fails correctly (two files, no genes)", async () => {
    let res;
    let err;
    try {
        res = await bakana.validateAnnotations(
            {
                "brain": new bakana.H5adDataset("files/datasets/zeisel-brain.h5ad"),
                "3k": new bakana.TenxMatrixMarketDataset("files/datasets/pbmc3k-matrix.mtx.gz")
            }
        );
    } catch (e) {
        err = e.toString();
    }

    expect(res).toBeUndefined();
    expect(err).toMatch("gene annotations");
})

test("annotation preflight works correctly (ADTs)", async () => {
    let res = await bakana.validateAnnotations(
        {
            default: new bakana.TenxMatrixMarketDataset(
                    "files/datasets/immune_3.0.0-matrix.mtx.gz",
                    "files/datasets/immune_3.0.0-features.tsv.gz",
                    "files/datasets/immune_3.0.0-barcodes.tsv.gz"
                ) 
        }
    );

    expect(Object.keys(res.annotations.default).length).toBeGreaterThan(0);
    expect(res.features.RNA.common).toBeGreaterThan(0);
    expect(res.features.ADT.common).toBeGreaterThan(0);

    // still works with multiple hits.
    let res2 = await bakana.validateAnnotations(
        {
            mtx: new bakana.TenxMatrixMarketDataset(
                    "files/datasets/immune_3.0.0-matrix.mtx.gz",
                    "files/datasets/immune_3.0.0-features.tsv.gz",
                    "files/datasets/immune_3.0.0-barcodes.tsv.gz"
                ),
            tenx: new bakana.TenxHdf5Dataset("files/datasets/immune_3.0.0-tenx.h5")
        }
    );

    expect(Object.keys(res.annotations.default).length).toBeGreaterThan(0);
    expect(res.features.RNA.common).toBeGreaterThan(0);
    expect(res.features.ADT.common).toBeGreaterThan(0);
})


