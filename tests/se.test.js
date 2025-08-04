import * as bakana from "../src/index.js";
import * as inputs from "../src/steps/inputs.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";
import * as bioc from "bioconductor";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

test("RDS dataset readers works for SCEs", async () => {
    const fpath = "files/datasets/zeisel-brain.rds";
    let ds = new bakana.SummarizedExperimentDataset(fpath);
    await utils.checkDatasetGeneral(ds);

    let summ = await utils.checkDatasetSummary(ds);
    expect(Object.keys(summ.modality_features)).toEqual(["endogenous", "ERCC", "repeat"]);

    ds.setOptions({ rnaExperiment: "endogenous" });
    let loaded = await utils.checkDatasetLoad(ds);
    expect(loaded.matrix.available()).toEqual(["RNA"]);

    let copy = await utils.checkDatasetSerialize(ds);
    utils.sameDatasetSummary(summ, await copy.summary());
    utils.sameDatasetLoad(loaded, await copy.load());

    {
        let rdshandle = scran.readRds(fpath);
        let rhandle = rdshandle.value();
        expect(rhandle.className()).toBe("SingleCellExperiment");

        let rrhandle = rhandle.attribute("rowRanges");
        expect(rrhandle.className()).toMatch("GRangesList");
        let nhandle = rrhandle.attribute("partitioning").attribute("NAMES"); // technically a memory leak, but finalizers should handle it.
        let simple_names = nhandle.values();

        let ahandle = rhandle.attribute("assays").attribute("data").attribute("listData").load(0); // again, technically memory leaks here.
        let simple = scran.initializeSparseMatrixFromRds(ahandle, { layered: false });

        ahandle.free();
        nhandle.free();
        rrhandle.free();
        rhandle.free();
        rdshandle.free();

        utils.checkMatrixContents(simple, simple_names, loaded.matrix.get("RNA"), loaded.primary_ids["RNA"]);
        simple.free();
    }

    ds.clear();
})

test("RDS dataset readers work for a base SummarizedExperiment", async () => {
    const fpath = "files/datasets/paul-hsc.rds";
    let ds = new bakana.SummarizedExperimentDataset(fpath);
    await utils.checkDatasetGeneral(ds);

    let summ = await utils.checkDatasetSummary(ds);
    expect(Object.keys(summ.modality_features)).toEqual([""]);

    ds.setOptions({ rnaExperiment: "" });
    let loaded = await utils.checkDatasetLoad(ds);
    expect(loaded.matrix.available()).toEqual(["RNA"]);

    // Check we're dealing with a pure SE.
    {
        let rdshandle = scran.readRds(fpath);
        let rhandle = rdshandle.value();
        expect(rhandle.className()).toBe("SummarizedExperiment");
        rhandle.free();
        rdshandle.free();
    }

    ds.clear();
})

test("RDS dataset readers work for GRanges SE", async () => {
    const fpath = "files/datasets/zeisel-brain-dense.rds";
    let ds = new bakana.SummarizedExperimentDataset(fpath);
    await utils.checkDatasetGeneral(ds);

    let summ = await utils.checkDatasetSummary(ds);
    expect(Object.keys(summ.modality_features)).toEqual(["endogenous", "ERCC", "repeat"]);

    ds.setOptions({ rnaExperiment: "endogenous" });
    let loaded = await utils.checkDatasetLoad(ds);
    expect(loaded.matrix.available()).toEqual(["RNA"]);

    // Check we're dealing with a GRanges-containing SCE.
    {
        let rdshandle = scran.readRds(fpath);
        let rhandle = rdshandle.value();
        expect(rhandle.className()).toBe("SingleCellExperiment");

        let rrhandle = rhandle.attribute("rowRanges");
        expect(rrhandle.className()).toBe("GRanges");

        rrhandle.free();
        rhandle.free();
        rdshandle.free();
    }

    ds.clear();
})

test("RDS dataset readers work for a SingleCellExperiment with altExps", async () => {
    const fpath = "files/datasets/immune_3.0.0-tenx.rds";
    let ds = new bakana.SummarizedExperimentDataset(fpath);
    await utils.checkDatasetGeneral(ds);

    let summ = await utils.checkDatasetSummary(ds);
    expect(Object.keys(summ.modality_features)).toEqual(["Gene Expression", "Antibody Capture"]);

    let loaded = await utils.checkDatasetLoad(ds);
    expect(loaded.matrix.available()).toEqual(["RNA", "ADT"]);
    let loaded_names = loaded.features["RNA"].rowNames();
    expect(loaded_names.some(x => x.match(/^ENSG/))).toBe(true);
    let loaded_names_adt = loaded.features["ADT"].rowNames();
    expect(loaded_names_adt.some(x => x.match(/^ENSG/))).toBe(false);

    // Works with integer codes for the experiment and names for the assays.
    ds.setOptions({ adtExperiment: 0, rnaCountAssay: "counts" });
    let reloaded = await utils.checkDatasetLoad(ds);
    expect(reloaded.matrix.available()).toEqual(["RNA", "ADT"]);

    ds.clear();
})

test("RDS result readers work", async () => {
    let res = new bakana.SummarizedExperimentResult("files/datasets/zeisel-brain-with-results.rds");

    let details = await utils.checkResultSummary(res);
    expect(Object.keys(details.modality_features)).toEqual(["endogenous", "ERCC", "repeat"]);
    expect(details.modality_assay_names["endogenous"]).toEqual(["counts", "logcounts"]);
    expect(details.reduced_dimension_names).toEqual(["PCA", "TSNE", "UMAP"]);

    // Already normalized.
    {
        res.setOptions({ 
            primaryAssay: { "endogenous": "logcounts" },
            isPrimaryNormalized: true,
            reducedDimensionNames: ["TSNE", "UMAP"]
        });

        let payload = await utils.checkResultLoad(res);
        expect(payload.matrix.available()).toEqual(["endogenous"]);
        expect(Object.keys(payload.reduced_dimensions)).toEqual(["TSNE", "UMAP"]);
        expect(payload.reduced_dimensions["TSNE"].length).toEqual(2);
        expect(payload.reduced_dimensions["UMAP"].length).toEqual(2);

        let col0 = payload.matrix.get("endogenous").column(0);
        expect(utils.hasNonInteger(col0)).toBe(true);
    }

    // Manual normalization.
    {
        res.setOptions({
            primaryAssay: "counts",
            isPrimaryNormalized: { "ERCC": false }
        });

        let payload = await utils.checkResultLoad(res);
        expect(payload.matrix.available()).toEqual(["endogenous", "ERCC", "repeat"]);

        // Checking that naked gene counts are loaded, as a control.
        {
            let col0 = payload.matrix.get("endogenous").column(0);
            expect(utils.hasNonInteger(col0)).toBe(false);
        }

        // Checking that the ERCC count matrix is correctly normalized.
        {
            let col0 = payload.matrix.get("ERCC").column(0);
            expect(utils.hasNonInteger(col0)).toBe(true);
        }
    }

    res.clear();
})

