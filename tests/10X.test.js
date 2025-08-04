import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";
import * as bioc from "bioconductor";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

test("10X HDF5 readers work in simple cases", async () => {
    const h5path = "files/datasets/pbmc4k-tenx.h5";
    let ds = new bakana.TenxHdf5Dataset(h5path);
    await utils.checkDatasetGeneral(ds);

    let summ = await utils.checkDatasetSummary(ds);
    expect(Object.keys(summ.modality_features)).toEqual([""]);

    let loaded = await utils.checkDatasetLoad(ds);
    expect(Object.keys(loaded.features)).toEqual(["RNA"]);

    let copy = await utils.checkDatasetSerialize(ds);
    utils.sameDatasetSummary(summ, await copy.summary());
    utils.sameDatasetLoad(loaded, await copy.load());

    // Input reorganization is done correctly.
    {
        let mat = loaded.matrix.get("RNA");
        let names = loaded.primary_ids["RNA"];
        let simple = scran.initializeSparseMatrixFromHdf5(h5path, "matrix", { layered: false });
        let simple_names = (new scran.H5File(h5path)).open("matrix").open("features").open("id", { load: true }).values;
        utils.checkMatrixContents(simple, simple_names, mat, names);
        simple.free();
    }

    // Trying what happens if we force it to be another modality.
    ds.setOptions({ featureTypeAdtName: "", featureTypeCrisprName: null, featureTypeRnaName: null });
    let loaded2 = await utils.checkDatasetLoad(ds);
    expect(Object.keys(loaded2.features)).toEqual(["ADT"]);

    ds.clear();
})

test("10X HDF5 readers work for multiple modalities", async () => {
    let ds = new bakana.TenxHdf5Dataset("files/datasets/crispr_6.0.0-tenx.h5");
    await utils.checkDatasetGeneral(ds);

    let summ = await utils.checkDatasetSummary(ds);
    expect(Object.keys(summ.modality_features)).toEqual(["Gene Expression", "CRISPR Guide Capture"]);

    let loaded = await utils.checkDatasetLoad(ds);
    expect(Object.keys(loaded.features)).toEqual(["RNA", "CRISPR"]);

    let copy = await utils.checkDatasetSerialize(ds);
    utils.sameDatasetSummary(summ, await copy.summary());
    utils.sameDatasetLoad(loaded, await copy.load());

    ds.clear();
})
