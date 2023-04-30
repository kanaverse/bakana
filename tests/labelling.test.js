import * as bakana from "../src/index.js";
import * as butils from "../src/steps/utils/general.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";
import * as bioc from "bioconductor";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

let hs_files = { 
    default: new bakana.TenxMatrixMarketDataset(
        "files/datasets/pbmc3k-matrix.mtx.gz",
        "files/datasets/pbmc3k-features.tsv.gz",
        "files/datasets/pbmc3k-barcodes.tsv.gz")
};

test("labelling works correctly for single human reference", async () => {
    let params = utils.baseParams();
    params.cell_labelling.references = [ "ImmGen", "BlueprintEncode" ];
    let state = await bakana.createAnalysis();
    await bakana.runAnalysis(state, hs_files, params);

    let markers = state.marker_detection.fetchResults().RNA;
    let nclust = markers.numberOfGroups();

    expect(state.cell_labelling.fetchNumberOfSharedFeatures().BlueprintEncode).toBeGreaterThan(0);

    {
        let res = state.cell_labelling.computeLabels(markers);
        expect(res.per_reference["BlueprintEncode"].length).toBe(nclust);

        let best = res.per_reference["BlueprintEncode"][0].best;
        expect(typeof best).toBe("string");
        expect(typeof res.per_reference["BlueprintEncode"][0].all[best]).toBe('number');

        expect("ImmGen" in res.per_reference).toBe(false);
        expect("integrated" in res).toBe(false);
    }

    // Forcing us to load the mice.
    params.cell_labelling.automatic = false;
    params.cell_labelling.species = ["10090"];
    params.cell_labelling.gene_id_column = "id";
    {
        await bakana.runAnalysis(state, hs_files, params);
        expect(state.marker_detection.changed).toBe(false);
        expect(state.cell_labelling.changed).toBe(true);

        let res = state.cell_labelling.computeLabels(markers);
        expect(res.per_reference["ImmGen"].length).toEqual(nclust);
        expect(state.cell_labelling.fetchNumberOfSharedFeatures().ImmGen).toBe(0); // well, there's no overlap with mouse IDs, obviously.
        expect("BlueprintEncode" in res.per_reference).toBe(false);
    }

    // Checking that reloading works.
    params.cell_labelling.automatic = true;
    {
        await bakana.runAnalysis(state, hs_files, params);
        expect(state.marker_detection.changed).toBe(false);
        expect(state.cell_labelling.changed).toBe(true);

        let res = state.cell_labelling.computeLabels(markers);
        expect("BlueprintEncode" in res.per_reference).toBe(true);
    }

    params.cell_labelling.species = ["9606"]; 
    {
        await bakana.runAnalysis(state, hs_files, params);
        expect(state.marker_detection.changed).toBe(false);
        expect(state.cell_labelling.changed).toBe(false); // no effect when automatic = true.
    }

    // Release me!
    await bakana.freeAnalysis(state);
})

let mm_files = { 
    default: new bakana.H5adDataset("files/datasets/zeisel-brain.h5ad")
};

test("labelling works correctly for multiple mice references", async () => {
    let params = utils.baseParams();
    params.cell_labelling.references = [ "ImmGen", "MouseRNAseq" ];
    let state = await bakana.createAnalysis();
    await bakana.runAnalysis(state, mm_files, params);

    let markers = state.marker_detection.fetchResults().RNA;
    let nclust = markers.numberOfGroups();

    {
        let res = await state.cell_labelling.computeLabels(markers);
        expect(res.per_reference["ImmGen"].length).toEqual(nclust);
        expect(res.per_reference["MouseRNAseq"].length).toEqual(nclust);
        expect(res.integrated.length).toEqual(nclust);
    }

    // Checking that automatic discovery actually has an effect.
    params.cell_labelling.automatic = false;
    params.cell_labelling.gene_id_type = "ENSEMBL";
    params.cell_labelling.species = ["10090"];
    {
        await bakana.runAnalysis(state, mm_files, params);
        expect(state.cell_labelling.fetchNumberOfSharedFeatures().ImmGen).toBe(0); // well, there's no overlap with Ensembl IDs, obviously.
        expect(state.cell_labelling.fetchNumberOfSharedFeatures().MouseRNAseq).toBe(0); 
    }

    params.cell_labelling.gene_id_type = "SYMBOL";
    {
        await bakana.runAnalysis(state, mm_files, params);
        expect(state.cell_labelling.fetchNumberOfSharedFeatures().ImmGen).toBeGreaterThan(0);
        expect(state.cell_labelling.fetchNumberOfSharedFeatures().MouseRNAseq).toBeGreaterThan(0);
    }

    // Release me!
    await bakana.freeAnalysis(state);
})
