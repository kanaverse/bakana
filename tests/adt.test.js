import * as bakana from "../src/index.js";
import * as utils from "./utils.js";
import * as scran from "scran.js";
import * as valkana from "valkana";
import * as fs from "fs";
import * as combine from "../src/steps/combine_embeddings.js";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

test("runAnalysis works correctly (MatrixMarket)", async () => {
    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    let mtx = "files/datasets/immune_3.0.0-matrix.mtx.gz";
    let feats = "files/datasets/immune_3.0.0-features.tsv.gz";
    let res = await bakana.runAnalysis(state, 
        { 
            default: {
                format: "MatrixMarket",
                mtx: mtx,
                genes: feats,
                annotations: "files/datasets/immune_3.0.0-barcodes.tsv.gz"
            }
        },
        params
    );

    // Checking that the ADTs were split out.
    expect(state.inputs.hasAvailable("RNA")).toBe(true);
    expect(state.inputs.hasAvailable("ADT")).toBe(true);

    // Input reorganization is done correctly. 
    {
        let simple = scran.initializeSparseMatrixFromMatrixMarket(mtx, { layered: false });
        let parsed = bakana.readTable(feats, { compression: "gz" });
        let simple_names = parsed.map(x => x[0]);
        let simple_types = parsed.map(x => x[2]);

        // For RNA.
        {
            let keep = [];
            simple_types.forEach((x, i) => {
                if (x.match(/Gene Expression/i)) {
                    keep.push(i);
                }
            });

            let simple_names_rna = keep.map(x => simple_names[x]);
            let simple_rna = scran.subsetRows(simple.matrix, keep);
            let simple_ids_rna = keep.slice();

            let loaded = state.inputs.fetchCountMatrix({ type: "RNA" });
            let loaded_ids = state.inputs.fetchCountMatrix({ type : "RNA" });
            let loaded_names = state.inputs.fetchGenes({ type: "RNA" }).id;

            expect(simple.matrix.numberOfRows()).toBeGreaterThan(loaded.numberOfRows());
            let sorted = Array.from(loaded_ids).sort((a, b) => a - b);
            expect(sorted).toEqual(keep);

            utils.checkReorganization(simple_rna, simple_ids_rna, simple_names_rna, loaded, loaded_ids, loaded_names, { referenceSubset: true }); 
            simple_rna.free();
        }

        // For ADT.
        {
            let keep = [];
            simple_types.forEach((x, i) => {
                if (x.match(/Antibody Capture/i)) {
                    keep.push(i);
                }
            });

            let simple_names_adt = keep.map(x => simple_names[x]);
            let simple_adt = scran.subsetRows(simple.matrix, keep);
            let simple_ids_adt = keep.slice();

            let loaded = state.inputs.fetchCountMatrix({ type: "ADT" });
            let loaded_ids = state.inputs.fetchCountMatrix({ type : "ADT" });
            let loaded_names = state.inputs.fetchGenes({ type: "ADT" }).id;

            expect(simple.matrix.numberOfRows()).toBeGreaterThan(loaded.numberOfRows());
            let sorted = Array.from(loaded_ids).sort((a, b) => a - b);
            expect(sorted).toEqual(keep);

            utils.checkReorganization(simple_adt, simple_ids_adt, simple_names_adt, loaded, loaded_ids, loaded_names, { referenceSubset: true }); 
            simple_adt.free();
        }

        simple.matrix.free();
    }

    // Checking all the computations.
    {
        // QC.
        let summ = state.adt_quality_control.summary();
        let positive_total = 0;
        summ.data.default.igg_total.forEach(x => { positive_total += (x > 0); });
        expect(positive_total).toBeGreaterThan(0);
        expect(summ.thresholds.default.detected).toBeGreaterThan(0);
        expect(summ.thresholds.default.igg_total).toBeGreaterThan(0);
    }

    {
        // Normalization.
        let norm = state.adt_normalization.summary();
        expect(norm.size_factors.length).toBeGreaterThan(0);
        let positive_total = 0;
        norm.size_factors.forEach(x => { positive_total += (x > 0); });
        expect(positive_total).toBeGreaterThan(0);
    }

    // Saving and loading.
    const path = "TEST_state_adt.h5";
    let collected = await bakana.saveAnalysis(state, path);
    utils.validateState(path);
    expect(collected.collected.length).toBe(3);
    expect(typeof(collected.collected[0])).toBe("string");

    let offsets = utils.mockOffsets(collected.collected);
    let reloaded = await bakana.loadAnalysis(
        path, 
        (offset, size) => offsets[offset]
    );

    let new_params = bakana.retrieveParameters(reloaded);

    {
        // Check that steps unserialize correctly.
        let old_keys = Object.keys(state);
        old_keys.sort();
        let new_keys = Object.keys(reloaded);
        new_keys.sort();
        expect(old_keys).toEqual(new_keys);

        for (const step of old_keys) {
            let qc_deets = reloaded[step].summary();
            let ref = state[step].summary();
            expect(ref).toEqual(ref);
        }

        // Check that we still get some markers.
        let reloaded_markers = reloaded.marker_detection.fetchGroupResults(0, "auc-min-rank", "ADT");
        let ref_markers = state.marker_detection.fetchGroupResults(0, "auc-min-rank", "ADT");
        expect(reloaded_markers).toEqual(ref_markers);

        // Check that the QC steps got dragged out.
        let adt_sums = state.adt_quality_control.fetchSums();
        let positive = 0;
        adt_sums.forEach(x => { positive += (x > 0); });
        expect(positive).toBe(adt_sums.length);
        expect(adt_sums).toEqual(state.adt_quality_control.fetchSums());
    }

    // Checking that the permutation is unchanged on reload, 
    // even when identities are subsetted.
    let old_adt_ids = state.inputs.summary()["genes"]["ADT"]["id"];
    let new_adt_ids = reloaded.inputs.summary()["genes"]["ADT"]["id"];
    expect(old_adt_ids.length).toBeGreaterThan(0);
    expect(old_adt_ids).toEqual(new_adt_ids);

    let old_rna_ids = state.inputs.summary()["genes"]["RNA"]["id"];
    let new_rna_ids = reloaded.inputs.summary()["genes"]["RNA"]["id"];
    expect(old_rna_ids.length).toBeGreaterThan(0);
    expect(old_rna_ids).toEqual(new_rna_ids);

    let old_res = state.feature_selection.summary();
    let new_res = reloaded.feature_selection.summary();
    expect("means" in old_res).toBe(true);
    expect(old_res["means"]).toEqual(new_res["means"]);

    // Release me!
    await bakana.freeAnalysis(state);
    await bakana.freeAnalysis(reloaded);
})

test("runAnalysis works correctly (10X)", async () => {
    let state = await bakana.createAnalysis();
    let params = utils.baseParams();
    let res = await bakana.runAnalysis(state, 
        { 
            default: {
                format: "10X",
                h5: "files/datasets/immune_3.0.0-tenx.h5"
            }
        },
        params
    );

    // Checking that the ADTs were split out.
    expect(state.inputs.hasAvailable("RNA")).toBe(true);
    expect(state.inputs.hasAvailable("ADT")).toBe(true);

    // What happens when one of the modalities has zero weight?
    {
        state.combine_embeddings.compute({ "RNA": 1, "ADT": 0 }, true);
        let pcs = state.combine_embeddings.fetchPCs();
        expect(pcs.pcs.owner !== null).toBe(true);
        expect(pcs.num_pcs).toBe(10);

        const path = "TEST_state_combine-embed.h5";
        let fhandle = scran.createNewHDF5File(path);
        state.combine_embeddings.serialize(fhandle);

        {
            let ncells = state.cell_filtering.summary().retained;
            let npcs_rna = state.pca.summary().var_exp.length;
            let npcs_adt = state.adt_pca.summary().var_exp.length;
            valkana.validateCombineEmbeddingsState(path, ncells, ["RNA"], npcs_rna + npcs_adt, bakana.kanaFormatVersion);
        }

        let reloaded = combine.unserialize(fhandle, {"RNA": state.pca, "ADT": state.adt_pca});
        let repcs = reloaded.fetchPCs();
        expect(repcs.pcs.owner !== null).toBe(true);
        expect(repcs.num_pcs).toBe(10);

        reloaded.free();
    }

    // Release me!
    await bakana.freeAnalysis(state);
})

test("runAnalysis works for ADTs with blocking", async () => {
    let state = await bakana.createAnalysis();
    let params = utils.baseParams();

    // Mocking up a blocking file with pretend batches.
    let exfile = "TEST_adt_block.tsv";
    {
        let previous = "files/datasets/immune_3.0.0-barcodes.tsv.gz";
        let f = fs.readFileSync(previous);
        let buff = f.buffer.slice(f.byteOffset, f.byteOffset + f.byteLength);
        let stuff = bakana.readTable(new Uint8Array(buff));

        let ncells = stuff.length;
        let per_block = Math.ceil(ncells / 3);
        let blocks = new Array(ncells);
        for (var c = 0; c < ncells; c++) {
            blocks[c] = 'A' + String(Math.floor(c / per_block));
        }

        fs.writeFileSync(exfile, blocks.join("\n"));
    }
    params.inputs.sample_factor = "A0";

    let res = await bakana.runAnalysis(state, 
        {
            "combined": {
                format: "MatrixMarket",
                mtx: "files/datasets/immune_3.0.0-matrix.mtx.gz",
                genes: "files/datasets/immune_3.0.0-features.tsv.gz",
                annotations: exfile
            }
        },
        params
    );

    let qcstate = state.adt_quality_control;
    let summ = qcstate.summary();

    // Checking that that there are multiple metrics.
    let positive_total = 0;
    summ.data["A0"].detected.forEach(x => { positive_total += (x > 0); });
    expect(positive_total).toBeGreaterThan(0);

    positive_total = 0;
    summ.data["A2"].igg_total.forEach(x => { positive_total += (x > 0); });
    expect(positive_total).toBeGreaterThan(0);

    // Checking that the thresholds are sensible.
    expect(summ.thresholds["A0"].detected).toBeGreaterThan(0);
    expect(summ.thresholds["A1"].detected).toBeGreaterThan(0);
    expect(summ.thresholds["A1"].igg_total).toBeGreaterThan(0);

    // Freeing everyone.
    await bakana.freeAnalysis(state);
})
