import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as bioc from "bioconductor";
import * as utils from "./utils.js";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

let h5path = "files/datasets/pbmc4k-tenx.h5";
let files = {
    default: new bakana.TenxHdf5Dataset(h5path)
};

test("Standalone marker detection works correctly", async () => {
    let loaded = files.default.load({ cache: true });
    let normed = new scran.MultiMatrix;
    for (const m of loaded.matrix.available()) {
        normed.add(m, scran.logNormCounts(loaded.matrix.get(m)));
    }

    let groups = new Int32Array(loaded.matrix.numberOfColumns());
    let half = groups.length / 2
    groups.fill(0, 0, half);
    groups.fill(1, half, groups.length);

    // Checks kick in.
    expect(() => new bakana.MarkerDetectionStandalone(normed, [])).toThrow("same number of columns");

    let markers = new bakana.MarkerDetectionStandalone(normed, groups);
    expect(markers.fetchGroupLevels()).toEqual([0, 1]);
    expect(markers._peekGroups().array()).toEqual(groups);

    markers.computeAll();
    let res = markers.fetchResults();
    expect(res["RNA"] instanceof scran.ScoreMarkersResults).toBe(true);

    let versus = markers.computeVersus(1, 0);
    expect(versus.results["RNA"] instanceof scran.ScoreMarkersResults).toBe(true);

    // Setting options works.
    {
        let markers2 = new bakana.MarkerDetectionStandalone(normed, groups);

        let defaults = bakana.MarkerDetectionState.defaults();
        defaults.compute_auc = false;
        markers2.setParameters(defaults);
        markers2.computeAll();
        let res2 = markers2.fetchResults();

        let rounded = res["RNA"].cohen(0).map(x => Math.round(x * 100));
        let rounded2 = res2["RNA"].cohen(0).map(x => Math.round(x * 100));
        expect(rounded2).toEqual(rounded); // avoid discrepancies due to numerical precision.
        expect(() => res2["RNA"].auc(0)).toThrow("no AUCs");

        markers2.free();
    }

    // Also works with blocking.
    let block = new Int32Array(loaded.matrix.numberOfColumns());
    for (var i = 0; i < block.length; i++) {
        block[i] = i % 2 == 0;
    }

    {
        // Checks kick in.
        expect(() => new bakana.MarkerDetectionStandalone(normed, groups, { block: [] })).toThrow("same length");

        let markers2 = new bakana.MarkerDetectionStandalone(normed, groups, { block: block } );
        expect(markers2.fetchBlockLevels()).toEqual([0, 1]);

        markers2.computeAll();
        let res2 = markers2.fetchResults();
        expect(res2["RNA"].cohen(0)).not.toEqual(res["RNA"].cohen(0));

        markers2.free();
    }

    // Sanitization works correctly.
    {
        let lev = ["A", "B"];
        let groups2 = Array.from(groups).map(i => lev[i]);
        groups2[0] = null;

        let partial = new bakana.MarkerDetectionStandalone(normed, groups2);
        expect(partial.fetchGroupLevels()).toEqual(lev);
        expect(partial._peekMatrices().get("RNA").row(0)).toEqual(normed.get("RNA").row(0).slice(1));
        expect(partial._peekGroups().array()).toEqual(groups.slice(1));
        expect(partial._peekBlock()).toBeNull();
        partial.free();

        let blev = ["x", "y"];
        let block2 = Array.from(block).map(i => blev[i]);
        let last = block2.length - 1;
        block2[last] = null;

        let bpartial = new bakana.MarkerDetectionStandalone(normed, groups2, { block: block2 });
        expect(bpartial._peekMatrices().get("RNA").row(0)).toEqual(normed.get("RNA").row(0).slice(1, last));
        expect(bpartial._peekGroups().array()).toEqual(groups.slice(1, last));
        let extracted = bpartial.fetchBlockLevels();
        expect(Array.from(bpartial._peekBlock().array()).map(i => extracted[i])).toEqual(Array.from(block.slice(1, last)).map(i => blev[i]));
        bpartial.free();
    }

    normed.free();
    markers.free();
})

test("Standalone custom selections work correctly", async () => {
    let loaded = files.default.load({ cache: true });
    let normed = new scran.MultiMatrix;
    for (const m of [ "foo", "bar" ]) { // mimicking multiple modalities.
        normed.add(m, scran.logNormCounts(loaded.matrix.get("RNA")));
    }

    let custom = new bakana.CustomSelectionsStandalone(normed);
    custom.addSelection("WHEEE", [ 1,2,3,4,5,6,7,8,9,10 ]);
    expect(custom.fetchSelectionIndices("WHEEE")).toEqual([ 1,2,3,4,5,6,7,8,9,10 ]);
    custom.addSelection("BLAH", [ 11,12,13,14,15,16,17 ]);
    expect(custom.fetchSelections()["BLAH"]).toEqual([ 11,12,13,14,15,16,17 ]);

    let res = custom.fetchResults("WHEEE");
    expect(res["foo"] instanceof scran.ScoreMarkersResults).toBe(true);
    expect(res["bar"] instanceof scran.ScoreMarkersResults).toBe(true);

    let vres = custom.computeVersus("WHEEE", "BLAH");
    expect(vres.results["foo"] instanceof scran.ScoreMarkersResults).toBe(true);
    expect(vres.results["bar"] instanceof scran.ScoreMarkersResults).toBe(true);

    // Setting options works.
    {
        let custom2 = new bakana.CustomSelectionsStandalone(normed);

        let defaults = bakana.CustomSelectionsState.defaults();
        defaults.compute_auc = false;
        custom2.setParameters(defaults);
        custom2.addSelection("WHEEE", [ 1,2,3,4,5,6,7,8,9,10 ]);
        let res2 = custom2.fetchResults("WHEEE");

        let rounded = res["foo"].cohen(0).map(x => Math.round(x * 100));
        let rounded2 = res2["foo"].cohen(0).map(x => Math.round(x * 100));
        expect(rounded2).toEqual(rounded); // avoid discrepancies due to numerical precision.
        expect(() => res2["foo"].auc(0)).toThrow("no AUCs");

        custom2.free();
    }

    // Also works with blocking.
    let block = new Int32Array(loaded.matrix.numberOfColumns());
    for (var i = 0; i < block.length; i++) {
        block[i] = i % 2 == 0;
    }

    {
        // Checks kick in.
        expect(() => new bakana.CustomSelectionsStandalone(normed, { block: [] })).toThrow("same length");

        let custom2 = new bakana.CustomSelectionsStandalone(normed, { block: block } );
        custom2.addSelection("WHEEE", [ 1,2,3,4,5,6,7,8,9,10 ]);

        let res2 = custom2.fetchResults("WHEEE");
        expect(res2["foo"].cohen(0)).not.toEqual(res["foo"].cohen(0));

        custom2.free();
    }

    // Sanitization works correctly.
    {
        let blev = ["x", "y"];
        let block2 = Array.from(block).map(i => blev[i]);
        block2[0] = null;

        let bpartial = new bakana.CustomSelectionsStandalone(normed, { block: block2 });
        expect(bpartial._peekMatrices().get("foo").row(0)).toEqual(normed.get("foo").row(0).slice(1));
        let extracted = bpartial.fetchBlockLevels();
        expect(Array.from(bpartial._peekBlock().array()).map(i => extracted[i])).toEqual(Array.from(block.slice(1)).map(i => blev[i]));

        bpartial.addSelection("odds", [ 1, 3, 5, 7, 9 ]);
        expect(bpartial.fetchSelectionIndices("odds")).toEqual([1, 3, 5, 7, 9]);
        bpartial.addSelection("evens", [ 0, 2, 4, 6, 8 ]); 
        expect(bpartial.fetchSelectionIndices("evens")).toEqual([2, 4, 6, 8]); // losing the first column as it's got a missing block.
        expect(bpartial.fetchSelections()["evens"]).toEqual([ 2, 4, 6, 8 ]); // same result.

        bpartial.free();
    }

    normed.free();
    custom.free();
})

test("Standalone feature set enrichment works correctly", async () => {
    let loaded = files.default.load({ cache: true });
    let normed = scran.logNormCounts(loaded.matrix.get("RNA"));

    let groups = new Int32Array(loaded.matrix.numberOfColumns());
    let half = groups.length / 2
    groups.fill(0, 0, half);
    groups.fill(1, half, groups.length);
    let markers = scran.scoreMarkers(normed, groups);

    // Checks kick in.
    {
        let subnorm = scran.subsetRows(normed, [1,2,3]);
        expect(() => new bakana.FeatureSetEnrichmentStandalone(loaded.features.RNA, { normalized: subnorm })).toThrow("number of rows");
        subnorm.free();
    }

    let enrich = new bakana.FeatureSetEnrichmentStandalone(loaded.features.RNA, { normalized: normed });
    await enrich.ready();

    let coldeets = enrich.fetchCollectionDetails();
    let setdeets = enrich.fetchSetDetails();
    expect(Array.from(new Set(coldeets.species))).toEqual(["9606"]);

    let res = enrich.computeEnrichment(markers, 1, "cohen", "min_rank");
    expect(res.num_markers).toBeGreaterThan(0);
    expect(res.counts.length).toBeGreaterThan(0);
    expect(res.counts.every(x => x > 0)).toBe(true);
    expect(res.pvalues.some(Number.isNaN)).toBe(false);

    let chosen = -1;
    for (var i = 0; i < setdeets.sizes.length; i++) {
        if (setdeets.sizes[i] > 20) {
            chosen = i;
            break;
        }
    }
    let scores = enrich.computePerCellScores(chosen);
    expect(scores.scores.some(Number.isNaN)).toBe(false);
    expect(scores.scores.length).toEqual(normed.numberOfColumns());
    expect(scores.weights.length).toEqual(setdeets.sizes[chosen]);

    // no-ops if the parameters haven't changed.
    {
        let params = bakana.FeatureSetEnrichmentState.defaults();
        enrich.setParameters(params);
        await enrich.ready(); 
    }

    // Also works with blocking.
    let block = new Int32Array(loaded.matrix.numberOfColumns());
    for (var i = 0; i < block.length; i++) {
        block[i] = i % 2 == 0;
    }

    {
        let enrich2 = new bakana.FeatureSetEnrichmentStandalone(loaded.features.RNA, { normalized: normed, block: block });
        await enrich2.ready();

        let scores2 = enrich2.computePerCellScores(chosen);
        expect(scores2.scores.length).toEqual(normed.numberOfColumns());
        expect(scores2).not.toEqual(scores);
        enrich2.free();
    }

    // Throws an error correctly if no normalized matrix is supplied.
    {
        let enrich2 = new bakana.FeatureSetEnrichmentStandalone(loaded.features.RNA);
        await enrich2.ready();
        expect(() => enrich2.computePerCellScores("human-GO", chosen)).toThrow("no normalized matrix");
    }

    // Sanitization works correctly.
    {
        let blev = ["x", "y"];
        let block2 = Array.from(block).map(i => blev[i]);
        block2[0] = null;

        let bpartial = new bakana.FeatureSetEnrichmentStandalone(loaded.features.RNA, { normalized: normed, block: block2 });
        expect(bpartial._peekMatrices().row(0)).toEqual(normed.row(0).slice(1));
        expect(bpartial._peekBlock().array()).toEqual(block.slice(1));

        await bpartial.ready();
        let scores = bpartial.computePerCellScores(chosen);
        expect(Number.isNaN(scores.scores[0])).toBe(true);
        expect(scores.scores.slice(1).every(Number.isNaN)).toBe(false);

        bpartial.free();
    }

    markers.free();
    normed.free();
    enrich.free();
})

test("Standalone cell labelling works correctly", async () => {
    let loaded = files.default.load({ cache: true });
    let labeller = new bakana.CellLabellingStandalone(loaded.features.RNA);

    let params = bakana.CellLabellingStandalone.defaults();
    params.references = [ "BlueprintEncode", "MouseRNAseq" ];
    labeller.setParameters(params);
    await labeller.ready();

    let sub = scran.subsetColumns(loaded.matrix.get("RNA"), [0,1,2,3,4,5,6,7,8,9,10]);
    let out = labeller.computeLabels(sub);
    expect(out.per_reference.BlueprintEncode.length).toEqual(sub.numberOfColumns());

    // no-ops if the parameters haven't changed.
    labeller.setParameters(params);
    await labeller.ready(); 

    // Checks kick in.
    {
        let sub2 = scran.subsetRows(loaded.matrix.get("RNA"), [0,1]);
        expect(() => labeller.computeLabels(sub2)).toThrow("unexpected number of genes");
        sub2.free();
    }

    // Sanitization works correctly.
    {
        let grouping = [];
        for (var i = 0; i < loaded.matrix.numberOfColumns(); i++) {
            grouping.push(Math.random() > 0.5 ? "A" : "B");
        }
        grouping[0] = null;
        grouping[1] = null;
        grouping[2] = null;

        let res = labeller.computeLabels(loaded.matrix.get("RNA"), { group: grouping });
        expect(res.groups.slice().sort()).toEqual(["A", "B"]);
        expect(res.per_reference.BlueprintEncode.length).toEqual(2);

        // Same results with manual removal of the first few cells.
        let keep = [];
        for (var i = 3; i < grouping.length; i++) {
            keep.push(i);
        }
        let sub = scran.subsetColumns(loaded.matrix.get("RNA"), keep);
        let subgroup = bioc.SLICE(grouping, keep);
        let subres = labeller.computeLabels(sub, { group: subgroup });
        expect(subres).toEqual(res);
    }

    sub.free();
    labeller.free();
})
