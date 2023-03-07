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

    // Also works with blocking.
    let block = new Int32Array(loaded.matrix.numberOfColumns());
    for (var i = 0; i < block.length; i++) {
        block[i] = i % 2 == 0;
    }

    {
        // Checks kick in.
        expect(() => new bakana.MarkerDetectionStandalone(normed, groups, { block: [] })).toThrow("same length");

        let markers2 = new bakana.MarkerDetectionStandalone(normed, groups, { block: block } );
        expect(markers2.fetchBlockLevels()).toEqual([1, 0]);

        markers2.computeAll();
        let res2 = markers2.fetchResults();
        expect(res2["RNA"].cohen(0)).not.toEqual(res["RNA"].cohen(0));

        markers2.free();
    }

    // Sanitization works correctly.
    {
        let lev = ["B", "A"];
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
    custom.addSelection("BLAH", [ 11,12,13,14,15,16,17 ]);

    let res = custom.fetchResults("WHEEE");
    expect(res["foo"] instanceof scran.ScoreMarkersResults).toBe(true);
    expect(res["bar"] instanceof scran.ScoreMarkersResults).toBe(true);

    let vres = custom.computeVersus("WHEEE", "BLAH");
    expect(vres.results["foo"] instanceof scran.ScoreMarkersResults).toBe(true);
    expect(vres.results["bar"] instanceof scran.ScoreMarkersResults).toBe(true);

    // Also works with blocking.
    let block = new Int32Array(loaded.matrix.numberOfColumns());
    for (var i = 0; i < block.length; i++) {
        block[i] = i % 2 == 0;
    }

    {
        // Checks kick in.
        expect(() => new bakana.CustomSelectionsStandalone(normed, { block: ["A"] })).toThrow("non-negative integer");
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
    let params = bakana.FeatureSetEnrichmentState.defaults();
    params.collections.push("human-GO");
    await enrich.setParameters(params);

    let deets = enrich.fetchCollectionDetails();
    expect("human-GO" in deets).toBe(true);
    let res = enrich.computeEnrichment(markers, 1, "cohen", "min_rank");
    expect(res["human-GO"].num_markers).toBeGreaterThan(0);
    expect(res["human-GO"].counts).not.toEqual(new Int32Array(normed.numberOfRows()));

    let chosen = -1;
    let human_sizes = deets["human-GO"].sizes;
    for (var i = 0; i < human_sizes.length; i++) {
        if (human_sizes[i] > 20) {
            chosen = i;
            break;
        }
    }
    let scores = enrich.computePerCellScores("human-GO", chosen);
    expect(scores.weights.length).toEqual(human_sizes[chosen]);

    // Also works with blocking.
    let block = new Int32Array(loaded.matrix.numberOfColumns());
    for (var i = 0; i < block.length; i++) {
        block[i] = i % 2 == 0;
    }

    {
        let enrich2 = new bakana.FeatureSetEnrichmentStandalone(loaded.features.RNA, { normalized: normed, block: block });
        await enrich2.setParameters(params);
        let scores2 = enrich2.computePerCellScores("human-GO", chosen);
        expect(scores2).not.toEqual(scores);
        enrich2.free();
    }

    // Throws an error correctly if no normalized matrix is supplied.
    {
        let enrich2 = new bakana.FeatureSetEnrichmentStandalone(loaded.features.RNA);
        await enrich2.setParameters(params);
        expect(() => enrich2.computePerCellScores("human-GO", chosen)).toThrow("no normalized matrix");
    }

    markers.free();
    normed.free();
    enrich.free();
})

