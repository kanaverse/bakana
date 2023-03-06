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
    expect(() => new bakana.MarkerDetectionStandalone(normed, ["A"])).toThrow("non-negative integer");
    expect(() => new bakana.MarkerDetectionStandalone(normed, [])).toThrow("same number of columns");

    let markers = new bakana.MarkerDetectionStandalone(normed, groups);
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
        expect(() => new bakana.MarkerDetectionStandalone(normed, groups, { block: ["A"] })).toThrow("non-negative integer");
        expect(() => new bakana.MarkerDetectionStandalone(normed, groups, { block: [] })).toThrow("same length");

        let markers2 = new bakana.MarkerDetectionStandalone(normed, groups, { block: block } );
        markers2.computeAll();
        let res2 = markers2.fetchResults();
        expect(res2["RNA"].cohen(0)).not.toEqual(res["RNA"].cohen(0));

        markers2.free();
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
        expect(() => new bakana.FeatureSetEnrichmentStandalone(bioc.SLICE(loaded.features.RNA, [1,2,3]), markers)).toThrow("number of rows");
        let subnorm = scran.subsetRows(normed, [1,2,3]);
        expect(() => new bakana.FeatureSetEnrichmentStandalone(loaded.features.RNA, markers, { normalized: subnorm })).toThrow("number of rows");
        subnorm.free();
    }

    let enrich = new bakana.FeatureSetEnrichmentStandalone(loaded.features.RNA, markers, { normalized: normed });
    let params = bakana.FeatureSetEnrichmentState.defaults();
    params.collections.push("human-GO");
    await enrich.compute(params);

    let deets = enrich.fetchCollectionDetails();
    expect("human-GO" in deets).toBe(true);
    let res = enrich.fetchGroupResults(1, "cohen", "min_rank");
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
    let scores = enrich.fetchPerCellScores("human-GO", chosen);
    expect(scores.weights.length).toEqual(human_sizes[chosen]);

    // Also works with blocking.
    let block = new Int32Array(loaded.matrix.numberOfColumns());
    for (var i = 0; i < block.length; i++) {
        block[i] = i % 2 == 0;
    }

    {
        let enrich2 = new bakana.FeatureSetEnrichmentStandalone(loaded.features.RNA, markers, { normalized: normed, block: block });
        await enrich2.compute(params);
        let scores2 = enrich2.fetchPerCellScores("human-GO", chosen);
        expect(scores2).not.toEqual(scores);
    }

    // Throws an error correctly if no normalized matrix is supplied.
    {
        let enrich2 = new bakana.FeatureSetEnrichmentStandalone(loaded.features.RNA, markers);
        await enrich2.compute(params);
        expect(() => enrich2.fetchPerCellScores("human-GO", chosen)).toThrow("no normalized matrix");
    }

    markers.free();
    normed.free();
})

