import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

let h5path = "files/datasets/pbmc4k-tenx.h5";
let files = {
    default: new bakana.TenxHdf5Dataset(h5path)
};

test("Standalone marker detection works correctly", async () => {
    let loaded = files.default.load({ cache: true });
    let normed = {};
    for (const m of loaded.matrix.available()) {
        normed[m] = scran.logNormCounts(loaded.matrix.get(m));
    }

    let groups = new Int32Array(loaded.matrix.numberOfColumns());
    let half = groups.length / 2
    groups.fill(0, 0, half);
    groups.fill(1, half, groups.length);

    // Checks kick in.
    expect(() => new bakana.MarkerDetectionStandalone(normed, ["A"])).toThrow("non-negative integer");
    expect(() => new bakana.MarkerDetectionStandalone(normed, [])).toThrow("same number of columns");

    let markers = new bakana.MarkerDetectionStandalone(normed, groups);
    markers.compute(bakana.MarkerDetectionStandalone.defaults());
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
        markers2.compute(bakana.MarkerDetectionStandalone.defaults());

        let res2 = markers2.fetchResults();
        expect(res2["RNA"].cohen(0)).not.toEqual(res["RNA"].cohen(0));

        markers2.free();
    }

    markers.free();
})

test("Standalone custom selections work correctly", async () => {
    let loaded = files.default.load({ cache: true });
    let normed = {};
    for (const m of [ "foo", "bar" ]) { // mimicking multiple modalities.
        normed[m] = scran.logNormCounts(loaded.matrix.get("RNA"));
    }

    let custom = new bakana.CustomSelectionsStandalone(normed);
    custom.compute(bakana.CustomSelectionsStandalone.defaults());
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
        custom2.compute(bakana.CustomSelectionsStandalone.defaults());
        custom2.addSelection("WHEEE", [ 1,2,3,4,5,6,7,8,9,10 ]);

        let res2 = custom2.fetchResults("WHEEE");
        expect(res2["foo"].cohen(0)).not.toEqual(res["foo"].cohen(0));

        custom2.free();
    }

    custom.free();
})

