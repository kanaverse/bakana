import * as bakana from "../../src/index.js";
import * as utils from "../../src/steps/utils/general.js";

beforeAll(async () => await bakana.initialize({ localFile: true }));
afterAll(async () => await bakana.terminate());

test("changedParameters works as expected for primitive types", () => {
    expect(utils.changedParameters(null, null)).toBe(false);
    expect(utils.changedParameters("A", null)).toBe(true);
    expect(utils.changedParameters(1, 2)).toBe(true);
    expect(utils.changedParameters(1, 1)).toBe(false);
    expect(utils.changedParameters(1, "1")).toBe(true); // type sensitive.
})

test("changedParameters works as expected for arrays", () => {
    expect(utils.changedParameters([1,2,3], [1,2,3])).toBe(false);
    expect(utils.changedParameters(1, [1,2,3])).toBe(true);
    expect(utils.changedParameters({}, [1,2,3])).toBe(true);
    expect(utils.changedParameters(null, [1,2,3])).toBe(true);
    expect(utils.changedParameters([1,3], [1,2,3])).toBe(true);
    expect(utils.changedParameters([2,3,1], [1,2,3])).toBe(true);

    // But fails for typedArrays.
    expect(() => utils.changedParameters(new Int32Array([1,2,3]), "ALPHABET")).toThrow("cannot contain");
    expect(() => utils.changedParameters(new Int32Array([1,2,3]), [1,2,3])).toThrow("cannot contain");
    expect(() => utils.changedParameters([1,2,3], new ArrayBuffer(24))).toThrow("cannot contain");
    expect(() => utils.changedParameters(null, new ArrayBuffer(24))).toThrow("cannot contain");
    expect(() => utils.changedParameters(new Int32Array([1,2,3]), { a: 2, b: 3})).toThrow("cannot contain");
})

test("changedParameters works as expected for objects", () => {
    expect(utils.changedParameters({a:1}, {a:1})).toBe(false);
    expect(utils.changedParameters({a:1, b:2}, {b:2, a:1})).toBe(false); // order doesn't matter.

    expect(utils.changedParameters({a:1}, {b:2, a:1})).toBe(true); 
    expect(utils.changedParameters({a:1}, {a:3})).toBe(true);
    expect(utils.changedParameters({a:1}, {A:3})).toBe(true);

    // still works with deep nesting.
    expect(utils.changedParameters({a:{b:{c:3}}}, {a:{b:{c:4}}})).toBe(true); 
})

test("changedParameters works as expected for subsets", () => {
    let subset = { field: "WHEE", values: ["A"] };
    expect(utils.changedParameters(subset, subset)).toBe(false);
    expect(utils.changedParameters(subset, null)).toBe(true);
    expect(utils.changedParameters(subset, { field: "BLAH", values: ["A"]})).toBe(true);
    expect(utils.changedParameters(subset, { field: "WHEE", values: ["A", "B"]})).toBe(true);
    expect(utils.changedParameters(subset, { field: "WHEE", values: ["B"]})).toBe(true);

    let subset2 = { field: "WHEE", ranges: [[0,1],[2,3]] };
    expect(utils.changedParameters(subset2, subset2)).toBe(false);
    expect(utils.changedParameters(subset2, subset)).toBe(true);
    expect(utils.changedParameters(subset2, null)).toBe(true);
    expect(utils.changedParameters(subset2, { field: "BLAH", ranges:[[0,1],[2,3]]})).toBe(true);
    expect(utils.changedParameters(subset2, { field: "WHEE", ranges:[[0,1]]})).toBe(true);
    expect(utils.changedParameters(subset2, { field: "WHEE", ranges:[[0,1],[2,4]]})).toBe(true);
})

test("changedParameters works as expected for weights", () => {
    expect(utils.changedParameters({RNA:1}, null)).toBe(true);
    expect(utils.changedParameters({RNA:1}, {RNA:1})).toBe(false);
    expect(utils.changedParameters({RNA:1}, {RNA:2})).toBe(true);

    expect(utils.changedParameters({RNA:1, ADT:2}, null)).toBe(true);
    expect(utils.changedParameters({RNA:1, ADT:2}, {RNA:1})).toBe(true);
    expect(utils.changedParameters({RNA:1, ADT:2}, {RNA:1,ADT:2})).toBe(false);
    expect(utils.changedParameters({RNA:1, ADT:2}, {RNA:1,ADT:0.5})).toBe(true);
})

test("changedParameters works as expected for input abbreviations", () => {
    let dump = {
        MatrixMarket: new bakana.TenxMatrixMarketDataset(
                "files/datasets/pbmc3k-matrix.mtx.gz",
                "files/datasets/pbmc3k-features.tsv.gz",
                "files/datasets/pbmc3k-barcodes.tsv.gz"
            ),
        tenx: new bakana.TenxHdf5Dataset("files/datasets/pbmc4k-tenx.h5"),
        h5ad: new bakana.H5adDataset("files/datasets/zeisel-brain.h5ad")
    };

    let abbreviations = {};
    for (const [k, v] of Object.entries(dump)) {
        abbreviations[k] = v.abbreviate();
    }

    expect(utils.changedParameters(abbreviations.MatrixMarket, abbreviations.MatrixMarket)).toBe(false);
    expect(utils.changedParameters(abbreviations.MatrixMarket, abbreviations.tenx)).toBe(true);
    expect(utils.changedParameters(abbreviations.MatrixMarket, abbreviations.h5ad)).toBe(true);
    expect(utils.changedParameters(abbreviations.tenx, abbreviations.tenx)).toBe(false);
    expect(utils.changedParameters(abbreviations.tenx, abbreviations.h5ad)).toBe(true);
    expect(utils.changedParameters(abbreviations.h5ad, abbreviations.h5ad)).toBe(false);
})
