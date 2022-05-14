import * as bakana from "./../../src/index.js";
import * as inputs from "./../../src/inputs.js";
import * as qc from "./../../src/adt/quality_control.js";

beforeAll(async () => await bakana.initialize({ localFile: true }));
afterAll(async () => await bakana.terminate());

test("quality control works as expected for ADTs", async () => {
    let istate = new inputs.InputsState;
    await istate.compute({
        default: {
            format: "MatrixMarket",
            mtx: "files/datasets/immune_3.0.0-matrix.mtx.gz",
            genes: "files/datasets/immune_3.0.0-features.tsv.gz",
            annotations: "files/datasets/immune_3.0.0-barcodes.tsv.gz"
        }
    }, null);

    let qcstate = new qc.AdtQualityControlState(istate);
    qcstate.compute("IgG", 3, 0.1);
    let summ = qcstate.summary();

    // Checking that that the metrics are sensibly computed.
    let positive_total = 0;
    summ.data.default.igg_total.forEach(x => { positive_total += (x > 0); });
    expect(positive_total).toBeGreaterThan(0);

    // Checking that the thresholds are sensible.
    expect(summ.thresholds.default.detected).toBeGreaterThan(0);
    expect(summ.thresholds.default.igg_total).toBeGreaterThan(0);

    istate.free();
    qcstate.free();
})

test("quality control skips in the absence of ADTs", async () => {
    let istate = new inputs.InputsState;
    await istate.compute({
        default: {
            format: "10X",
            h5: "files/datasets/pbmc4k-tenx.h5"
        }
    }, null);

    let qcstate = new qc.AdtQualityControlState(istate);
    qcstate.compute("IgG", 3, 0.1);

    let summ = qcstate.summary();
    expect(summ).toBeNull();
})

test("quality control works as expected for ADTs with blocking", async () => {
    let istate = new inputs.InputsState;
    await istate.compute({
        "combined": {
            format: "MatrixMarket",
            mtx: "files/datasets/pbmc-combined-matrix.mtx.gz",
            genes: "files/datasets/pbmc-combined-features.tsv.gz",
            annotations: "files/datasets/pbmc-combined-barcodes.tsv.gz"
        }
    }, "3k");

    // Let's pretend the main count matrix is the ADT matrix.
    let qcstate = new qc.AdtQualityControlState(istate);
    qcstate.useMainMatrix();
    qcstate.compute("IgG", 3, 0.1);
    let summ = qcstate.summary();

    // Checking that that there are multiple metrics.
    let positive_total = 0;
    summ.data["3k"].detected.forEach(x => { positive_total += (x > 0); });
    expect(positive_total).toBeGreaterThan(0);

    // This test should not have any IgGs.
    let zero_total = 0;
    summ.data["4k"].igg_total.forEach(x => { zero_total += (x > 0); });
    expect(zero_total).toBe(0);

    // Checking that the thresholds are sensible.
    expect(summ.thresholds["3k"].detected).toBeGreaterThan(0);
    expect(summ.thresholds["4k"].detected).toBeGreaterThan(0);

    istate.free();
    qcstate.free();
})
