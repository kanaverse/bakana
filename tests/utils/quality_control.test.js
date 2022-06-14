import * as bakana from "./../../src/index.js";
import * as qcutils from "./../../src/utils/quality_control.js";

beforeAll(async () => await bakana.initialize({ localFile: true }));
afterAll(async () => await bakana.terminate());

test("splitting metrics by block works as expected", () => {
    let a1 = bakana.callScran(module => module.createFloat64WasmArray(4));
    a1.set([1,2,3,4]);
    let a2 = bakana.callScran(module => module.createInt32WasmArray(4));
    a2.set([5,6,7,8]);

    let metrics = { "foo": a1, "bar": a2 };
    let blocks = [ "odd", "even" ];
    let block_ids = bakana.callScran(module => module.createInt32WasmArray(4));
    block_ids.set([0,1,0,1]);

    let split = qcutils.splitMetricsByBlock(metrics, blocks, block_ids);
    expect(split.odd.foo instanceof Float64Array).toBe(true);
    expect(split.odd.bar instanceof Int32Array).toBe(true);
    expect(Array.from(split.odd.foo)).toEqual([1,3]);
    expect(Array.from(split.odd.bar)).toEqual([5,7]);
    expect(Array.from(split.even.foo)).toEqual([2,4]);
    expect(Array.from(split.even.bar)).toEqual([6,8]);
})

test("splitting thresholds by block works as expected", () => {
    let a1 = bakana.callScran(module => module.createFloat64WasmArray(2));
    a1.set([1.2,3.4]);
    let a2 = bakana.callScran(module => module.createFloat64WasmArray(2));
    a2.set([5.6,7.8]);

    let thresholds = { "foo": a1, "bar": a2 };
    let blocks = [ "odd", "even" ];

    let split = qcutils.splitThresholdsByBlock(thresholds, blocks);
    expect(split.odd.foo).toEqual(1.2);
    expect(split.odd.bar).toEqual(5.6);
    expect(split.even.foo).toEqual(3.4);
    expect(split.even.bar).toEqual(7.8);
})
