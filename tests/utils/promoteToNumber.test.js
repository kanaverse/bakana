import * as utils from "../../src/readers/utils/extract.js";

test("promoteToNumber works as expected", () => {
    let promoted = utils.promoteToNumber(["1.0", "1e6", "-2.2E-5", "-2.0", "NaN", "NA", "Inf", "-Inf", "inf" ]);
    expect(promoted.constructor.name).toBe("Float64Array");

    expect(promoted[0]).toBe(1);
    expect(promoted[1]).toBe(1e6);
    expect(promoted[2]).toBe(-2.2e-5);
    expect(promoted[3]).toBe(-2);
    expect(promoted[4]).toBe(NaN);
    expect(promoted[5]).toBe(NaN);
    expect(promoted[6]).toBe(Infinity);
    expect(promoted[7]).toBe(-Infinity);
    expect(promoted[8]).toBe(Infinity);
})
