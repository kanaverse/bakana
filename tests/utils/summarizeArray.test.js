import * as utils from "../../src/readers/utils/extract.js";

test("summarizeArray works as expected for categorical variables", () => {
    let output = utils.summarizeArray(["A", "B", "A", "C"], 5);
    expect(output.type).toBe("categorical");
    expect(output.values).toEqual(["A", "B", "C"]);
    expect(output.truncated).toBe(false);

    // Works with truncation.
    output = utils.summarizeArray(["C", "D", "E", "B", "F", "A", "C", "F"], { limit: 5 });
    expect(output.type).toBe("categorical");
    expect(output.values).toEqual(["A", "B", "C", "D", "E"]);
    expect(output.truncated).toBe(true);
})

test("summarizeArray works as expected for continuous variables", () => {
    let y = new Float64Array([1,2,3.5,-5]);
    let output = utils.summarizeArray(y);
    expect(output.type).toBe("continuous");
    expect(output.min).toBe(-5);
    expect(output.max).toBe(3.5);

    // Works with a few Infs and NaNs involved.
    y = new Float64Array([1,2,3.5,-Infinity]);
    output = utils.summarizeArray(y);
    expect(output.min).toBe(-Infinity);
    expect(output.max).toBe(3.5);

    y[3] = Infinity;
    output = utils.summarizeArray(y);
    expect(output.min).toBe(1);
    expect(output.max).toBe(Infinity);

    y[3] = NaN;
    output = utils.summarizeArray(y);
    expect(output.min).toBe(1);
    expect(output.max).toBe(3.5);

    // Works if there's only NaNs.
    y = new Float64Array([NaN, NaN, NaN]);
    output = utils.summarizeArray(y);
    expect(output.min > output.max).toBe(true);
})

