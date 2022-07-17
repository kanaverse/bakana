import * as utils from "../../src/readers/utils/extract.js";

test("summarizeValues works as expected for categorical variables", () => {
    let output = utils.summarizeValues(["A", "B", "A", "C"], 5);
    expect(output.type).toBe("categorical");
    expect(output.values).toEqual(["A", "B", "C"]);
    expect(output.truncated).toBe(false);

    // Works with truncation.
    output = utils.summarizeValues(["C", "D", "E", "B", "F", "A", "C", "F"], 5);
    expect(output.type).toBe("categorical");
    expect(output.values).toEqual(["A", "B", "C", "D", "E"]);
    expect(output.truncated).toBe(true);
})

test("summarizeValues works as expected for continuous variables", () => {
    let y = new Float64Array([1,2,3.5,-5]);
    let output = utils.summarizeValues(y);
    expect(output.type).toBe("continuous");
    expect(output.min).toBe(-5);
    expect(output.max).toBe(3.5);

    // Works with a few Infs and NaNs involved.
    y = new Float64Array([1,2,3.5,-Infinity]);
    output = utils.summarizeValues(y);
    expect(output.min).toBe(-Infinity);
    expect(output.max).toBe(3.5);

    y[3] = Infinity;
    output = utils.summarizeValues(y);
    expect(output.min).toBe(1);
    expect(output.max).toBe(Infinity);

    y[3] = NaN;
    output = utils.summarizeValues(y);
    expect(output.min).toBe(1);
    expect(output.max).toBe(3.5);
})

