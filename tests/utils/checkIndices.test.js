import * as utils from "../../src/steps/utils/general.js";

test("checkIndices fails as expected", () => {
    utils.checkIndices([0,1,2,3], 4); // fine.

    expect(() => utils.checkIndices([0,1,2,3], 3)).toThrow("out of range");
    expect(() => utils.checkIndices([-1,1,2], 3)).toThrow("out of range");

    expect(() => utils.checkIndices([0,2,1], 3)).toThrow("sorted");
    expect(() => utils.checkIndices([0,1,1,2], 3)).toThrow("unique");
})
