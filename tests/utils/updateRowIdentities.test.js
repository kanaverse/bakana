import * as input from "../../src/steps/inputs.js";

test("updating the row identities works for 1:1 mappings", () => {
    {
        let fun = input.updateRowIdentities([1,2,3], [3,2,1]);
        expect(fun(["aaron", "alice", "bob"])).toEqual(["bob", "alice", "aaron"]);
        expect(fun(new Int32Array([8, 5, 10]))).toEqual(new Int32Array([10, 5, 8]));
    }

    {
        let fun = input.updateRowIdentities([400, 300, 200, 100], [300, 400, 100, 200]);
        expect(fun(["alicia", "athena", "akira", "aria"])).toEqual(["athena", "alicia", "aria", "akira"]);
        expect(fun(new Float64Array([-1, -2, -3, -4]))).toEqual(new Float64Array([-2, -1, -4, -3]));
    }
})

test("updating the row identities works for no-ops", () => {
    let noop = input.updateRowIdentities([1,2,3], [1,2,3]);
    let str_input = ["aaron", "alice", "bob"];
    expect(noop(str_input)).toBe(str_input); // use toBe as this should be exactly the same object.

    let finput = new Float64Array([1.2, 2.3, -1]);
    expect(noop(finput)).toBe(finput);
})

test("updating the row identities handles missing values", () => {
    let missing = input.updateRowIdentities([1,2,3,4,5,6], [2,4,6]);
    expect(missing(["Yuuko", "Mai", "Mio"])).toEqual([null, "Yuuko", null, "Mai", null, "Mio"]);
    expect(missing(new Float64Array([201, 302, 404]))).toEqual(new Float64Array([Number.NaN, 201, Number.NaN, 302, Number.NaN, 404]));
    expect(missing(new Int32Array([20, 40, 60]))).toEqual(new Int32Array([-1, 20, -1, 40, -1, 60]));
})
