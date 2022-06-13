import * as bakana from "../../src/index.js";
import * as pako from "pako";

test("unpacking text works as expected", () => {
    const foo = new TextEncoder;
    let target = "Aaron was here";
    let out = foo.encode(target);
    let res = bakana.unpackText(new Uint8Array(out));
    expect(res).toBe("Aaron was here");

    // With compression.
    let gz = pako.gzip(target);
    let res2 = bakana.unpackText(new Uint8Array(gz));
    expect(res2).toBe("Aaron was here");

    // Explicitly:
    let res3 = bakana.unpackText(new Uint8Array(gz), { compression: "gz" });
    expect(res3).toBe("Aaron was here");
})

test("unpacking text lines works as expected", () => {
    let arr = ["Aaron was here", "FOOBLE", "BAR"];
    let target = arr.join("\n");

    const foo = new TextEncoder;
    let out = foo.encode(target);
    let res = bakana.readLines(new Uint8Array(out));
    expect(res).toEqual(arr);

    // With the terminating newline.
    let out2 = foo.encode(target + "\n");
    let res2 = bakana.readLines(new Uint8Array(out2));
    expect(res2).toEqual(arr);
})

test("unpacking DSV files works as expected", () => {
    let target = "a,b,c\n1,2,3\n\"x\",\"y\",\"z\""; 
    let expected = [["a","b","c"],["1","2","3"],["x","y","z"]];

    const foo = new TextEncoder;
    let out = foo.encode(target);
    let res = bakana.readTable(new Uint8Array(out), {delim:","});
    expect(res).toEqual(expected);

    // Removes the trailing newline.
    let out2 = foo.encode(target + "\n");
    let res2 = bakana.readTable(new Uint8Array(out2), {delim:","});
    expect(res2).toEqual(expected);
})

test("promoting numbers works as expected", () => {
    let out = bakana.promoteToNumber(["1","2","3"]);
    expect(out instanceof Float64Array).toBe(true);
    expect(Array.from(out)).toEqual([1,2,3]);

    let out2 = bakana.promoteToNumber(["1","2","na",""]);
    expect(out2 instanceof Float64Array).toBe(true);
    expect(out2[2]).toBe(NaN);
    expect(out2[3]).toBe(NaN);

    let out3 = bakana.promoteToNumber(["1","2","foo"]);
    expect(out3 instanceof Float64Array).toBe(false);
})
