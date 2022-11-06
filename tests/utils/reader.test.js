import * as bakana from "../../src/index.js";
import * as pako from "pako";
import * as fs from "fs";

test("unpacking text lines works as expected", async () => {
    let arr = ["Aaron was here", "FOOBLE", "BAR"];
    let target = arr.join("\n");

    const foo = new TextEncoder;
    let out = foo.encode(target);
    let res = await bakana.readLines2(new Uint8Array(out));
    expect(res).toEqual(arr);

    // With the terminating newline.
    let out2 = foo.encode(target + "\n");
    let res2 = await bakana.readLines2(new Uint8Array(out2));
    expect(res2).toEqual(arr);

    // Plus some compression.
    let outgz = pako.gzip(target);
    let resgz = await bakana.readLines2(new Uint8Array(outgz));
    expect(resgz).toEqual(arr);

    let resgz2 = await bakana.readLines2(new Uint8Array(outgz), { compression: "gz" }); // explicitly
    expect(resgz2).toEqual(arr);
})

test("unpacking Gzipped text lines works with several chunk sizes", async () => {
    let arr = [
        "Kimi dake o kimi dake o", "Suki de ita yo", "Kaze de me ga nijinde", "Tooku naru yo",
        "Itsu made mo oboeteru", "Nani mo ka mo kawatte mo", "Hitotsu dake hitotsu dake", 
        "Arifureta mono da kedo", "Misete yaru kagayaki ni michita sono hitotsu dake", 
        "Itsu made mo itsu made mo mamotte iku"
    ];
    let target = arr.join("\n");

    let outgz = pako.gzip(target);
    let resgz10 = await bakana.readLines2(new Uint8Array(outgz), { chunkSize: 10 });
    expect(resgz10).toEqual(arr);

    let resgz20 = await bakana.readLines2(new Uint8Array(outgz), { chunkSize: 20 });
    expect(resgz20).toEqual(arr);

    let resgz50 = await bakana.readLines2(new Uint8Array(outgz), { chunkSize: 50 });
    expect(resgz50).toEqual(arr);
})

test("unpacking text lines works from file", async () => {
    let lines = [];
    for (var i = 0; i < 100000; i++) {
        lines.push(String(Math.random()));
    }
    let target = lines.join("\n");

    {
        let path = "TEST_unpack_lines.txt";
        fs.writeFileSync(path, target);
        let res = await bakana.readLines2(path);
        expect(res).toEqual(lines);
    }

    // Trying again with gzip.
    {
        let path = "TEST_unpack_lines.txt.gz";
        fs.writeFileSync(path, pako.gzip(target));
        let res = await bakana.readLines2(path);
        expect(res).toEqual(lines);
    }
})

test("unpacking DSV buffer works as expected", async () => {
    let target = "a,b,c\n1,2,3\n\"x\",\"y\",\"z\""; 
    let expected = [["a","b","c"],["1","2","3"],["x","y","z"]];

    const foo = new TextEncoder;
    let out = foo.encode(target);
    let res = await bakana.readTable2(new Uint8Array(out), {delim:","});
    expect(res).toEqual(expected);

    // Removes the trailing newline.
    let out2 = foo.encode(target + "\n");
    let res2 = await bakana.readTable2(new Uint8Array(out2), {delim:","});
    expect(res2).toEqual(expected);
})

test("unpacking DSV file works as expected", async () => {
    let lines = [];
    let expected = [];
    for (var i = 0; i < 100000; i++) {
        let current = [ Math.random(), Math.random(), Math.random() ];
        expected.push(current.map(String));
        lines.push(current.join("\t"));
    }

    let target = lines.join("\n");

    {
        let path = "TEST_unpack_lines.txt";
        fs.writeFileSync(path, target);
        let res = await bakana.readTable2(path);
        expect(res).toEqual(expected);
    }

    // With a smaller buffer
    {
        let path = "TEST_unpack_lines.txt";
        fs.writeFileSync(path, target);
        let res = await bakana.readTable2(path, { chunkSize: 65536 });
        expect(res).toEqual(expected);
    }

    // Trying again with gzip.
    {
        let path = "TEST_unpack_lines.txt.gz";
        fs.writeFileSync(path, pako.gzip(target));
        let res = await bakana.readTable2(path);
        expect(res).toEqual(expected);
    }
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
