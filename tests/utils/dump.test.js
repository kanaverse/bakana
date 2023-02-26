import * as scran from "scran.js";
import * as fs from "fs";
import * as abnode from "../../src/dump/abstract/dump_node.js";
import * as abbrowser from "../../src/dump/abstract/dump_node.js";

test("all MD5 sum calculations are consistent", async () => {
    let test = new Uint8Array(100000);
    test.forEach((x, i) => {
        test[i] = Math.random() * 255;
    });

    let nfiles = [ { contents: test, metadata: {} } ];
    await abnode.attachMd5sums(nfiles);
    expect(typeof nfiles[0].metadata.md5sum).toBe("string");

    let nfiles2 = [ { contents: test, metadata: {} } ];
    await abbrowser.attachMd5sums(nfiles2);
    expect(nfiles2[0].metadata.md5sum).toBe(nfiles[0].metadata.md5sum);

    let tmp = scran.chooseTemporaryPath();
    fs.writeFileSync(tmp, test);
    let nfiles3 = [ { contents: tmp, metadata: {} } ];
    await abnode.attachMd5sums(nfiles3);
    expect(nfiles3[0].metadata.md5sum).toBe(nfiles[0].metadata.md5sum);
})
