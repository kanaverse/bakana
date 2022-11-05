import * as fs from "fs";

export function stream(file, chunkSize, callback, resolve, reject) {
    let x = fs.createReadStream(file, { highWaterMark: chunkSize });
    x.on("data", callback);
    if (resolve !== null) {
        x.on("end", () => resolve(null));
    }
    x.on('error', e => reject("file streaming failed; " + e.message));
}

export function peek(file, n) {
    let handle = fs.openSync(file, "r");
    let output = new Uint8Array(n);
    let read = fs.readSync(handle, output, 0, n);
    return output.slice(0, read);
}
