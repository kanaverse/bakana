/**
 * 'obj' is expected to be a string containing the file path.
 */

import * as fs from "fs";

export function size(obj) {
    return fs.statSync(obj).size;
}

export function buffer(obj) {
    let f = fs.readFileSync(obj);
    return f.buffer.slice(f.byteOffset, f.byteOffset + f.byteLength);
}
