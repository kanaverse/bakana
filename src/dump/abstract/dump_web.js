import { md5 } from 'hash-wasm';
import * as scran from "scran.js";

export async function attachMd5sums(files) {
    for (const x of files) {
        if (!("contents" in x)) {
            continue;
        }
        x.metadata.md5sum = await md5(x.contents);
    }
}

export function realizeDirectory(files, directory) {
    throw new Error("cannot realize files into a directory in a web context");
    return;
}

export function loadFilePath(p) {
    return scran.readFile(p);
}
