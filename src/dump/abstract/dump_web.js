import { md5 } from 'hash-wasm';

export async function attachMd5sums(files) {
    for (const x of files) {
        x.metadata.md5sum = await md5(x.contents);
    }
}

export function realizeDirectory(files, directory) {
    throw new Error("cannot realize files into a directory in a web context");
    return;
}

export function loadFilePath(p) {
    throw new Error("cannot load file paths in a web context");
}
