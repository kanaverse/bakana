/**
 * 'obj' is expected to be a File object.
 */

import * as scran from "scran.js";

export function rawSize(obj) {
    return obj.size;
}

export function rawName(obj) {
    return obj.name;
}

export class LoadedFile {
    #buffer;

    constructor(obj) {
        if (obj instanceof File) {
            let reader = new FileReaderSync();
            this.#buffer = reader.readAsArrayBuffer(obj);
        } else if (obj instanceof ArrayBuffer) {
            this.#buffer = obj; // assumed to already be an ArrayBuffer.
        } else {
            throw "unknown type '" + typeof(obj) + "' for LoadedFile constructor";
        }
    }

    buffer() {
        return this.#buffer;
    }

    size() {
        return this.#buffer.byteLength;
    }

    serialized() {
        return this.#buffer;
    }
};

export function realizeMatrixMarket(loaded) {
    return new Uint8Array(loaded.buffer());
}

export function realizeH5(loaded) {
    let tmppath;

    do {
        tmppath = "temp_" + String(Number(new Date())) + "_" + String(Math.round(Math.random() * 10000)) + ".h5";
    } while (scran.fileExists(tmppath));

    scran.writeFile(tmppath, new Uint8Array(loaded.buffer()));
    return tmppath;
}

// The default is allowed = true here; we are always allowed to delete stuff
// from the VFS, given that we created it with realizeH5. 
export function removeH5(path, { allowed = true } = {}) {
    if (scran.fileExists(path)) {
        scran.removeFile(path);
    }
    return;
}
