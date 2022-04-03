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
    constructor(obj) {
        let reader = new FileReaderSync();
        this.buffer = reader.readAsArrayBuffer(obj);
    }

    buffer() {
        return this.buffer;
    }

    size() {
        return this.buffer.byteLength;
    }

    serialized() {
        return this.buffer;
    }
};

export function removeH5(path) {
    if (scran.fileExists(path)) {
        scran.removeFile(path);
    }
    return;
}

export function realizeH5(loaded) {  
    const tmppath = "temp_" + String(Number(new Date())) + ".h5";
    scran.writeFile(tmppath, new Uint8Array(loaded.buffer));
    return tmppath;
}