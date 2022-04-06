import * as sutils from "../utils/serialize.js";
import * as scran from "scran.js";

/**
 * This contains a function to create and load a kana file with the browser.
 */
export async function createKanaFileInternal(path, statePath, inputFiles) {
    let embedded = (inputFiles === null):
    let state = scran.readFile(statePath);

    let preamble = sutils.createPreamble(embedded, state.byteLength);

    let total = preamble.length + state.byteLength;
    for (const ibuf of inputFiles) {
        total += ibuf.byteLength;
    }

    let output = new ArrayBuffer(total);
    let arr = new Uint8Array(output);
    arr.set(new Uint8Array(preamble));

    let offset = preamble.byteLength;
    arr.set(new Uint8Array(state), offset);
    offset += state.byteLength;

    if (embedded) {
        for (const ibuf of inputFiles) {
            arr.set(new Uint8Array(ibuf), offset);
            offset += ibuf.byteLength;
        }
    }

    return output;
}

export async function parseKanaFileInternal(input, statePath) {
    const reader = new FileReaderSync;
    let buffer = reader.readAsArrayBuffer(input);
    let arr = new Uint8Array(buffer);

    let parsed = sutils.parsePreamble(buffer);
    let delta = parsed.offset + parsed.state;
    let statebuffer = buffer.slice(parsed.offset, delta);
    scran.writeFile(statePath, statebuffer);

    if (parsed.embedded) {
        return (offset, size) => buffer.slice(delta + offset, delta + offset + size);
    } else {
        return null;
    }
}
