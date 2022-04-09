import * as sutils from "../utils/serialize.js";
import * as scran from "scran.js";
import * as v0 from "../legacy/from_v0.js";
import * as pako from "pako";

/**
 * This contains a function to create and load a kana file with the browser.
 */
export function createKanaFileInternal(statePath, inputFiles) {
    let embedded = (inputFiles !== null);
    let state = scran.readFile(statePath);

    let preamble = sutils.createPreamble(embedded, state.byteLength);

    let total = preamble.byteLength + state.byteLength;
    if (embedded) {
        for (const ibuf of inputFiles) {
            total += ibuf.byteLength;
        }
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

export function parseKanaFileInternal(input, statePath) {
    let parsed = sutils.parsePreamble(input);
    let delta = parsed.offset + parsed.state;
    let statebuffer = input.slice(parsed.offset, delta);

    if (parsed.version < 1000000) {
        var contents = pako.ungzip(new Uint8Array(statebuffer), { "to": "string" });
        let state = JSON.parse(contents);
        v0.convertFromVersion0(state, statePath);
    } else {
        scran.writeFile(statePath, new Uint8Array(statebuffer));
    }

    if (parsed.embedded) {
        return (offset, size) => input.slice(delta + offset, delta + offset + size);
    } else {
        return null;
    }
}
