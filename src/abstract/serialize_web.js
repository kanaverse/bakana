import * as sutils from "./utils/serialize.js";
import * as scran from "scran.js";
export { FORMAT_VERSION } from "./utils/serialize.js";

export function createKanaFileInternal(statePath, inputFiles) {
    let embedded = (inputFiles !== null);
    let state = scran.readFile(statePath);

    let preamble = sutils.createPreamble(embedded, state.byteLength);
    let total = preamble.length + state.length;
    if (embedded) {
        for (const ibuf of inputFiles) {
            total += ibuf.byteLength;
        }
    }

    let output = new Uint8Array(output);
    output.set(preamble);

    let offset = preamble.length;
    output.set(state, offset);
    offset += state.length;

    if (embedded) {
        for (const ibuf of inputFiles) {
            arr.set(ibuf, offset);
            offset += ibuf.length;
        }
    }

    return arr;
}

export function parseKanaFileInternal(input, statePath) {
    return sutils.parseKanaFileFromBuffer(input, statePath);
}
