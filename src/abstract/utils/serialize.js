import * as from_v0 from "./legacy_v0.js";
import * as scran from "scran.js";
import * as pako from "pako";

// Must be integers!
const FORMAT_EMBEDDED = 0;
const FORMAT_LINKED = 1;
export const FORMAT_VERSION = 3000000;

function numberToBuffer(number) {
    // Store as little-endian. Probably safer
    // than trying to cast it from a Uint64Array;
    // not sure that endianness is strictly defined.
    var output = new Uint8Array(8);

    var i = 0;
    while (number > 0) {
        output[i] = number % 256;
        number = Math.floor(number / 256);
        i++;
    }

    return output;
}

function bufferToNumber(buffer) {
    var output = 0;
    var multiplier = 1;
    for (const x of buffer) {
        output += multiplier * x;
        multiplier *= 256;
    }
    return output;
}

export function createPreamble(embedded, stateSize) {
    var combined = new Uint8Array(24);
    var offset = 0;

    let format = numberToBuffer(embedded ? FORMAT_EMBEDDED : FORMAT_LINKED);
    combined.set(format, offset); 
    offset += format.length;

    let version = numberToBuffer(FORMAT_VERSION);
    combined.set(version, offset); 
    offset += version.length;

    let state_len = numberToBuffer(stateSize);
    combined.set(state_len, offset); 
    offset += state_len.length;

    if (offset != 24) {
        throw "oops - accounting error in the serialization code!";
    }

    return combined;
}

export function parsePreamble(buffer) {
    var offset = 0;
    var format = bufferToNumber(buffer.subarray(offset, offset + 8));
    offset += 8;

    var version = bufferToNumber(buffer.subarray(offset, offset + 8));
    offset += 8;

    var state_len = bufferToNumber(buffer.subarray(offset, offset + 8));
    offset += 8;

    return {
        "embedded": (format == FORMAT_EMBEDDED),
        "version": version,
        "state": state_len,
        "offset": offset
    };
}

export function parseKanaFileFromBuffer(input, statePath) {
    let parsed = parsePreamble(input);
    let delta = parsed.offset + parsed.state;
    let statebuffer = input.subarray(parsed.offset, parsed.offset + delta);

    if (parsed.version < 1000000) {
        var contents = pako.ungzip(statebuffer, { "to": "string" });
        let state = JSON.parse(contents);
        from_v0.convertFromVersion0(state, statePath);
    } else {
        scran.writeFile(statePath, statebuffer);
    }

    if (parsed.embedded) {
        // The buffer takes ownership of the embedded file to allow eventual
        // garbage collection of the input buffer.
        return (offset, size) => input.slice(delta + offset, delta + offset + size);
    } else {
        return null;
    }
}
