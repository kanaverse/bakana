import * as from_v0 from "./../legacy/from_v0.js";
import * as scran from "scran.js";
import * as pako from "pako";

// Must be integers!
const FORMAT_EMBEDDED = 0;
const FORMAT_LINKED = 1;
const FORMAT_VERSION = 1001000;

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

export function embedFiles() {
    let output = { collected: [], total: 0 };
    output.saver = (obj) => {
        output.collected.push(obj.buffer);
        let current = output.total;
        let size = obj.buffer.byteLength;
        output.total += size;
        return {
            "offset": current,
            "size": size
        };
    }
    return output;
}

function save_internal(format_type, state, extras) {
    var combined = new ArrayBuffer(24 + state.length + extras);
    var combined_arr = new Uint8Array(combined);
    var offset = 0;

    let format = numberToBuffer(format_type);
    combined_arr.set(format, offset); 
    offset += format.length;

    let version = numberToBuffer(FORMAT_VERSION);
    combined_arr.set(version, offset); 
    offset += version.length;

    let state_len = numberToBuffer(state.length);
    combined_arr.set(state_len, offset); 
    offset += state_len.length;

    if (offset != 24) {
        throw "oops - accounting error in the serialization code!";
    }

    combined_arr.set(state, offset);
    offset += state.length;

    return {
        "offset": offset,
        "combined": combined                
    }
}

export function createEmbeddedKanaFile(state, collected) {
    let total_len = 0;
    for (const buf of collected) {
        total_len += buf.byteLength;
    }

    let saved = save_internal(FORMAT_EMBEDDED, state, total_len);
    let offset = saved.offset;
    let combined_arr = new Uint8Array(saved.combined);

    for (const buf of collected) {
        const tmp = new Uint8Array(buf);
        combined_arr.set(tmp, offset);
        offset += tmp.length;
    }

    return saved.combined;
}

export function createLinkedKanaFile(state) {
    let saved = save_internal(FORMAT_LINKED, state, 0);
    return saved.combined;
}

export async function loadKanaFile(buffer, statePath) {
    var offset = 0;
    var format = bufferToNumber(new Uint8Array(buffer, offset, 8));
    offset += 8;

    var version = bufferToNumber(new Uint8Array(buffer, offset, 8));
    offset += 8;

    var state_len = bufferToNumber(new Uint8Array(buffer, offset, 8));
    offset += 8;

    let state = new Uint8Array(buffer, offset, state_len);
    offset += state_len;
    if (version < 1000000) {
        let contents = pako.ungzip(state, { "to": "string" });
        let values = JSON.parse(contents);
        from_v0.convertFromVersion0(values, statePath);
    } else {
        scran.writeFile(statePath, state);
    }

    let bundle = {};
    if (format == FORMAT_EMBEDDED) {
        bundle.remaining = new Uint8Array(buffer, offset, buffer.byteLength - offset);
        bundle.loader = (start, size) => bundle.remaining.slice(start, start + size);
        bundle.embedded = true;
    } else if (format == FORMAT_LINKED) {
        bundle.embedded = false;
    } else {
        throw "unsupported format type";
    }

    return bundle;
}
