import * as pako from "pako";
import ppp from "papaparse";
import * as astream from "./abstract/stream.js";
import * as afile from "../abstract/file.js";

export function extractHdf5Strings(handle, name) {
    if (!(name in handle.children)) {
        return null;
    }

    if (handle.children[name] !== "DataSet") {
        return null;
    }

    let content = handle.open(name);
    if (content.type !== "String") {
        return null;
    }

    return content.load();
}

/**
 * Summarize an array, typically corresponding to a single column of per-cell annotation.
 * This can be used as part of a preflight response in a Reader.
 *
 * @param {Array|TypedArray} array - Per-cell annotation array of length equal to the number of cells for a given matrix.
 * An Array is treated as categorical data and should contain strings, while TypedArrays are treated as continuous data.
 * @param {object} [options] - Optional parameters.
 * @param {number} [options.limit=50] - Maximum number of unique values to report for categorical `x`.
 *
 * @return {object} Object containing `type`, a string indicating whether `array` was categorical or continuous.
 *
 * If `"categorical"`, the object will contain `values`, an array of unique values up to the length specified by `limit`.
 * It will also contain `truncated`, a boolean indicating whether the actual number of unique values exceeds `limit`.
 *
 * If `"continuous"`, the object will contain the numbers `min` and `max` specifying the minimum and maximum value in `x`, respectively.
 * `min` or `max` may be negative or positive infinity, respectively, if there is no bound on one or both ends.
 * If `min > max`, all values in `array` are `NaN`s such that no bound can be found.
 */
export function summarizeArray(array, { limit = 50 } = {}) {
    if (array instanceof Array) {
        let chosen = Array.from(new Set(array));
        chosen.sort();
        let truncated = false;
        if (chosen.length > limit) {
            chosen = chosen.slice(0, limit);
            truncated = true;
        }
        return {
            "type": "categorical",
            "values": chosen,
            "truncated": truncated
        };
    } else {
        let min = Number.POSITIVE_INFINITY, max = Number.NEGATIVE_INFINITY;
        array.forEach(x => {
            if (x < min) {
                min = x;
            }
            if (x > max) {
                max = x;
            }
        });

        return { 
            "type": "continuous",
            "min": min, 
            "max": max 
        };
    }
}

function guess_compression(x, compression) {
    if (compression !== null) {
        return compression;
    }

    let buffer;
    if (x instanceof Uint8Array) {
        buffer = x;
    } else {
        buffer = astream.peek(x, 3);
    }

    // Compare against magic words for auto-detection.
    if (buffer.length >= 3 && buffer[0] == 0x1F && buffer[1] == 0x8B && buffer[2] == 0x08) {
        return 'gz';
    }

    return 'none';
}

export function unpackText(buffer, { compression = null } = {}) {
    compression = guess_compression(buffer, compression);
    let txt = (compression === "gz" ? pako.ungzip(buffer) : buffer);
    const dec = new TextDecoder();
    return dec.decode(txt);
}

// Soft-deprecated as of 1.1.0.
export function readLines(buffer, { compression = null } = {}) {
    let decoded = unpackText(buffer, { compression: compression });
    let lines = decoded.split("\n");
    if (lines.length > 0 && lines[lines.length - 1] == "") { // ignoring the trailing newline.
        lines.pop();
    }
    return lines;    
}

function merge_bytes(leftovers, decoder) {
    let total = 0;
    for (const x of leftovers) {
        total += x.length;
    }

    let combined = new Uint8Array(total);
    total = 0;
    for (const x of leftovers) {
        combined.set(x, total);
        total += x.length;
    }

    return decoder.decode(combined);
}

async function stream_callback(x, compression, chunkSize, callback) {
    // Force the input to be either a Uint8Array or a file path string.
    if (typeof x == "string") {
        ;
    } else if (x instanceof Uint8Array) {
        ;
    } else if (x instanceof afile.SimpleFile) {
        x = x.content();
    } else {
        x = (new afile.SimpleFile(x, { name: "dummy" })).content();
    }

    if (guess_compression(x, compression) == "gz") {
        await (new Promise((resolve, reject) => {
            let gz = new pako.Inflate({ chunkSize: chunkSize });
            gz.onData = callback;
            gz.onEnd = status => {
                if (status) {
                    reject("gzip decompression failed; " + gz.msg);
                } else {
                    resolve(null);
                }
            };

            if (typeof x == "string") {
                astream.stream(x, chunkSize, chunk => gz.push(chunk), null, reject);
            } else {
                gz.push(x);
            }
        }));
        return;
    }

    // Remaining possibilities are uncompressed.
    if (typeof x == "string") {
        await (new Promise((resolve, reject) => astream.stream(x, chunkSize, callback, resolve, reject)));
        return;
    }

    callback(x);
    return;
}

/**
 * Read lines of text from a file, possibly with decompression.
 *
 * @param {string|Uint8Array|SimpleFile|File} x - Contents of the file to be read.
 * On Node.js, this may be a string containing a path to a file;
 * on browsers, this may be a File object.
 * @param {object} [options={}] - Optional parameters.
 * @param {?string} [options.compression=null] - Compression of `buffer`, either `"gz"` or `"none"`.
 * If `null`, it is determined automatically from the `buffer` header.
 * @param {number} [options.chunkSize=65536] - Chunk size in bytes to use for file reading (if `x` is a file path) and decompression (if `compression="gz"`).
 * Larger values improve speed at the cost of memory.
 *
 * @return {Array} Array of strings where each entry contains a line in `buffer`.
 * The newline itself is not included in each string.
 * @async 
 */
export async function readLines2(x, { compression = null, chunkSize = 65536 } = {}) {
    const dec = new TextDecoder;
    let leftovers = [];
    let lines = [];

    let callback = (chunk) => {
        let last = 0;
        for (var i = 0; i < chunk.length; i++) {
            if (chunk[i] == 10) { // i.e., ASCII newline.
                let current = chunk.subarray(last, i);
                if (leftovers.length) {
                    leftovers.push(current);
                    lines.push(merge_bytes(leftovers, dec));
                    leftovers = [];
                } else {
                    lines.push(dec.decode(current));
                }
                last = i + 1; // skip past the newline.
            }
        }

        if (last != chunk.length) {
            leftovers.push(chunk.slice(last)); // copy to avoid problems with ownership as chunk gets deref'd.
        }
    };

    await stream_callback(x, compression, chunkSize, callback);

    if (leftovers.length) {
        lines.push(merge_bytes(leftovers, dec));
    }

    return lines;    
}

// Soft-deprecated as of 1.1.0.
export function readTable(buffer, { compression = null, delim = "\t", firstOnly = false } = {}) {
    let decoded = unpackText(buffer, { compression: compression });
    let res = ppp.parse(decoded, { delimiter: delim, preview: (firstOnly ? 1 : 0) });

    // Handle terminating newlines.
    let last = res.data[res.data.length - 1];
    if (last.length === 1 && last[0] === "") {
        res.data.pop();
    }

    return res.data;
}

/**
 * Read a delimiter-separated table from a buffer, possibly with decompression.
 * This assumes that newlines represent the end of each row of the table, i.e., there cannot be newlines inside quoted strings.
 *
 * @param {string|Uint8Array|SimpleFile|File} x - Contents of the file to be read.
 * On Node.js, this may be a string containing a path to a file;
 * on browsers, this may be a File object.
 * @param {object} [options={}] - Optional parameters.
 * @param {?string} [options.compression=null] - Compression of `buffer`, either `"gz"` or `"none"`.
 * If `null`, it is determined automatically from the `buffer` header.
 * @param {string} [options.delim="\t"] - Delimiter between fields.
 * @param {number} [options.chunkSize=1048576] - Chunk size in bytes to use for file reading (if `x` is a path), parsing of rows, and decompression (if `compression="gz"`).
 * Larger values improve speed at the cost of memory.
 *
 * @return {Array} Array of length equal to the number of lines in `buffer`.
 * Each entry is an array of strings, containing the `delim`-separated fields for its corresponding line.
 *
 * @async
 */
export async function readTable2(x, { compression = null, delim = "\t", chunkSize = 1048576 } = {}) {
    const dec = new TextDecoder;

    let rows = [];
    let parse = (str) => {
        let out = ppp.parse(str, { delimiter: delim });
        if (out.meta.aborted) {
            let msg = "failed to parse delimited file";
            for (const e of out.errors) {
                msg += "; " + e.message;
            }
            throw new Error(msg);
        }
        for (const x of out.data) {
            rows.push(x);
        }
    };

    let leftovers = [];
    let size_left = 0;
    let callback = (chunk) => {
        let last = 0;
        for (var i = 0; i < chunk.length; i++) {
            // We assume that all newlines are end-of-rows, i.e., there are no
            // newlines inside quoted strings. Under this assumption, we can
            // safely chunk the input stream based on newlines, parse each
            // chunk, and then combine the parsing results together. To avoid
            // too many parsing calls, we accumulate buffers until we hit 
            // the chunkSize and then we decode + parse them altogether.
            if (chunk[i] == 10 && (i - last) + size_left >= chunkSize) {
                let current = chunk.subarray(last, i);
                if (leftovers.length) {
                    leftovers.push(current);
                    parse(merge_bytes(leftovers, dec));
                    leftovers = [];
                } else {
                    parse(dec.decode(current));
                }
                last = i + 1; // skip past the newline.
                size_left = 0;
            }
        }

        if (last != chunk.length) {
            leftovers.push(chunk.slice(last)); // copy to avoid problems with ownership as chunk gets deref'd.
            size_left += chunk.length - last;
        }
    };

    await stream_callback(x, compression, chunkSize, callback);

    if (leftovers.length) {
        let combined = merge_bytes(leftovers, dec);
        parse(combined);
        if (combined[combined.length - 1] == "\n") { // guaranteed to have non-zero length, by virtue of how 'leftovers' is filled.
            rows.pop();            
        }
    }

    return rows;    
}

/**
 * Detect if an array contains only stringified numbers and, if so, convert it into a TypedArray.
 * Conversion will still be performed for non-number strings corresponding to missing values or explicit not-a-number entries.
 *
 * @param {Array} x Array of strings, usually corresponding to a column in a table read by {@linkcode readDSVFromBuffer}.
 *
 * @return {?Float64Array} A Float64Array is returned if `x` contains stringified numbers.
 * Otherwise, `null` is returned if the conversion could not be performed.
 */
export function promoteToNumber(x) {
    let as_num = new Float64Array(x.length);

    for (const [i, v] of Object.entries(x)) {
        // See discussion at https://stackoverflow.com/questions/175739/how-can-i-check-if-a-string-is-a-valid-number.
        let opt1 = Number(v);
        let opt2 = parseFloat(v);
        if (!isNaN(opt1) && !isNaN(opt2)) {
            as_num[i] = opt1;
        } else if (v === "" || v === "NA" || v == "na" || v == "NaN" || v == "nan") {
            as_num[i] = NaN;
        } else if (v == "Inf" || v == "inf") {
            as_num[i] = Number.POSITIVE_INFINITY;
        } else if (v == "-Inf" || v == "-inf") {
            as_num[i] = Number.NEGATIVE_INFINITY;
        } else {
            return null;
        }
    }

    return as_num;
}
