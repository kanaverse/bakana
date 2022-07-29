import * as pako from "pako";
import ppp from "papaparse";

export function extractHDF5Strings(handle, name) {
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
 * Summarize annotation columns for a preflight request.
 *
 * @param {Array|TypedArray} values - Column of annotation values, of length equal to the number of cells for a given matrix.
 * Array are treated as categorical data and should typically contain strings;
 * TypedArrays are treated as continuous data.
 * @param {number} limit - Maximum number of unique values to report for Array `values`.
 *
 * @return {object} Object containing `type`, a string indicating whether `values` was categorical or continuous.
 *
 * If `"categorical"`, the object will contain `values`, an array of unique values up to the length specified by `limit`.
 * It will also contain `truncated`, a boolean indicating whether the actual number of unique values exceeds `limit`.
 *
 * If `"continuous"`, the object will contain `min` and `max`.
 * These are numbers specifying the minimum and maximum value in `values`, respectively.
 */
export function summarizeValues(values, limit) {
    if (values instanceof Array) {
        let chosen = Array.from(new Set(values));
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
        let min = Infinity, max = -Infinity;
        values.forEach(x => {
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

/**
 * Unpack a buffer to text, possibly with decompression.
 *
 * @param {Uint8Array} buffer - Contents to be unpacked.
 * @param {object} [options] - Optional parameters.
 * @param {?string} [options.compression] - Compression mode used for `buffer`.
 * This can be `"none"` or `"gz"`, with `null` performing automatic detection based on the magic words.
 *
 * @return {string} Text encoded by `buffer`.
 */
export function unpackText(buffer, { compression = null } = {}) {
    let txt = buffer;
    
    // Compare against magic words for auto-detection.
    if (compression === null) {
        if (buffer.length >= 3 && buffer[0] == 0x1F && buffer[1] == 0x8B && buffer[2] == 0x08) {
            compression = 'gz';
        }
    }

    if (compression === "gz") {
        txt = pako.ungzip(buffer);
    }
    
    const dec = new TextDecoder();
    return dec.decode(txt);
}

/**
 * Read lines of text from a buffer, possibly with decompression.
 *
 * @param {Uint8Array} buffer - Content to be read.
 * @param {object} [options] - Optional parameters.
 * @param {?string} [options.compression] - See {@linkcode unpackText} for details.
 *
 * @return {Array} Array of strings where each entry contains a line in `buffer`.
 * The newline itself is not included in each string. 
 * A trailing newline at the end of each file is removed.
 */
export function readLines(buffer, { compression = null } = {}) {
    let decoded = unpackText(buffer, { compression: compression });
    let lines = decoded.split("\n");
    if (lines.length > 0 && lines[lines.length - 1] == "") { // ignoring the trailing newline.
        lines.pop();
    }
    return lines;    
}

/**
 * Read a delimiter-separated table from a buffer, possibly with decompression.
 *
 * @param {Uint8Array} buffer - Content to be read.
 * @param {object} [options] - Optional parameters.
 * @param {?string} [options.compression] - See {@linkcode unpackText} for details.
 * @param {string} [options.delim] - Delimiter between fields.
 * @param {boolean} [options.firstOnly] - Whether to only read the first line, presumably containing the header.
 *
 * @return {Array} Array of arrays where each entry corresponds to a line in `buffer` and contains the `delim`-separated fields.
 * If `firstOnly = true`, the output array only contains one element corresponding to the contents of the first line.
 */
export function readTable(content, { compression = null, delim = "\t", firstOnly = false } = {}) {
    let decoded = unpackText(content, { compression: compression });
    let res = ppp.parse(decoded, { delimiter: delim, preview: (firstOnly ? 1 : 0) });

    // Handle terminating newlines.
    let last = res.data[res.data.length - 1];
    if (last.length === 1 && last[0] === "") {
        res.data.pop();
    }

    return res.data;
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
