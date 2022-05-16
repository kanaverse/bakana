import * as pako from "pako";
import * as afile from "../abstract/file.js";
import * as scran from "scran.js";
import * as utils from "./general.js";
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
    let must_string = false;

    for (const [i, v] of Object.entries(x)) {
        // See discussion at https://stackoverflow.com/questions/175739/how-can-i-check-if-a-string-is-a-valid-number.
        let opt1 = Number(v);
        let opt2 = parseFloat(v);
        if (!isNaN(opt1) && !isNaN(opt2)) {
            as_num[i] = opt1;
        } else if (v === "" || v === "NA" || v == "na" || v == "NaN" || v == "nan") {
            as_num[i] = NaN;
        } else {
            return null;
        }
    }

    return as_num;
}

export class MultiMatrix {
    #store;
    #nrows;
    #ncols;

    constructor({ store = {} } = {}) {
        this.#store = store;
        this.#ncols = null;
        this.#nrows = 0;

        let keys = Object.keys(store);
        if (keys.length) {
            // We ignore numberOfColumns here, as everyone should have the same number of cells.
            for (var k = 0; k < keys.length; k++) {
                let current = store[keys[k]];
                if (k == 0) {
                    this.#ncols = current.numberOfColumns();
                } else if (current.numberOfColumns() != this.#ncols) {
                    throw new Error("all matrices should have the same number of columns");
                }
            }
            keys.forEach(x => { this.#nrows += store[x].numberOfRows(); });
        }
    }

    numberOfColumns() {
        return this.#ncols;
    }

    numberOfRows() {
        return this.#nrows;
    }

    available() {
        return Object.keys(this.#store);
    }

    has(i) {
        return (i in this.#store);
    }

    get(i) {
        return this.#store[i];
    }

    add(i, matrix) {
        if (this.#ncols === null) {
            this.#ncols = matrix.numberOfColumns();
        } else if (matrix.numberOfColumns() != this.#ncols) {
            throw new Error("all matrices should have the same number of columns");
        }

        if (i in this.#store) {
            let old = this.#store[i];
            this.#nrows -= old.numberOfRows();
            utils.freeCache(old);
        }

        this.#store[i] = matrix;
        this.#nrows += matrix.numberOfRows();
    }

    remove(i) {
        utils.freeCache(this.#store[i]);
        delete this.#store[i];
    }

    rename(from, to) {
        if (from !== to) {
            this.#store[to] = this.#store[from];
            delete this.#store[from];
        }
    }

    free() {
        for (const [x, v] of Object.entries(this.#store)) {
            utils.freeCache(v);
        }
        return;
    }
}

export function reorganizeGenes(matrix, geneInfo) {
    if (geneInfo === null) {
        let genes = [];
        if (matrix.isReorganized()) {
            let ids = matrix.identities();
            for (const i of ids) {
                genes.push(`Gene ${i + 1}`);
            }
        } else {
            for (let i = 0; i < matrix.numberOfRows(); i++) {
                genes.push(`Gene ${i + 1}`);
            }
        }
        geneInfo = { "id": genes };
    } else {
        if (matrix.isReorganized()) {
            scran.matchFeatureAnnotationToRowIdentities(matrix, geneInfo);
        }
    }
    return geneInfo;
}

var cache = {
    file2link: null, 
    link2file: null
};

/**
 * Specify a function to create links for data files.
 * By default, this only affects files for the MatrixMarket, H5AD and 10X formats, and is only used when linking is requested.
 *
 * @param {function} fun - Function that returns a linking idenfier to a data file.
 * The function should accept the following arguments:
 *
 * - A string specifying the type of the file, e.g., `"mtx"`, `"h5"`.
 * - A string containing the name of the file.
 * - An ArrayBuffer containing the file content (for browsers) or a string containing the file path (for Node.js),
 *
 * The function is expected to return a string containing some unique identifier to the file.
 * This is most typically used to register the file with some user-specified database system for later retrieval.
 *
 * @return `fun` is set as the global link creator for this step. 
 * The _previous_ value of the creator is returned.
 */
export function setCreateLink(fun) {
    let previous = cache.file2link;
    cache.file2link = fun;
    return previous;
}

/**
 * Specify a function to resolve links for data files.
 * By default, this only affects files for the MatrixMarket, H5AD and 10X formats, and is only used when links are detected.
 *
 * @param {function} fun - Function that accepts a string containing a linking idenfier and returns an ArrayBuffer containing the file content (for browsers) or a string containing the file path (for Node.js),
 * This is most typically used to retrieve a file from some user-specified database system.
 *
 * @return `fun` is set as the global resolver for this step. 
 * The _previous_ value of the resolver is returned.
 */
export function setResolveLink(fun) {
    let previous = cache.link2file;
    cache.link2file = fun;
    return previous;
}

export async function standardSerialize(details, type, embeddedSaver) {
    let output = { 
        "type": type,
        "name": details.name 
    };
    let serialized = details.content.serialized();

    if (embeddedSaver !== null) {
        let eout = await embeddedSaver(serialized, details.content.size());
        output.offset = eout.offset;
        output.size = eout.size;
    } else {
        let fun = cache.file2link;
        if (fun === null) {
            throw new Error("link-creating function has not been set by 'setCreateLink'");
        }
        output.id = await fun(type, details.name, serialized);
    }

    return output;
}

export async function standardUnserialize(details, embeddedLoader) {
    let output = { name: details.name };

    if ("id" in details) {
        let fun = cache.link2file;
        if (fun === null) {
            throw new Error("link-resolving function has not been set by 'setResolveLink'");
        }
        output.content = new afile.LoadedFile(await fun(details.id));
    } else {
        output.content = new afile.LoadedFile(await embeddedLoader(details.offset, details.size));
    }

    return output;
}

export function formatFile(file, sizeOnly) {
    let output = { "name": afile.rawName(file) };
    if (sizeOnly) {
        output.size = afile.rawSize(file);
    } else {
        output.content = new afile.LoadedFile(file);
    }
    return output;
}

export function splitByFeatureType(matrix, genes) { 
    if (!("type" in genes)) {
        return null;
    }

    let types = {};
    genes.type.forEach((x, i) => {
        if (!(x in types)) {
            types[x] = [];
        }
        types[x].push(i);
    });

    let out_mats = new MultiMatrix;
    let out_genes = {};
    try {
        for (const [k, v] of Object.entries(types)) {
            // Standardizing the names to something the rest of the pipeline
            // will recognize. By default, we check the 10X vocabulary here.
            let name = k;
            if (name.match(/gene expression/i)) {
                name = "RNA";
            } else if (name.match(/antibody capture/i)) {
                name = "ADT";
            }

            out_mats.add(name, scran.subsetRows(matrix, v));

            let curgenes = {};
            for (const [k2, v2] of Object.entries(output.genes)) {
                // Skipping 'type', as it's done its purpose now.
                if (k !== "type") {
                    curgenes[k2] = scran.matchVectorToRowIdentities(current.matrix, v2);
                }
            }
            out_genes[name] = curgenes;
        }

    } catch (e) {
        utils.freeCache(out_mats);
        throw e;
    }

    return { 
        matrices: out_mats, 
        genes: out_genes 
    };
}
