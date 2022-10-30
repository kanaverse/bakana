import * as aserialize from "./abstract/serialize.js";

/**
 * Current version of the kana format.
 * This is encoded as an `XXXYYYZZZ` integer for version `XXX.YYY.ZZZ`, e.g., 1001002 for version 1.1.2.
 * @type {number}
 */
export const kanaFormatVersion = aserialize.FORMAT_VERSION;

/**
 * Create a `*.kana` file from the HDF5 state file and the various data files.
 *
 * @param {string} statePath - String containing a file path to an existing HDF5 state file.
 * This should be the same as the `path` used in {@linkcode saveAnalysis}.
 * On browsers, the path should exist inside the virtual file system of the **scran.js** module.
 * @param {Array} inputFiles - Array of files to be embedded into the `*.kana` file.
 * On Node.js, this should be an array of file paths; on browsers, this should be an array of ArrayBuffers.
 * Typically this is obtained as the resolved return value of {@linkcode saveAnalysis}.
 *
 * If `null`, it is assumed that files are linked instead of embedded.
 * @param {object} [options] - Further options. 
 * For Node.js, callers can specify `outputPath`, a string containing the output path for the newly created `*.kana` file.
 *
 * @return 
 * For Node.js, a promise is returned that resolves to a path to a new `*.kana` file. 
 * This is equal to `outputPath` if supplied, otherwise a path to a file in a temporary directory is returned.
 *
 * In browsers, an ArrayBuffer is returned containing the full contents of the new `*.kana` file.
 */
export function createKanaFile(statePath, inputFiles, options = {}) {
    return aserialize.createKanaFileInternal(statePath, inputFiles, options);
}

/**
 * Parse a `*.kana` file by extracting the HDF5 state file and returning a function to extract embeddded data files.
 *
 * @param {string|ArrayBuffer} input - The input `*.kana` file.
 * For Node.js, this should be a string containing a path to the file.
 * On browsers, this should be an ArrayBuffer containing the full file contents.
 * @param {string} statePath - String containing a file path to save the HDF5 state file.
 * This will also be the path supplied to {@linkcode loadAnalysis} to load the state into memory.
 * On browsers, this will exist inside the virtual file system of the **scran.js** module.
 * @param {object} [options] - Further options. 
 * For Node.js, callers can specify `stageDir`, a string containing a path to a staging directory for the extracted data files.
 *
 * @return The HDF5 state file is written to `statePath`.

 * In the browser, if `input` contains embedded files, a function is returned that extracts each data file given its offset and size.
 * (This should be used as `loadFun` in {@linkcode loadAnalysis}.)
 * If `input` contains linked files, `null` is returned.
 *
 * For Node.js, a promise is returned that evaluates to the aforementioned function or `null`. 
 */
export function parseKanaFile(input, statePath, options = {}) {
    return aserialize.parseKanaFileInternal(input, statePath, options);
}
