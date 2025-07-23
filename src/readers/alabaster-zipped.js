import JSZip from "jszip";
import * as jsp from "jaspagate";
import * as adb from "./alabaster-abstract.js";
import * as afile from "./abstract/file.js";

class ZippedProjectNavigator {
    #zipfile;
    #ziphandle;
    #prefix;

    constructor(zipfile, ziphandle, prefix) {
        this.#zipfile = zipfile;
        this.#ziphandle = null;
        this.#prefix = prefix;
    }

    async get(path, asBuffer) {
        if (this.#ziphandle == null) {
            this.#ziphandle = await JSZip.loadAsync(this.#zipfile.buffer());
        }
        // We always return a buffer.
        return this.#ziphandle.file(jsp.joinPath(this.#prefix, path)).async("uint8array");
    }

    exists(path) {
        return this.#ziphandle.file(jsp.joinPath(this.#prefix, path)) !== null; 
    }

    clean(path) {}
};

/**
 * Search a ZIP file for SummarizedExperiments to use in {@linkplain ZippedAlabasterDataset} or {@linkplain ZippedAlabasterResult}.
 *
 * @param {JSZip} handle - A handle into the ZIP file, generated using the [**JSZip**](https://stuk.github.io/jszip/) package.
 * 
 * @return {Map} Object where the keys are the paths/names of possible SummarizedExperiment objects,
 * and each value is a 2-element array of dimensions.
 */
export async function searchZippedAlabaster(handle) {
    let candidates = [];
    for (const name of Object.keys(handle.files)) {
        if (name == "OBJECT") {
            candidates.push("");
        } else if (name.endsWith("/OBJECT")) {
            candidates.push(name.slice(0, name.length - 6));
        }
    }
    candidates.sort();

    const nonchildren = new Map;
    let counter = 0;
    while (counter < candidates.length) {
        let prefix = candidates[counter];
        counter++;
        let contents = await handle.file(prefix).async("string");

        let meta = {};
        try {
            meta = JSON.parse(contents);
        } catch(e) {
            ;
        }

        if ("summarized_experiment" in meta) {
            nonchildren.set(prefix, meta.summarized_experiment.dimensions);            
            while (counter < candidates) {
                if (candidates[counter].startsWith(prefix)) {
                    break;
                }
                counter++;
            }
        }
    }

    return nonchildren;
}

/************************
 ******* Dataset ********
 ************************/

/**
 * Dataset as a ZIP file containing a SummarizedExperiment in the **alabaster** representation.
 * Specifically, the ZIP file should contain the contents of an **alabaster** project directory.
 * This project directory may contain multiple objects; the SummarizedExperiment of interest is identified in the constructor.
 *
 * @extends AbstractAlabasterDataset
 */
export class ZippedAlabasterDataset extends adb.AbstractAlabasterDataset {
    #zipfile;
    #prefix;

    /**
     * @param {string} prefix - Name of the SummarizedExperiment object inside the ZIP file.
     * This should be `.` if the `OBJECT` file is stored at the root of the ZIP file,
     * otherwise it should be the relative path to the object directory containing the `OBJECT` file.
     * @param {SimpleFile|string|File} zipfile - Contents of the ZIP file containing the project directory.
     * On browsers, this may be a File object.
     * On Node.js, this may also be a string containing a file path.
     * @param {object} [options={}] - Optional parameters. 
     * @param {?JSZip} [options.existingHandle=null] - An existing handle into the ZIP file, generated using the [**JSZip**](https://stuk.github.io/jszip/) package.
     * If an existing handle already exists, passing it in here will allow it to be re-used for greater efficiency.
     * If `null`, a new handle is created for this ZippedAlabasterDataset instance.
     */
    constructor(prefix, zipfile, options={}) {
        let ziphandle = null;
        if ("existingHandle" in options) {
            ziphandle = options.existingHandle;
            delete options.existingHandle;
        } else {
            if (!(zipfile instanceof afile.SimpleFile)) {
                zipfile = new afile.SimpleFile(zipfile);
            }
        }

        let nav = new ZippedProjectNavigator(zipfile, ziphandle, prefix);
        super(nav);
        this.#zipfile = zipfile;
        this.#prefix = prefix;
    }

    /**
     * @return {string} String specifying the format for this dataset.
     */
    static format() {
        return "alabaster-zipped";
    }

    #dump_summary(fun) {
        let files = [ { type: "zip", file: fun(this.#zipfile) } ]; 
        let opt = this.options();
        opt.datasetPrefix = this.#prefix; // storing the name as a special option... can't be bothered to store it as a separate file.
        return { files: files, options: opt };
    }

    /**
     * @return {object} Object containing the abbreviated details of this dataset.
     */
    abbreviate() {
        return this.#dump_summary(f => { 
            return { size: f.size(), name: f.name() }
        });
    }

    /**
     * @return {object} Object describing this dataset, containing:
     *
     * - `files`: Array of objects representing the files used in this dataset.
     *   Each object corresponds to a single file and contains:
     *   - `type`: a string denoting the type.
     *   - `file`: a {@linkplain SimpleFile} object representing the file contents.
     * - `options`: An object containing additional options to saved.
     */
    serialize() {
        return this.#dump_summary(f => f);
    }

    /**
     * @param {Array} files - Array of objects like that produced by {@linkcode ZippedAlabasterDataset#serialize serialize}.
     * @param {object} options - Object containing additional options to be passed to the constructor.
     * @return {ZippedAlabasterDataset} A new instance of this class.
     * @static
     */
    static unserialize(files, options) {
        if (files.length != 1 || files[0].type != "zip") {
            throw new Error("expected exactly one file of type 'zip' for Zipped alabaster unserialization");
        }

        let prefix = options.datasetPrefix;
        delete options.datasetPrefix;

        let output = new ZippedAlabasterDataset(prefix, files[0].file);
        output.setOptions(output);
        return output;
    }
}

/***********************
 ******* Result ********
 ***********************/

/**
 * Result as a ZIP file containing a SummarizedExperiment in the **alabaster** representation,
 * Specifically, the ZIP file should contain the contents of an **alabaster** project directory.
 * This project directory may contain multiple objects; the SummarizedExperiment of interest is identified in the constructor.
 *
 * @extends AbstractAlabasterResult
 */
export class ZippedAlabasterResult extends adb.AbstractAlabasterResult {
    /**
     * @param {string} prefix - Name of the SummarizedExperiment object inside the project directory.
     * @param {SimpleFile|string|File} zipfile - Contents of the ZIP file containing the project directory.
     * On browsers, this may be a File object.
     * On Node.js, this may also be a string containing a file path.
     * @param {object} [options={}] - Optional parameters. 
     * @param {?JSZip} [options.existingHandle=null] - An existing handle into the ZIP file, generated using the [**JSZip**](https://stuk.github.io/jszip/) package.
     * If an existing handle already exists, passing it in here will allow it to be re-used for greater efficiency.
     * If `null`, a new handle is created for this ZippedAlabasterDataset instance.
     */
    constructor(prefix, zipfile, options={}) {
        let ziphandle = null;
        if ("existingHandle" in options) {
            ziphandle = options.existingHandle;
            delete options.existingHandle;
        } else {
            if (!(zipfile instanceof afile.SimpleFile)) {
                zipfile = new afile.SimpleFile(zipfile);
            }
        }

        let nav = new ZippedProjectNavigator(zipfile, ziphandle, prefix);
        super(nav);
    }
}
