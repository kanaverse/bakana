import * as adb from "./ArtifactDB-abstract.js";
import JSZip from "jszip";
import * as afile from "./abstract/file.js";

class ZippedProjectNavigator {
    #zipfile;
    #ziphandle;

    constructor(zipfile, ziphandle) {
        this.#zipfile = zipfile;
        this.#ziphandle = null;
    }

    async file(path) {
        if (this.#ziphandle == null) {
            this.#ziphandle = await JSZip.loadAsync(this.#zipfile.buffer());
        }
        return await this.#ziphandle.file(path).async("uint8array");
    }

    async metadata(path) {
        if (this.#ziphandle == null) {
            this.#ziphandle = await JSZip.loadAsync(this.#zipfile.buffer());
        }

        while (1) {
            if (!path.endsWith(".json")) { 
                path += ".json";
            }

            let contents = await this.#ziphandle.file(path).async("string");
            let values = JSON.parse(contents);

            if (values["$schema"].startsWith("redirection/")){
                path = values.redirection.targets[0].location;
            } else {
                return values;
            }
        }
    }

    clear() {
        this.#ziphandle = null;
    }
};

/**
 * Search a ZIP file for SummarizedExperiments to use in {@linkplain ZippedArtifactdbDataset} or {@linkplain ZippedArtifactdbResult}.
 *
 * @param {JSZip} handle - A handle into the ZIP file, generated using the [**JSZip**](https://stuk.github.io/jszip/) package.
 * 
 * @return {Map} Object where the keys are the paths/names of possible SummarizedExperiment objects,
 * and each value is a 2-element array of dimensions.
 */
export async function searchZippedArtifactdb(handle) {
    // Sorting by the number of slashes.
    let all_json = [];
    for (const name of Object.keys(handle.files)) {
        if (name.endsWith(".json")) {
            all_json.push({ name: name, path: name.split("/") });
        }
    }
    all_json.sort((a, b) => a.path.length - b.path.length);

    let found_se = new Map;
    let nonchildren = new Map;
    let redirects = new Map;

    for (const x of all_json) {
        // Avoid loading JSONs for files in subdirectories of known SEs.
        let current = found_se;
        let already_found = false;

        for (const comp of x.path) {
            let val = current.get(comp);
            if (typeof val === "undefined") {
                val = new Map;
                current.set(comp, val);
            } else if (val === null) {
                already_found = true;
                break;
            }
            current = val;
        }

        // Otherwise, we load it in and peel out some information.
        if (!already_found) {
            let contents = await handle.file(x.name).async("string");
            let values = JSON.parse(contents);
            if ("summarized_experiment" in values) {
                nonchildren.set(values.path, values.summarized_experiment.dimensions);            
            } else if (values["$schema"].startsWith("redirection/")) {
                redirects.set(values.path, values.redirection.targets[0].location);
            }
        }
    }

    for (const [rr, loc] of redirects) {
        let found = nonchildren.get(loc);
        if (typeof found !== "undefined") {
            nonchildren.delete(loc);
            nonchildren.set(rr, found);
        }
    }

    return nonchildren;
}

/************************
 ******* Dataset ********
 ************************/

/**
 * Dataset as a ZIP file containing a SummarizedExperiment in the **ArtifactDB** representation,
 * e.g., as produced by {@linkcode saveSingleCellExperiment}.
 * Specifically, the ZIP file should contain the contents of an **ArtifactDB** project directory.
 * This project directory may contain multiple objects; the SummarizedExperiment of interest is identified in the constructor.
 *
 * @extends AbstractArtifactdbDataset
 */
export class ZippedArtifactdbDataset extends adb.AbstractArtifactdbDataset {
    #zipfile;
    #name;

    /**
     * @param {string} name - Name of the SummarizedExperiment object inside the project directory.
     * @param {SimpleFile|string|File} zipfile - Contents of the ZIP file containing the project directory.
     * On browsers, this may be a File object.
     * On Node.js, this may also be a string containing a file path.
     * @param {object} [options={}] - Optional parameters. 
     * @param {?JSZip} [options.existingHandle=null] - An existing handle into the ZIP file, generated using the [**JSZip**](https://stuk.github.io/jszip/) package.
     * If an existing handle already exists, passing it in here will allow it to be re-used for greater efficiency.
     * If `null`, a new handle is created for this ZippedArtifactdbDataset instance.
     */
    constructor(name, zipfile, options={}) {
        let ziphandle = null;
        if ("existingHandle" in options) {
            ziphandle = options.existingHandle;
            delete options.existingHandle;
        } else {
            if (!(zipfile instanceof afile.SimpleFile)) {
                zipfile = new afile.SimpleFile(zipfile);
            }
        }

        let nav = new ZippedProjectNavigator(zipfile, ziphandle);
        super(name, nav);
        this.#zipfile = zipfile;
        this.#name = name;
    }

    /**
     * @return {string} String specifying the format for this dataset.
     */
    static format() {
        return "ArtifactDB-zipped";
    }

    #dump_summary(fun) {
        let files = [ { type: "zip", file: fun(this.#zipfile) } ]; 
        let opt = this.options();
        opt.datasetName = this.#name; // storing the name as a special option... can't be bothered to store it as a separate file.
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
     * @param {Array} files - Array of objects like that produced by {@linkcode ZippedArtifactdbDataset#serialize serialize}.
     * @param {object} options - Object containing additional options to be passed to the constructor.
     * @return {ZippedArtifactdbDataset} A new instance of this class.
     * @static
     */
    static unserialize(files, options) {
        if (files.length != 1 || files[0].type != "zip") {
            throw new Error("expected exactly one file of type 'zip' for Zipped ArtifactDB unserialization");
        }

        let name = options.datasetName;
        delete options.datasetName;

        let output = new ZippedArtifactdbDataset(name, files[0].file);
        output.setOptions(output);
        return output;
    }
}

/***********************
 ******* Result ********
 ***********************/

/**
 * Result as a ZIP file containing a SummarizedExperiment in the **ArtifactDB** representation,
 * e.g., as produced by {@linkcode saveSingleCellExperiment}.
 * Specifically, the ZIP file should contain the contents of an **ArtifactDB** project directory.
 * This project directory may contain multiple objects; the SummarizedExperiment of interest is identified in the constructor.
 *
 * @extends AbstractArtifactdbResult
 */
export class ZippedArtifactdbResult extends adb.AbstractArtifactdbResult {
    /**
     * @param {string} name - Name of the SummarizedExperiment object inside the project directory.
     * @param {SimpleFile|string|File} zipfile - Contents of the ZIP file containing the project directory.
     * On browsers, this may be a File object.
     * On Node.js, this may also be a string containing a file path.
     * @param {object} [options={}] - Optional parameters. 
     * @param {?JSZip} [options.existingHandle=null] - An existing handle into the ZIP file, generated using the [**JSZip**](https://stuk.github.io/jszip/) package.
     * If an existing handle already exists, passing it in here will allow it to be re-used for greater efficiency.
     * If `null`, a new handle is created for this ZippedArtifactdbDataset instance.
     */
    constructor(name, zipfile, options={}) {
        let ziphandle = null;
        if ("existingHandle" in options) {
            ziphandle = options.existingHandle;
            delete options.existingHandle;
        } else {
            if (!(zipfile instanceof afile.SimpleFile)) {
                zipfile = new afile.SimpleFile(zipfile);
            }
        }

        let nav = new ZippedProjectNavigator(zipfile, ziphandle);
        super(name, nav);
    }
}
