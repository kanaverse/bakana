import * as adb from "./ArtifactDB.js";
import JSZip from "jszip";

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

/************************
 ******* Dataset ********
 ************************/

/**
 * Dataset as a ZIP file containing a SummarizedExperiment in the **ArtifactDB** representation,
 * e.g., as produced by {@linkcode saveSingleCellExperiment}.
 * Specifically, the ZIP file should contain the contents of an **ArtifactDB** project directory.
 * This project directory may contain multiple objects; the SummarizedExperiment of interest is identified in the constructor.
 *
 * @extends ArtifactDbSummarizedExperimentDatasetBase
 */
export class ZippedArtifactdbDataset extends adb.ArtifactDbSummarizedExperimentDatasetBase {
    #zipfile;
    #name;

    /**
     * @param {string} name - Name of the SummarizedExperiment object inside the project directory.
     * @param {SimpleFile} zipfile - A {@linkplain SimpleFile} object representing the ZIP file containing the project directory.
     * @param {object} [options={}] - Optional parameters, including those to be passed to the {@linkplain ArtifactdbSummarizedExperimentDatasetBase} constructor.
     * @param {?JSZip} [options.existingHandle=null] - An existing handle into the ZIP file, generated using the [**JSZip**](https://stuk.github.io/jszip/) package.
     * If an existing handle already exists, passing it in here will allow it to be re-used for greater efficiency.
     * If `null`, a new handle is created for this ZippedArtifactdbDataset instance.
     */
    constructor(name, zipfile, options={}) {
        let ziphandle = null;
        if ("existingHandle" in options) {
            ziphandle = options.existingHandle;
            delete options.existingHandle;
        }

        let nav = new ZippedProjectNavigator(zipfile, ziphandle);
        super(name, nav, options);
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
        return new ZippedArtifactdbDataset(name, files[0].file, options);
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
 * @extends ArtifactDbSummarizedExperimentResultBase
 */
export class ZippedArtifactdbResult extends adb.ArtifactDbSummarizedExperimentResultBase {
    /**
     * @param {string} name - Name of the SummarizedExperiment object inside the project directory.
     * @param {SimpleFile} zipfile - A {@linkplain SimpleFile} object representing the ZIP file containing the project directory.
     * @param {object} [options={}] - Optional parameters, including those to be passed to the {@linkplain ArtifactdbSummarizedExperimentDatasetBase} constructor.
     * @param {?JSZip} [options.existingHandle=null] - An existing handle into the ZIP file, generated using the [**JSZip**](https://stuk.github.io/jszip/) package.
     * If an existing handle already exists, passing it in here will allow it to be re-used for greater efficiency.
     * If `null`, a new handle is created for this ZippedArtifactdbDataset instance.
     */
    constructor(name, zipfile, options={}) {
        let ziphandle = null;
        if ("existingHandle" in options) {
            ziphandle = options.existingHandle;
            delete options.existingHandle;
        }

        let nav = new ZippedProjectNavigator(zipfile, ziphandle);
        super(name, nav, options);
    }
}
