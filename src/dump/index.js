import * as sce from "./SingleCellExperiment.js";
import * as markers from "./markers.js";
import * as adump from "./abstract/dump.js";
import JSZip from "jszip";

/**
 * Save the analysis results into an [**ArtifactDB** representation](https://github.com/ArtifactDB/BiocObjectSchemas) of a SingleCellExperiment.
 * This uses a language-agnostic format mostly based on HDF5 and JSON, which can be read into a variety of frameworks like R and Python.
 * The aim is to facilitate downstream analysis procedures that are not supported by **bakana** itself; 
 * for example, a bench scientist can do a first pass with **kana** before passing the results to a computational collaborator for deeper inspection.
 *
 * @param {object} state - Existing analysis state containing results, after one or more runs of {@linkcode runAnalysis}.
 * @param {string} name - Name of the SingleCellExperiment to be saved.
 * @param {object} [options={}] - Optional parameters.
 * @param {boolean} [options.reportOneIndex=false] - Whether to report 1-based indices, for use in 1-based languages like R.
 * Currently, this refers to the column indices of the custom selections reported in the SingleCellExperiment's metadata.
 * @param {boolean} [options.forceBuffer=true] - Whether to force all files to be loaded into memory as Uint8Array buffers.
 * @param {?directory} [options.directory=null] - Project directory in which to save the file components of the SingleCellExperiment.
 * Only used for Node.js; if supplied, it overrides any setting of `forceBuffer`.
 *
 * @return {Array} Array of objects where each object corresponds to a file in the SingleCellExperiment's representation.
 * Each inner object contains `metadata`, another object containing the metadata for the corresponding file.
 *
 * For files that are not purely metadata, the inner object will also contain  `contents`, the contents of the file.
 * The type of `contents` depends on the context and optional parameters.
 *
 * - If `forceBuffer = true`, `contents` is a Uint8Array of the file contents.
 *   This is probably the most advisable setting for browser environments.
 * - If `forceBuffer = false`, `contents` may be either a Uint8Array or a string to a temporary file path. 
 *   In a browser context, the temporary path is located on the **scran.js** virtual file system,
 *   see [here](https://kanaverse.github.io/scran.js/global.html#readFile) to extract its contents and to clean up afterwards.
 * - If the function is called from Node.js and `directory` is supplied, `contents` is a string to a file path inside `directory`.
 *   This overrides any expectations from the setting of `forceBuffer` and is the most efficient approach for Node.js.
 * 
 * The SingleCellExperiment itself will be accessible from a top-level object at the path defined by `name`.
 */
export async function saveSingleCellExperiment(state, name, { reportOneIndex = false, forceBuffer = null, directory = null } = {}) {
    if (directory !== null) {
        forceBuffer = false;
    } 

    let { metadata, files } = await sce.dumpSingleCellExperiment(state, name, { forceBuffer, reportOneIndex });
    files.push({ metadata });

    // Added a redirection document.
    files.push({
        "metadata": {
            "$schema": "redirection/v1.json",
            "path": name,
            "redirection": {
                "targets": [
                    {
                        "type": "local",
                        "location": metadata.path
                    }
                ]
            }
        }
    });

    await adump.attachMd5sums(files);

    // Either dumping everything to file or returning all the buffers.
    if (!forceBuffer && directory !== null) {
        await adump.realizeDirectory(files, directory, name);
    }

    return files;
}

/**
 * Save the per-gene analysis results into [**ArtifactDB** data frames](https://github.com/ArtifactDB/BiocObjectSchemas).
 * This includes the marker tables for the clusters and custom selections, as well as the variance modelling statistics from the feature selection step.
 * Each data frames has the same order of rows as the SingleCellExperiment saved by {@linkcode saveSingleCellExperiment};
 * for easier downstream use, we set the row names of each data frame to row names of the SingleCellExperiment
 * (or if no row names are available, we set each data frame's row names to the first column of the `rowData`).
 *
 * @param {object} state - Existing analysis state containing results, after one or more runs of {@linkcode runAnalysis}.
 * @param {string} prefix - Prefix to attach to the path for each result.
 * This can contain subdirectories, or it can be an empty string.
 * @param {object} [options={}] - Optional parameters.
 * @param {boolean} [options.includeMarkerDetection=true] - Whether to save the marker detection results.
 * @param {boolean} [options.includeCustomSelections=true] - Whether to save the custom selection results.
 * @param {boolean} [options.includeFeatureSelection=true] - Whether to save the feature selection results.
 * @param {boolean} [options.forceBuffer=true] - Whether to force all files to be loaded into memory as Uint8Array buffers.
 * @param {?directory} [options.directory=null] - Project directory in which to save the file components of the SingleCellExperiment.
 * Only used for Node.js; if supplied, it overrides any setting of `forceBuffer`.
 *
 * @return {Array} Array of objects where each object corresponds to a data frame of per-gene statistics.
 * These data frames are grouped by analysis step, i.e., `marker_detection`, `custom_selections` and `feature_selection`.
 * See {@linkcode saveSingleCellExperiment} for more details on the contents of each object inside the returned array.
 */
export async function saveGenewiseResults(state, name, { includeMarkerDetection = true, includeCustomSelections = true, includeFeatureSelection = true, forceBuffer = null, directory = null } = {}) {
    let modalities = {};
    let anno = state.inputs.fetchFeatureAnnotations();
    for (const [k, v] of Object.entries(anno)) {
        let rn = v.rowNames();
        if (rn === null && v.numberOfColumns() > 0) {
            rn = v.column(0);
        }
        modalities[k] = rn;
    }

    // Stripping trailing slashes.
    if (name.endsWith("/")) {
        name = name.replace(/\/+$/, "");
    }

    let files = [];

    if (includeMarkerDetection) {
        markers.dumpMarkerDetectionResults(state, modalities, name, files, forceBuffer);
    }

    if (includeCustomSelections) {
        markers.dumpCustomSelectionResults(state, modalities, name, files, forceBuffer);
    }

    if (state.feature_selection.valid() && includeFeatureSelection) {
        markers.dumpFeatureSelectionResults(state, modalities.RNA, name, files, forceBuffer);
    }

    await adump.attachMd5sums(files);

    // Either dumping everything to file or returning all the buffers.
    if (!forceBuffer && directory !== null) {
        await adump.realizeDirectory(files, directory, name);
    }

    return files;
}

/**
 * Zip a set of files, typically corresponding to the contents of an **ArtifactDB** project directory.
 *
 * @param {Array} files - Array of objects describing files in the project directory.
 * See the output of {@linkcode saveSingleCellExperiment} for more details on the expected format.
 *
 * @return {Promise<Uint8Array>} Uint8Array with the contents of the ZIP file containing all `files`.
 */
export function zipFiles(files) {
    var zip = new JSZip();

    for (const x of files) {
        let suffix;
        if ("contents" in x) {
            suffix = ".json";
            if (x.contents instanceof Uint8Array) {
                zip.file(x.metadata.path, x.contents);
            } else {
                zip.file(x.metadata.path, adump.loadFilePath(x.contents));
            }
        } else if (x.metadata["$schema"].startsWith("redirection/")) {
            suffix = ".json";
        } else {
            suffix = "";
        }

        // Add a trailing newline to avoid no-newline warnings. 
        zip.file(x.metadata.path + suffix, JSON.stringify(x.metadata, null, 2) + "\n");
    }

    return zip.generateAsync({ type: "uint8array" });
}
