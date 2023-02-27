import * as sce from "./SingleCellExperiment.js";
import * as adump from "./abstract/dump.js";
import JSZip from "jszip";

/**
 * Save the analysis results into an [**ArtifactDB** representation](https://github.com/ArtifactDB/BiocObjectSchemas) of a SingleCellExperiment.
 * This uses a language-agnostic format mostly based on HDF5 and JSON, which can be read into a variety of frameworks like R and Python.
 * The aim is to facilitate downstream analysis procedures that are not supported by **bakana** itself; 
 * for example, a bench scientist can do a first pass with **kana** before passing the results to a computational collaborator for deeper inspection.
 *
 * Note that all indices are 0-based and should be incremented by 1 before use in 1-based languages like R.
 * This involves the block assignments (indexing the block levels), the cluster assignments (indexing the marker results), and the custom selections (indexing the matrix columns).
 *
 * @param {object} state - Existing analysis state containing results, after one or more runs of {@linkcode runAnalysis}.
 * @param {string} name - Name of the SingleCellExperiment to be saved.
 * @param {object} [options={}] - Optional parameters.
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
 */
export async function saveSingleCellExperiment(state, name, { forceBuffer = null, directory = null } = {}) {
    if (directory !== null) {
        forceBuffer = false;
    } 

    let { metadata, files } = await sce.dumpSingleCellExperiment(state, name, { forceBuffer });
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
