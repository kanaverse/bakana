import * as sce from "./SingleCellExperiment.js";
import * as adump from "./abstract/dump.js";

/**
 * Save the results into a [language-agnostic representation](https://github.com/ArtifactDB/BiocObjectSchemas) of a SingleCellExperiment.
 * Note that all row/column indices are 0-based and should be incremented by 1 before use in 1-based languages like R.
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
 *   see [here](https://jkanche.com/scran.js/global.html#readFile) to extract its contents and to clean up afterwards.
 * - If the function is called from Node.js and `directory` is supplied, `contents` is a string to a file path inside `directory.
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

    // Either dumping everything to file or returning all the buffers.
    if (!forceBuffer && directory !== null) {
        adump.realizeDirectory(files, directory, name);
    }

    return files;
}
