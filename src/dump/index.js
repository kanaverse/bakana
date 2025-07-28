import * as sce from "./SingleCellExperiment.js";
import * as markers from "./markers.js";
import * as ass from "./assays.js";
import * as rd from "./reducedDimensions.js";
import * as internal from "./abstract/dump.js";

/**
 * Save the analysis results into a SingleCellExperiment using the [**takane**](https://github.com/ArtifactDB/takane) format.
 * This uses a language-agnostic format mostly based on HDF5 and JSON, which can be read into a variety of frameworks like R and Python.
 * The aim is to facilitate downstream analysis procedures that are not supported by **bakana** itself; 
 * for example, a bench scientist can do a first pass with **kana** before passing the results to a computational collaborator for deeper inspection.
 *
 * @param {object} state - Existing analysis state containing results, after one or more runs of {@linkcode runAnalysis}.
 * @param {string} name - Name of the SingleCellExperiment.
 * This may also contain forward slashes, in which case it is treated as a local path.
 * If a local filesystem is present and `directory` is supplied, `name` defines the subdirectory within `directory` in which the SingleCellExperiment is to be saved.
 * @param {object} [options={}] - Optional parameters.
 * @param {boolean} [options.reportOneIndex=false] - Whether to report 1-based indices, for use in 1-based languages like R.
 * Currently, this only refers to the column indices of the custom selections reported in the SingleCellExperiment's metadata.
 * @param {boolean} [options.storeModalityColumnData=false] - Whether to store modality-specific per-cell statistics (e.g., QC metrics, size factors) in the column data of the associated alternative Experiment.
 * This can yield a cleaner SingleCellExperiment as statistics are assigned to the relevant modalities.
 * That said, the default is still `false` as many applications (including **bakana** itself, via the {@linkcode AbstractArtifactdbDataset} and friends) will not parse the column data of alternative Experiments. 
 * In such cases, all modality-specific metrics are stored in the main experiment's column data with a modality-specific name, e.g., `kana::ADT::quality_control::sums`.
 * @param {?string} [options.directory=null] - Directory in which to save the components of the SingleCellExperiment.
 * If `null` or if no local file system exists, files are stored in memory as Uint8Arrays.
 *
 * @return {?Object} If `directory` is supplied and a local filesystem exists, the components are saved to disk and `null` is returned.
 * Otherwise, if `directory = null` or a local filesystem does not exist, an object is returned where each key is a local path to a component file and each value is a Uint8Array with the file contents.
 */
export async function saveSingleCellExperiment(state, name, { reportOneIndex = false, storeModalityColumnData = false, directory = null } = {}) {
    let my_sce = await sce.formatSingleCellExperiment(state, { forceBuffer, reportOneIndex, storeModalityColumnData });

    if (!(internal.fsexists())) {
        directory = null;
    }
    let files = {};
    let globals = new AlabasterGlobalsInterface(directory, files);

    jsp.saveObjectRegistry.push([ ass.MockSparseMatrix, ass.saveSparseMatrix ]);
    jsp.saveObjectRegistry.push([ ass.MockNormalizedMatrix, ass.saveNormalizedMatrix ]); 
    jsp.saveObjectRegistry.push([ rd.MockReducedDimensionMatrix, rd.saveReducedDimensionMatrix]); 
    try {
        await jsp.saveObject(my_sce, name, globals);
    } finally {
        jsp.saveObjectRegistry.pop();
        jsp.saveObjectRegistry.pop();
        jsp.saveObjectRegistry.pop();
    }

    if (directory === null) {
        return files;
    } else {
        return null;
    }
}

/**
 * Save the per-gene analysis results as data frames in the [**takane**](https://github.com/ArtifactDB/takane) format.
 * This includes the marker tables for the clusters and custom selections, as well as the variance modelling statistics from the feature selection step.
 * Each data frames has the same order of rows as the SingleCellExperiment saved by {@linkcode saveSingleCellExperiment};
 * for easier downstream use, we set the row names of each data frame to row names of the SingleCellExperiment
 * (or if no row names are available, we set each data frame's row names to the first column of the `rowData`).
 *
 * @param {object} state - Existing analysis state containing results, after one or more runs of {@linkcode runAnalysis}.
 * @param {?string} path - Local path in which to save all results.
 * If a local filesystem is present and `directory` is supplied, `path` defines the subdirectory within `directory` in which the SingleCellExperiment is to be saved.
 * If `null`, results are directly saved to `directory` itself.
 * @param {object} [options={}] - Optional parameters.
 * @param {boolean} [options.includeMarkerDetection=true] - Whether to save the marker detection results.
 * @param {boolean} [options.includeCustomSelections=true] - Whether to save the custom selection results.
 * @param {boolean} [options.includeFeatureSelection=true] - Whether to save the feature selection results.
 * @param {boolean} [options.forceBuffer=true] - Whether to force all files to be loaded into memory as Uint8Array buffers.
 * @param {?string} [options.directory=null] - Directory in which to save the data frames. 
 * If `null` or if no local file system exists, files are stored in memory as Uint8Arrays.
 *
 * @return {?Object} If `directory` is supplied and a local filesystem exists, the data frame components are saved to disk and `null` is returned.
 * Otherwise, if `directory = null` or a local filesystem does not exist, an object is returned where each key is a local path to a component file and each value is a Uint8Array with the file contents.
 */
export async function saveGenewiseResults(state, path, { includeMarkerDetection = true, includeCustomSelections = true, includeFeatureSelection = true, directory = null } = {}) {
    let modalities = {};
    let anno = state.inputs.fetchFeatureAnnotations();
    for (const [k, v] of Object.entries(anno)) {
        let rn = v.rowNames();
        if (rn === null && v.numberOfColumns() > 0) {
            rn = v.column(0);
        }
        modalities[k] = rn;
    }

    if (!(internal.fsexists())) {
        directory = null;
    }
    let files = {};
    let globals = new AlabasterGlobalsInterface(directory, forceBuffer, files);

    if (includeMarkerDetection) {
        let dir = "marker_detection";
        if (path !== null) {
            dir = jsp.joinPath(dir, path);
        }
        let all = markers.formatMarkerDetectionResults(state, modalities);
        for (const [k, v] of Object.entries(all)) {
            await jsp.saveObject(v, jsp.joinPath(dir, k), globals);
        }
    }

    if (includeCustomSelections) {
        let dir = "custom_selections";
        if (path !== null) {
            dir = jsp.joinPath(dir, path);
        }
        let all = markers.formatCustomSelectionResults(state, modalities);
        for (const [k, v] of Object.entries(all)) {
            await jsp.saveObject(v, jsp.joinPath(dir, k), globals);
        }
    }

    if (state.feature_selection.valid() && includeFeatureSelection) {
        let df = markers.dumpFeatureSelectionResults(state, modalities.RNA);
        await jsp.saveObject(df, jsp.joinPath(path, "feature_selection"), globals);
    }

    if (directory === null) {
        return files;
    } else {
        return null;
    }
}
