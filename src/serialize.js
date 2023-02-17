import * as aserialize from "./abstract/serialize.js";
import * as readers from "./readers/index.js";
import * as anal from "./analysis.js";

/**
 * Format a collection of {@linkplain Dataset} objects so that they can be saved to file.
 *
 * @param {object} datasets - Object containing Dataset instances, just like that used in {@linkcode InputsState#compute InputsState.compute}.
 * @param {function} saver - Function that converts a {@linkplain SimpleFile} instance into an identifier string.
 * Specifically, it should accept two arguments:
 *
 * 1. A string containing the name of the Dataset.
 * 2. A string containing the format of the Dataset, e.g., `"10X"`, `"MatrixMarket"`.
 * 3. A SimpleFile object representing one of the files of that Dataset.
 *
 * It should then return a string that uniquely identifies this file within `datasets`.
 * The nature of this string is left to the application, e.g., it may be a file path for Node.js, a virtual file path in the browser, or some database identifier.
 * This function may be async.
 *
 * @return {object} Object containing information about the files and datasets in `datasets`.
 * @async
 */
export async function serializeDatasets(datasets, saver) {
    let output = {};

    for (const [key, val] of Object.entries(datasets)) {
        let dformat = val.constructor.format();
        let { files, options } = await val.serialize();

        let current = {
            format: dformat,
            options: options,
            files: []
        };

        for (const obj of files) {
            current.files.push({
                type: obj.type,
                name: obj.file.name(),
                id: await saver(key, dformat, obj.file)
            });
        }

        output[key] = current;
    }

    return output;
}

/**
 * Unserialize dataset information into their corresponding {@linkplain Dataset} instances.
 * This assumes that {@linkcode availableReaders} has been configured for all dataset formats that might be present.
 *
 * @param {object} serialized - Object containing the output of {@linkcode serializeDatasets}.
 * @param {function} loader - Function that accepts a single argument, the identifier string produced by `saver` in {@linkcode serializeDatasets};
 * and returns any value that can be used in the {@linkplain SimpleFile} constructor.
 * This may be async.
 *
 * @return {object} An object containing {@linkplain Dataset} instances that can be directly used in {@linkcode InputsState#compute InputsState.compute}.
 * @async
 */
export async function unserializeDatasets(serialized, loader) {
    let output = {};
    let known = readers.availableReaders;

    for (const [key, val] of Object.entries(serialized)) {
        if (!(val.format in known)) {
            throw new Error("unknown dataset format '" + val.format + "'");
        }
        let cls = readers.availableReaders[val.format];

        let handles = [];
        for (const obj of val.files) {
            let b = await loader(obj.id);
            let handle = new readers.SimpleFile(b, { name: obj.name }) 
            handles.push({ type: obj.type, file: handle });
        }

        output[key] = await cls.unserialize(handles, val.options);
    }

    return output;
}

/**
 * Save the analysis configuration to file, including the parameters and datasets.
 * This can be stringified and saved to file, or it can be used in {@linkcode unserializeConfiguration}.
 *
 * @param {object} state - State object produced by {@linkcode createAnalysis} and run through {@linkcode runAnalysis}.
 * @param {function} saver - Function to save files, see {@linkcode serializeDatasets} for more details.
 *
 * @return {object} Object containing the serialized analysis configuration, with the following properties:
 *
 * - `parameters`, an object containing parameters that can be used in {@linkcode runAnalysis}.
 * - `datasets`, an object containing serialized datasets that can be used in {@linkcode unserializeDatasets}.
 * - `other`, an object containing more parameters that need special handling outside of `parameters`.
 *   This typically involves calling setter functions directly on the State objects:
 *   - `inputs.direct_subset` contains a direct subset that can be used in {@linkcode InputsState#setDirectSubset InputsState.setDirectSubset} before calling {@linkcode runAnalysis}.
 *   - `custom_selections.selections` contains selections that can be used in {@linkcode CustomSelectionsState#addSelection CustomSelectionsState.addSelection} after {@linkcode runAnalysis}.
 *
 * @async
 */
export async function serializeConfiguration(state, saver) {
    let parameters = anal.retrieveParameters(state);
    let datasets = await serializeDatasets(state.inputs.fetchDatasets(), saver);

    let isub = state.inputs.fetchDirectSubset({ copy: false });
    if (isub !== null) {
        isub = Array.from(isub);
    }

    return {
        parameters: parameters,
        datasets: datasets,

        // Other parameters that need special handling.
        other: {
            inputs: {
                direct_subset: isub,
            },
            custom_selections: {
                selections: state.custom_selections.fetchSelections({ force: "Array" })
            }
        }
    };
}

/**
 * Load the analysis configuration from its serialized format.
 * This is effectively the reverse operation of {@linkcode serializeConfiguration}.
 *
 * @param {object} serialized - Configuration object produced by {@linkcode serializeConfiguration}.
 * @param {function} loader - Function to load files, see {@linkcode unserializeDatasets} for more details.
 * @param {object} [options={}] - Optional parameters.
 * @param {?function} [options.startFun=null] - Passed directly to {@linkcode runAnalysis}.
 * @param {?function} [options.finishFun=null] - Passed directly to {@linkcode runAnalysis}.
 *
 * @return {object} State object containing analysis results.
 * This is identical to the `state` passed into {@linkcode serializeConfiguration}.
 * @async
 */
export async function unserializeConfiguration(serialized, loader, { startFun = null, finishFun = null } = {}) {
    let state = await anal.createAnalysis();

    // Set this before running the analysis.
    if ("other" in serialized && "inputs" in serialized.other && "direct_subset" in serialized.other.inputs) {
        if (serialized.other.inputs.direct_subset !== null) {
            state.inputs.setDirectSubset(new Int32Array(serialized.other.inputs.direct_subset));
        }
    }

    let datasets = await unserializeDatasets(serialized.datasets, loader);
    await anal.runAnalysis(state, datasets, serialized.parameters, { startFun, finishFun });

    // Set this after the analysis is done, as the markers get computed directly.
    if ("other" in serialized && "custom_selections" in serialized.other && "selections" in serialized.other.custom_selections) {
        for (const [k, v] of Object.entries(serialized.other.custom_selections.selections)) {
            state.custom_selections.addSelection(k, new Int32Array(v), { copy: false });
        }
    }

    return state;
}

/****************************
 ********* LEGACY ***********
 ****************************/

/**
 * Parse a `*.kana` file by extracting the HDF5 state file and returning a function to extract embeddded data files.
 *
 * @param {string|Uint8Array} input - The input `*.kana` file.
 * In general, this should be a Uint8Array containing the full file contents.
 * For Node.js, this may also be a string containing a path to the file.
 * @param {string} statePath - String containing a file path to save the HDF5 state file.
 * This will also be the path supplied to {@linkcode loadAnalysis} to load the state into memory.
 * On browsers, this will exist inside the virtual file system of the **scran.js** module.
 * @param {object} [options] - Further options. 
 * For Node.js, callers can specify `stageDir`, a string containing a path to a staging directory for the extracted data files.
 *
 * @return {object|Promise<object>}
 * An object containing:
 *
 * - `version`: the version number of the `kana` format, as XXXYYYZZZ for "X.Y.Z".
 * - `embedded`: whether the data files are embedded in `input`. 
 * - `loader`: a function that extracts each data file given its offset and size.
 *   This should be used as `loadFun` in {@linkcode loadAnalysis}.
 *   If `embedded = false`, this is set to `null` instead.
 *
 * For Node.js, a promise is returned that evaluates to the above object.
 *
 * The HDF5 state file is also written to `statePath`.
 */
export function parseKanaFile(input, statePath, options = {}) {
    return aserialize.parseKanaFileInternal(input, statePath, options);
}
