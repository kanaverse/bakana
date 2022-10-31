import * as scran from "scran.js";
import * as bioc from "bioconductor";
import * as afile from "./abstract/file.js";
import * as eutils from "./utils/extract.js";
import * as futils from "./utils/features.js";

/**
 * Dataset in the 10X HDF5 feature-barcode matrix format, see [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices) for details.
 */
export class TenxHdf5Dataset { 
    #h5_file;
    #h5_path;
    #h5_flush;

    #raw_features;
    #raw_cells;

    /**
     * @param {SimpleFile|string|Uint8Array|File} h5File - Contents of a HDF5 file.
     * On browsers, this may be a File object.
     * On Node.js, this may also be a string containing a file path.
     */
    constructor(h5File) {
        if (h5File instanceof afile.SimpleFile) {
            this.#h5_file = h5File;
        } else {
            this.#h5_file = new afile.SimpleFile(h5File);
        }
        this.clear();
    }

    #instantiate() {
        if (this.#h5_path !== null) {
            return;
        }

        let info = scran.realizeFile(this.#h5_file.content());
        this.#h5_path = info.path;
        this.#h5_flush = info.flush;
    }

    /**
     * Destroy caches if present, releasing the associated memory.
     * This may be called at any time but only has an effect if `cache = true` in {@linkcode TenxHdf5Dataset#load load} or {@linkcodeTenxHdf5Dataset#annotations annotations}. 
     */
    clear() {
        if (typeof this.#h5_flush == "function") {
            this.#h5_flush();
        }
        this.#h5_flush = null;
        this.#h5_path = null;

        this.#raw_features = null;
        this.#raw_cells = null;
    }

    /**
     * @return {string} Format of this dataset class.
     * @static
     */
    static format() {
        return "10X";
    }

    /**
     * @return {object} Object containing the abbreviated details of this dataset.
     */
    abbreviate() {
        return { 
            "h5": {
                "name": this.#h5_file.name(),
                "size": this.#h5_file.size()
            }
        };
    }

    #features() {
        if (this.#raw_features !== null) {
            return;
        }

        this.#instantiate();
        let handle = new scran.H5File(this.#h5_path);
        if (!("matrix" in handle.children) || handle.children["matrix"] != "Group") {
            throw new Error("expected a 'matrix' group at the top level of the file");
        }
        let mhandle = handle.open("matrix");

        if (!("features" in mhandle.children) || mhandle.children["features"] != "Group") {
            throw new Error("expected a 'matrix/features' group containing the feature annotation");
        }
        let fhandle = mhandle.open("features");

        let ids = eutils.extractHDF5Strings(fhandle, "id");
        if (ids == null) {
            throw new Error("expected a 'matrix/features/id' string dataset containing the feature IDs");
        }
        let feats = { id: ids };

        let names = eutils.extractHDF5Strings(fhandle, "name");
        if (names !== null) {
            feats.name = names;
        }

        let ftype = eutils.extractHDF5Strings(fhandle, "feature_type");
        if (ftype !== null) {
            feats.type = ftype;
        }

        this.#raw_features = new bioc.DataFrame(feats);
    }

    #cells() {
        if (this.#raw_cells !== null) {
            return;
        }

        this.#instantiate();
        let details = scran.extractHDF5MatrixDetails(this.#h5_path, "matrix");
        this.#raw_cells = new bioc.DataFrame({}, { numberOfRows: details.columns });
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.cache=false] - Whether to cache the results for re-use in subsequent calls to this method or {@linkcode TenxHdf5Dataset#load load}.
     * If `true`, users should consider calling {@linkcode TenxHdf5Dataset#clear clear} to release the memory once this dataset instance is no longer needed.
     *
     * @return {object} Object containing the per-feature and per-cell annotations.
     * This has the following properties:
     *
     * - `features`: an object where each key is a modality name and each value is a {@linkplain external:DataFrame DataFrame} of per-feature annotations for that modality.
     * - `cells`: an object containing:
     *   - `number`: the number of cells in this dataset.
     *   - `summary`: an object where each key is the name of a per-cell annotation field and its value is a summary of that annotation field,
     *     following the same structure as returned by {@linkcode summarizeArray}.
     */
    annotations({ cache = false } = {}) {
        this.#features();
        this.#cells();

        let output = {
            "features": futils.reportFeatures(this.#raw_features, "type", cache),
            "cells": {
                "number": this.#raw_cells.numberOfColumns(),
                "summary": {}
            }
        };

        if (!cache) {
            this.clear();
        }
        return output;
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.cache=false] - Whether to cache the results for re-use in subsequent calls to this method or {@linkcode TenxHdf5Dataset#annotations annotations}.
     * If `true`, users should consider calling {@linkcode TenxHdf5Dataset#clear clear} to release the memory once this dataset instance is no longer needed.
     *
     * @return {object} Object containing the per-feature and per-cell annotations.
     * This has the following properties:
     *
     * - `features`: an object where each key is a modality name and each value is a {@linkplain external:DataFrame DataFrame} of per-feature annotations for that modality.
     * - `cells`: a {@linkplain external:DataFrame DataFrame} containing per-cell annotations.
     * - `matrix`: a {@linkplain external:MultiMatrix MultiMatrix} containing one {@linkplain external:ScranMatrix ScranMatrix} per modality.
     * - `row_ids`: an object where each key is a modality name and each value is an integer array containing the feature identifiers for each row in that modality.
     */
    load({ cache = false } = {}) {
        this.#features();
        this.#cells();

        let loaded = scran.initializeSparseMatrixFromHDF5(this.#h5_path, "matrix"); // collection gets handled inside splitScranMatrixAndFeatures.
        let output = futils.splitScranMatrixAndFeatures(loaded, this.#raw_features, "type", cache);
        output.cells = futils.cloneCached(this.#raw_cells, cache);

        if (!cache) {
            this.clear();
        }
        return output;
    }

    /**
     * @return {Array} Array of objects representing the files used in this dataset.
     * Each object corresponds to a single file and contains:
     * - `type`: a string denoting the type.
     * - `file`: a {@linkplain SimpleFile} object representing the file contents.
     */
    serialize() {
        return [ { type: "h5", file: this.#h5_file } ];
    }

    /**
     * @param {Array} files - Array of objects like that produced by {@linkcode TenxHdf5Dataset#serialize serialize}.
     * @return {TenxHdf5Dataset} A new instance of this class.
     * @static
     */
    static async unserialize(files) {
        if (files.length != 1 || files[0].type != "h5") {
            throw new Error("expected exactly one file of type 'h5' for 10X HDF5 unserialization");
        }
        return new TenxHdf5Dataset(files[0].file);
    }
}
