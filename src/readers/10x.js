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

    #options;

    #dump_summary(fun) {
        let files = [{ type: "h5", file: fun(this.#h5_file) }];
        let options = this.options();
        return { files, options };
    }

    /**
     * @param {SimpleFile|string|Uint8Array|File} h5File - Contents of a HDF5 file in the 10X feature-barcode format.
     * On browsers, this may be a File object.
     * On Node.js, this may also be a string containing a file path.
     */
    constructor(h5File) {
        if (h5File instanceof afile.SimpleFile) {
            this.#h5_file = h5File;
        } else {
            this.#h5_file = new afile.SimpleFile(h5File);
        }

        this.#options = TenxHdf5Dataset.defaults();
        this.clear();
    }

    /**
     * @return {object} Default options, see {@linkcode TenxHdf5Dataset#setOptions setOptions} for more details.
     */
    static defaults() {
        return {
            featureTypeRnaName: "Gene Expression", 
            featureTypeAdtName: "Antibody Capture", 
            featureTypeCrisprName: "CRISPR Guide Capture", 
            primaryRnaFeatureIdColumn: 0, 
            primaryAdtFeatureIdColumn: 0,
            primaryCrisprFeatureIdColumn: 0
        };
    }

    /**
     * @return {object} Object containing all options used for loading.
     */
    options() {
        return { ...(this.#options) };
    }

    /**
     * @param {object} options - Optional parameters that affect {@linkcode TenxHdf5Dataset#load load} (but not {@linkcode TenxHdf5Dataset#summary summary}.
     * @param {?string} [options.featureTypeRnaName] - Name of the feature type for gene expression.
     * Alternatively `null`, to indicate that no RNA features are to be loaded.
     * @param {?string} [options.featureTypeAdtName] - Name of the feature type for ADTs.
     * Alternatively `null`, to indicate that no ADT features are to be loaded.
     * @param {?string} [options.featureTypeCrisprName] - Name of the feature type for CRISPR guides.
     * Alternatively `null`, to indicate that no guides are to be loaded.
     * @param {string|number} [options.primaryRnaFeatureIdColumn] - Name or index of the column of the `features` {@linkplain external:DataFrame DataFrame} that contains the primary feature identifier for gene expression.
     * If `i` is invalid (e.g., out of range index, unavailable name), it is ignored and the primary identifier is treated as undefined.
     * @param {string|number} [options.primaryAdtFeatureIdColumn] - Name or index of the column of the `features` {@linkplain external:DataFrame DataFrame} that contains the primary feature identifier for the ADTs.
     * If `i` is invalid (e.g., out of range index, unavailable name), it is ignored and the primary identifier is treated as undefined.
     * @param {string|number} [options.primaryCrisprFeatureIdColumn] - Name or index of the column of the `features` {@linkplain external:DataFrame DataFrame} that contains the primary feature identifier for the CRISPR guides.
     * If `i` is invalid (e.g., out of range index, unavailable name), it is ignored and the primary identifier is treated as undefined.
     */
    setOptions(options) {
        for (const [k, v] of Object.entries(options)) {
            this.#options[k] = v;
        }
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
     * This may be called at any time but only has an effect if `cache = true` in {@linkcode TenxHdf5Dataset#load load} or {@linkcodeTenxHdf5Dataset#summary summary}. 
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
     * @return {object} Object containing the abbreviated details of this dataset,
     * in a form that can be cheaply stringified.
     */
    abbreviate() {
        return this.#dump_summary(f => { return { name: f.name(), size: f.size() }; });
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
        let feats = new bioc.DataFrame({ id: ids }); // build it piece-by-piece for a well-defined ordering.

        let names = eutils.extractHDF5Strings(fhandle, "name");
        if (names !== null) {
            feats.$setColumn("name", names);
        }

        let ftype = eutils.extractHDF5Strings(fhandle, "feature_type");
        if (ftype !== null) {
            feats.$setColumn("type", ftype);
        }

        this.#raw_features = feats;
        return;
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
     * - `modality_features`: an object where each key is a modality name and each value is a {@linkplain external:DataFrame DataFrame} of per-feature annotations for that modality.
     *   Unlike {@linkcode TenxMatrixMarketDataset#load load}, modality names are arbitrary.
     * - `cells`: a {@linkplain external:DataFrame DataFrame} of per-cell annotations.
     */
    summary({ cache = false } = {}) {
        this.#features();
        this.#cells();

        let output = {
            "modality_features": futils.reportFeatures(this.#raw_features, "type"),
            "cells": this.#raw_cells
        };

        if (!cache) {
            this.clear();
        }
        return output;
    }

    /**
     * @param {object} [options={}] - Optional parameters.
     * @param {boolean} [options.cache=false] - Whether to cache the results for re-use in subsequent calls to this method or {@linkcode TenxHdf5Dataset#summary summary}.
     * If `true`, users should consider calling {@linkcode TenxHdf5Dataset#clear clear} to release the memory once this dataset instance is no longer needed.
     *
     * @return {object} Object containing the per-feature and per-cell annotations.
     * This has the following properties:
     *
     * - `features`: an object where each key is a modality name and each value is a {@linkplain external:DataFrame DataFrame} of per-feature annotations for that modality.
     * - `cells`: a {@linkplain external:DataFrame DataFrame} containing per-cell annotations.
     * - `matrix`: a {@linkplain external:MultiMatrix MultiMatrix} containing one {@linkplain external:ScranMatrix ScranMatrix} per modality.
     * - `primary_ids`: an object where each key is a modality name and each value is an array of strings containing the feature identifiers for each row in that modality.
     *
     * Modality names are guaranteed to be one of `"RNA"`, `"ADT"` or `"CRISPR"`.
     * It is assumed that an appropriate mapping from the feature types inside the `featureFile` was previously declared,
     * either in the constructor or in setters like {@linkcode setFeatureTypeRnaName}.
     *
     * If the feature annotation lacks information about the feature types, it is assumed that all features are genes, i.e., only the RNA modality is present.
     */
    load({ cache = false } = {}) {
        this.#features();
        this.#cells();

        let loaded = scran.initializeSparseMatrixFromHDF5(this.#h5_path, "matrix"); // collection gets handled inside splitScranMatrixAndFeatures.

        let mappings = { 
            RNA: this.#options.featureTypeRnaName, 
            ADT: this.#options.featureTypeAdtName,
            CRISPR: this.#options.featureTypeCrisprName
        };
        let output = futils.splitScranMatrixAndFeatures(loaded, this.#raw_features, "type", mappings, "RNA");
        output.cells = this.#raw_cells;

        let primaries = { 
            RNA: this.#options.primaryRnaFeatureIdColumn, 
            ADT: this.#options.primaryAdtFeatureIdColumn,
            CRISPR: this.#options.primaryCrisprFeatureIdColumn
        };
        output.primary_ids = futils.extractPrimaryIds(output.features, primaries);

        if (!cache) {
            this.clear();
        }
        return output;
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
     * @param {Array} files - Array of objects like that produced by {@linkcode TenxHdf5Dataset#serialize serialize}.
     * @param {object} options - Object containing additional options to be passed to the constructor.
     * @return {TenxHdf5Dataset} A new instance of this class.
     * @static
     */
    static async unserialize(files, options) {
        if (files.length != 1 || files[0].type != "h5") {
            throw new Error("expected exactly one file of type 'h5' for 10X HDF5 unserialization");
        }
        return new TenxHdf5Dataset(files[0].file, options);
    }
}
