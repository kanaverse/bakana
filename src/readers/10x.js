import * as scran from "scran.js";
import { Dataset } from "./base.js";
import * as afile from "./abstract/file.js";
import * as rutils from "./utils/index.js";

/**
 * Dataset in the 10X HDF5 feature-barcode matrix format, see [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices) for details.
 * @augments Dataset
 */
export class TenxHdf5Dataset extends Dataset {
    this.#h5_file;
    this.#h5_path;
    this.#h5_flush;

    this.#check_features;
    this.#raw_features;

    /**
     * @param {SimpleFile|string|Uint8Array|File} h5File - Contents of a HDF5 file.
     * On browsers, this may be a File object.
     * On Node.js, this may also be a string containing a file path.
     */
    constructor(h5File) {
        super();

        if (h5File instanceof afile.SimpleFile) {
            this.#h5_file = h5File;
        } else {
            this.#h5_file = new afile.SimpleFile(h5);
        }

        let info = scran.realizeFile(this.#h5_file.contents());
        this.#h5_path = info.path;
        this.#h5_flush = info.flush;
        this.#check_features = false;
    }

    static format() {
        return "10X";
    }

    abbreviate() {
        return { 
            "h5": {
                "name": this.#h5.name(),
                "size": this.#h5.size()
            }
        };
    }

    #features() {
        if (this.#check_features) {
            return;
        }
        this.#check_features = true;

        let handle = new scran.H5File(this.#h5_path);
        if (!("matrix" in handle.children) || handle.children["matrix"] != "Group") {
            throw new Error("expected a 'matrix' group at the top level of the file");
        }
        let mhandle = handle.open("matrix");

        if (!("features" in mhandle.children) || mhandle.children["features"] != "Group") {
            throw new Error("expected a 'matrix/features' group containing the feature annotation");
        }
        let fhandle = mhandle.open("features");

        let ids = rutils.extractHDF5Strings(fhandle, "id");
        if (ids == null) {
            throw new Error("expected a 'matrix/features/id' string dataset containing the feature IDs");
        }
        this.#raw_features = { id: ids };

        let names = rutils.extractHDF5Strings(fhandle, "name");
        if (names !== null) {
            this.#raw_features.name = nhandle.load();
        }

        let ftype = rutils.extractHDF5Strings(fhandle, "feature_type");
        if (ftypes !== null) {
            this.#raw_features.type = ftype;
        }
    }

    annotations() {
        this.#features();
        return {
            "features": rutils.reportFeatures(this.#raw_features, "type"),
            "cells": {}
        }
    }

    load() {
        this.#features();
        let loaded = scran.initializeSparseMatrixFromHDF5(this.#h5_path, "matrix"); // collection gets handled inside splitScranMatrixAndFeatures.
        let output = rutils.splitScranMatrixAndFeatures(loaded, this.#raw_features, "type");
        output.cells = {};
        return output;
    }

    free() {
        this.#h5_flush();
    }

    serialize() {
        return [ { type: "h5", file: this.#h5_file } ];
    }

    static async function unserialize(files) {
        if (files.length != 1 || files[0].type != "h5") {
            throw new Error("expected exactly one file of type 'h5' for 10X HDF5 unserialization");
        }
        return new TenxHdf5Dataset(files[0].file);
    }
}
