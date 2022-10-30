import * as scran from "scran.js";
import * as bioc from "bioconductor";
import { Dataset } from "./base.js";
import * as afile from "./abstract/file.js";
import * as eutils from "./utils/extract.js";
import * as futils from "./utils/features.js";

/**
 * Dataset in the 10X HDF5 feature-barcode matrix format, see [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices) for details.
 * @augments Dataset
 */
export class TenxHdf5Dataset extends Dataset {
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
        super();

        if (h5File instanceof afile.SimpleFile) {
            this.#h5_file = h5File;
        } else {
            this.#h5_file = new afile.SimpleFile(h5File);
        }

        let info = scran.realizeFile(this.#h5_file.content());

        this.#h5_path = info.path;
        this.#h5_flush = info.flush;

        this.#raw_features = null;
        this.#raw_cells = null;
    }

    static format() {
        return "10X";
    }

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
        let details = scran.extractHDF5MatrixDetails(this.#h5_path, "matrix");
        this.#raw_cells = new bioc.DataFrame({}, { numberOfRows: details.columns });
    }

    annotations() {
        this.#features();
        this.#cells();

        return {
            "features": futils.reportFeatures(this.#raw_features, "type"),
            "cells": {
                "number": this.#raw_cells.numberOfColumns(),
                "summary": {}
            }
        }
    }

    load() {
        this.#features();
        this.#cells();

        let loaded = scran.initializeSparseMatrixFromHDF5(this.#h5_path, "matrix"); // collection gets handled inside splitScranMatrixAndFeatures.
        let output = futils.splitScranMatrixAndFeatures(loaded, this.#raw_features, "type");
        output.cells = bioc.CLONE(this.#raw_cells);
        return output;
    }

    free() {
        this.#h5_flush();
    }

    serialize() {
        return [ { type: "h5", file: this.#h5_file } ];
    }

    static async unserialize(files) {
        if (files.length != 1 || files[0].type != "h5") {
            throw new Error("expected exactly one file of type 'h5' for 10X HDF5 unserialization");
        }
        return new TenxHdf5Dataset(files[0].file);
    }
}
