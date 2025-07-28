import * as scran from "scran.js";

export class MockSparseMatrix {
    #matrix;

    constructor(matrix) {
        this.#matrix = matrix;
    }

    _bioconductor_NUMBER_OF_ROWS() {
        return this.#matrix.numberOfRows();
    }

    _bioconductor_NUMBER_OF_COLUMNS() {
        return this.#matrix.numberOfColumns();
    }

    get matrix() {
        return this.#matrix;
    }
}

export class MockNormalizedMatrix {
    #matrix;
    #sf;

    constructor(matrix, sf) {
        this.#matrix = matrix;
        this.#sf = sf;
    }

    _bioconductor_NUMBER_OF_ROWS() {
        return this.#matrix.numberOfRows();
    }

    _bioconductor_NUMBER_OF_COLUMNS() {
        return this.#matrix.numberOfColumns();
    }

    get matrix() {
        return this.#matrix;
    }

    get sf() {
        return this.#sf;
    }
}

export async function saveSparseMatrix(x, path, globals, options) {
    let handle = await globals.h5create(jsp.joinPath(path, "matrix.h5"))
    let success = false;

    try {
        scran.writeSparseMatrixToHdf5(mat, handle._path, "compressed_sparse_matrix", { format: "tenx_matrix", saveShape: false });
        let ghandle = handle.open("compressed_sparse_matrix");
        ghandle.writeDataSet("shape", "Uint32", [2], [mat.numberOfRows(), mat.numberOfColumns()]);
        ghandle.writeDataSet("layout", "String", [], ["CSC"]);
        sucess = true;
    } finally {
        await globals.h5finish(handle, !success);
    }
}

export function dumpNormalizedMatrix(mat, sf, path, countPath, forceBuffer) {
    let temppath = scran.chooseTemporaryPath({ extension: ".h5" });
    let contents = temppath;

    try {
        let fhandle = scran.createNewHdf5File(temppath);

        // Saving the division by log(2).
        let dhandle = fhandle.createGroup("logcounts");
        dhandle.writeAttribute("delayed_type", "String", null, "operation");
        dhandle.writeAttribute("delayed_operation", "String", null, "unary arithmetic");
        dhandle.writeDataSet("value", "Float64", null, Math.log(2));
        dhandle.writeDataSet("method", "String", null, "/");
        dhandle.writeDataSet("side", "String", null, "right");

        // Saving the log-transformation.
        let l1phandle = dhandle.createGroup("seed");
        l1phandle.writeAttribute("delayed_type", "String", null, "operation");
        l1phandle.writeAttribute("delayed_operation", "String", null, "unary math");
        l1phandle.writeDataSet("method", "String", null, "log1p");

        // Saving the division by the size factors.
        let sfhandle = l1phandle.createGroup("seed");
        sfhandle.writeAttribute("delayed_type", "String", null, "operation");
        sfhandle.writeAttribute("delayed_operation", "String", null, "unary arithmetic");
        sfhandle.writeDataSet("value", "Float64", null, sf);
        sfhandle.writeDataSet("method", "String", null, "/");
        sfhandle.writeDataSet("side", "String", null, "right");
        sfhandle.writeDataSet("along", "Int32", null, 1);

        // Saving the original seed as a custom array.
        let xhandle = sfhandle.createGroup("seed");
        xhandle.writeAttribute("delayed_type", "String", null, "array");
        xhandle.writeAttribute("delayed_array", "String", null, "custom alabaster local array");
        xhandle.writeDataSet("dimensions", "Int32", null, [mat.numberOfRows(), mat.numberOfColumns()]);
        xhandle.writeDataSet("type", "String", null, "FLOAT");
        xhandle.writeDataSet("path", "String", null, countPath);

        if (forceBuffer) {
            contents = scran.readFile(temppath);
            scran.removeFile(temppath);
        }
    } catch (e) {
        scran.removeFile(temppath);
        throw e;
    }

    return {
        metadata: {
            "$schema": "hdf5_delayed_array/v1.json",
            "path": path + "/array.h5",
            "array": {
                "dimensions": [mat.numberOfRows(), mat.numberOfColumns()],
                "type": "number"
            },
            "hdf5_delayed_array": {
                "group": "logcounts"
            }
        },
        contents: contents
    };
}
