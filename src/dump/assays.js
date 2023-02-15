import * as scran from "scran.js";

export function dumpCountMatrix(mat, path, forceBuffer) {
    let temppath = scran.chooseTemporaryPath({ extension: ".h5" });
    let contents = temppath;

    try {
        scran.writeSparseMatrixToHdf5(mat, temppath, "matrix", { format: "tenx_matrix" });
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
            "$schema": "hdf5_sparse_matrix/v1.json",
            "path": path + "/matrix.h5",
            "array": {
                "dimensions": [mat.numberOfRows(), mat.numberOfColumns()]
            },
            "hdf5_sparse_matrix": {
                "group": "matrix",
                "format": "tenx_matrix"
            }
        },
        contents: contents
    };
}
