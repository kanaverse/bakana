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
                "dimensions": [mat.numberOfRows(), mat.numberOfColumns()],
                "type": "integer"
            },
            "hdf5_sparse_matrix": {
                "group": "matrix",
                "format": "tenx_matrix"
            }
        },
        contents: contents
    };
}

export function dumpDelayedCountMatrix(state, modality, handlers, path, forceBuffer) {
    let temppath = scran.chooseTemporaryPath({ extension: ".h5" });
    let contents = temppath;

    // Obtaining the subset.
    let mat = state.cell_filtering.fetchFilteredMatrix().get(modality);
    let ncells = mat.numberOfColumns();
    let current = new Int32Array(ncells);
    for (var i = 0; i < ncells; i++) {
        current[i] = i;
    }
    state.cell_filtering.undoFilter(current);
    state.inputs.undoSubset(current);

    // Obtaining the order of blocks.
    let block_order = [];
    let all_ds = state.inputs.fetchDatasetDetails();
    let dsnames = all_ds.keys();
    if (dsnames.length == 1) {
        block_order.push(dsnames[0]);
    } else {
        block_order = state.inputs.fetchBlockLevels();
    }

    let save_seed = (handle, name) => {
        handle.writeAttribute("delayed_type", "String", null, "array");
        handle.writeDataSet("type", "String", null, "INTEGER");
        handle.writeDataSet("dimensions", "Int32", null, [all_ds[name].features[modality], all_ds[name].cells]);
        if (!(name in handlers)) {
            throw new Error("could not find dataset '" + name + "' in delayed override handlers");
        }
        handlers[name](handle);
    };

    try {
        let fhandle = scran.createNewHDF5File(temppath);

        // Saving the subsetting operation to get to the QC'd matrix.
        let dhandle = fhandle.createGroup("counts");
        dhandle.writeAttribute("delayed_type", "String", null, "operation");
        dhandle.writeAttribute("delayed_operation", "String", null, "subset");
        let ghandle = dhandle.createGroup("index");
        ghandle.writeAttribute("delayed_type", "String", null, "list");
        ghandle.writeAttribute("delayed_length", "Int32", null, 2);
        ghandle.writeDataSet("1", "Int32", null, current);

        // Saving the seed, or the bound matrix.
        let xhandle = sfhandle.createGroup("seed");
        if (block_order.length == 1) {
            save_seed(xhandle, block_order[0]);
        } else {
            xhandle.writeAttribute("delayed_type", "String", null, "operation");
            xhandle.writeAttribute("delayed_operation", "String", null, "combine");
            xhandle.writeDataSet("along", "Int32", null, 1);

            let shandle = dhandle.createGroup("seeds");
            for (var b = 0; b < block_order.length; b++) {
                let curhandle = shandle.createGroup(String(b));
                save_seed(curhandle, block_order[b]);
            }
        }

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
                "type": "integer"
            },
            "hdf5_delayed_array": {
                "group": "counts"
            }
        },
        contents: contents
    };
}

export function dumpNormalizedMatrix(mat, sf, path, countPath, forceBuffer) {
    let temppath = scran.chooseTemporaryPath({ extension: ".h5" });
    let contents = temppath;

    try {
        let fhandle = scran.createNewHDF5File(temppath);

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
