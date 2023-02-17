import * as scran from "scran.js";

export function dumpPcaResultsToHdf5(pcs, path, forceBuffer) {
    let ncells = pcs.numberOfCells();
    let npcs = pcs.numberOfPCs();

    let temppath = scran.chooseTemporaryPath({ extension: ".h5" });
    let content = temppath;

    let fhandle = scran.createNewHDF5File(temppath);
    let buffer;

    try {
        buffer = scran.createFloat64WasmArray(ncells * npcs);

        // Transposing so that cells are the fastest-changing dimension.
        let comp = pcs.principalComponents({ copy: false });
        let arr = buffer.array();
        for (var c = 0; c < ncells; c++) {
            for (var p = 0; p < npcs; p++) {
                arr[c + p * ncells] = comp[p + c * npcs];
            }
        }

        fhandle.writeDataSet(
            "data", 
            "Float64", 
            [npcs, ncells],
            buffer
        ); 

        if (forceBuffer) {
            content = scran.readFile(temppath);
            scran.removeFile(temppath);
        }
    } catch (e) {
        scran.removeFile(temppath);
        throw e;
    } finally {
        scran.free(buffer);
    }

    return { 
        metadata: {
            "$schema": "hdf5_dense_array/v1.json",
            "path": path + "/matrix.h5",
            "array": {
                "dimensions": [ncells, npcs]
            },
            "hdf5_dense_array": {
                "dataset": "data",
            }
        },
        contents: content
    };
}

export function dumpOtherReducedDimensionsToHdf5(dimensions, path, forceBuffer) {
    let ncells = dimensions[0].length;
    let ndims = dimensions.length;

    let temppath = scran.chooseTemporaryPath({ extension: ".h5" });
    let content = temppath;

    let fhandle = scran.createNewHDF5File(temppath);
    let buffer;

    try {
        buffer = scran.createFloat64WasmArray(ncells * ndims);

        for (var d = 0; d < ndims; d++) {
            if (dimensions[d].length !== ncells) {
                throw new Error("all dimensions must have the same length");
            }
            buffer.set(dimensions[d], ncells * d);            
        }

        fhandle.writeDataSet(
            "data", 
            "Float64", 
            [ndims, ncells],
            buffer
        ); 

        if (forceBuffer) {
            content = scran.readFile(temppath);
            scran.removeFile(temppath);
        }
    } catch (e) {
        scran.removeFile(temppath);
        throw e;
    } finally {
        scran.free(buffer);
    }

    return { 
        metadata: {
            "$schema": "hdf5_dense_array/v1.json",
            "path": path + "/matrix.h5",
            "array": {
                "dimensions": [ncells, ndims]
            },
            "hdf5_dense_array": {
                "dataset": "data",
            }
        },
        contents: content
    };
}



