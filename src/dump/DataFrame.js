import * as scran from "scran.js";
import * as wa from "wasmarrays.js";
import * as bioc from "bioconductor";

// Monkey-patching these methods so that we can use these WasmArrays
// as columns in a bioc.DataFrame.
wa.Uint8WasmArray.prototype._bioconductor_LENGTH = function() { return this.length; };
wa.Int32WasmArray.prototype._bioconductor_LENGTH = function() { return this.length; };
wa.Float64WasmArray.prototype._bioconductor_LENGTH = function() { return this.length; };

export function writeHdf5DataFrame(x, path, { group = "data", forceBuffer = false } = {}) {
    let metadata = {
        "path": path + "/simple.h5",
        "$schema": "hdf5_data_frame/v1.json",
        "data_frame": {
            "dimensions": [ x.numberOfRows(), x.numberOfColumns() ],
            "columns": [],
            "row_names": false
        },
        "hdf5_data_frame": {
            "group": group
        }
    };

    let temppath = scran.chooseTemporaryPath({ extension: ".h5" });
    let contents = temppath;
    let children = [];

    let fhandle = scran.createNewHdf5File(temppath);
    try {
        let ghandle = fhandle.createGroup(group);

        ghandle.writeDataSet("column_names", "String", null, x.columnNames());
        let rn = x.rowNames();
        if (rn !== null) {
            metadata.data_frame.row_names = true;
            ghandle.writeDataSet("row_names", "String", null, rn);
        }

        let dhandle = ghandle.createGroup("data");
        let coltypes = metadata.data_frame.columns;

        for (var i = 0; i < x.numberOfColumns(); i++) {
            const curcol = x.column(i);
            const colname = x.columnNames()[i];

            if (curcol instanceof Array) {
                let all_types = new Set;
                let has_null = false;
                for (const y of curcol) {
                    if (y === null) {
                        has_null = true;
                    } else {
                        all_types.add(typeof y);
                    }
                }

                if (all_types.size > 1) {
                    throw new Error("column '" + colname + "' has multiple types");
                }

                if (all_types.has("string")) {
                    coltypes.push({ name: colname, type: "string" });
                    if (!has_null) {
                        dhandle.writeDataSet(String(i), "String", null, curcol);
                    } else {
                        let contents = new Set(curcol);
                        let placeholder = "NA";
                        while (contents.has(placeholder)) {
                            placeholder = placeholder + "_";
                        }

                        let copy = curcol.slice();
                        for (var i = 0; i < copy.length; i++) {
                            if (copy[i] === null) {
                                copy[i] = placeholder;
                            }
                        }

                        let shandle = dhandle.writeDataSet(String(i), "String", null, copy);
                        shandle.writeAttribute("missing-value-placeholder", "String", null, placeholder);
                    }

                } else if (all_types.has("number")) {
                    coltypes.push({ name: colname, type: "number" });
                    let temp = scran.createFloat64WasmArray(curcol.length);
                    try {
                        if (!has_null) {
                            temp.set(curcol);
                        } else {
                            let temparr = temp.array();
                            for (var i = 0; i < curcol.length; i++) {
                                if (curcol[i] === null) {
                                    temparr[i] = Number.NaN;
                                } else {
                                    temparr[i] = curcol[i];
                                }
                            }
                        }
                        dhandle.writeDataSet(String(i), "Float64", null, temp);
                    } finally {
                        temp.free();
                    }

                } else if (all_types.size == 0 || all_types.has("boolean")) {
                    coltypes.push({ name: colname, type: "boolean" });
                    if (has_null) {
                        let temp = scran.createInt32WasmArray(curcol.length);
                        try {
                            for (var i = 0; i < curcol.length; i++) {
                                if (curcol[i] === null) {
                                    temparr[i] = -2147483648;
                                } else {
                                    temparr[i] = curcol[i];
                                }
                            }
                            dhandle.writeDataSet(String(i), "Int32", null, temp);
                        } finally {
                            temp.free();
                        }
                    } else {
                        let temp = scran.createUint8WasmArray(curcol.length);
                        try {
                            temp.fill(curcol);
                            dhandle.writeDataSet(String(i), "Uint8", null, temp);
                        } finally {
                            temp.free();
                        }
                    }

                } else {
                    throw new Error("unknown type '" + Array.from(all_types)[0] + "' for column '" + colname + "'");
                }

            } else if (curcol instanceof Uint8Array || curcol instanceof wa.Uint8WasmArray) {
                coltypes.push({ name: colname, type: "boolean" });
                dhandle.writeDataSet(String(i), "Uint8", null, curcol);
                
            } else if (curcol instanceof Int32Array || curcol instanceof wa.Int32WasmArray) {
                coltypes.push({ name: colname, type: "integer" });
                dhandle.writeDataSet(String(i), "Int32", null, curcol);

            } else if (curcol instanceof Float64Array || curcol instanceof wa.Float64WasmArray) {
                coltypes.push({ name: colname, type: "number" });
                dhandle.writeDataSet(String(i), "Float64", null, curcol);

            } else if (curcol instanceof bioc.DataFrame) {
                let subpath = path + "/column" + String(i);
                let child = writeHdf5DataFrame(curcol, subpath, { group, forceBuffer });
                coltypes.push({ name: colname, type: "other", resource: { type: "local", path: child.self.metadata.path } });
                child.self.metadata.is_child = true;
                children.push(child.self);
                for (const x of child.children) {
                    children.push(x);
                }

            } else {
                throw new Error("unknown type for column '" + colname + "'");
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
        self: {
            metadata: metadata,
            contents: contents
        },
        children: children
    };
}
