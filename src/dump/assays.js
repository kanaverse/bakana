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
    let has_ext_store = "external_Matrix_store" in options;
    if (has_ext_store) {
        for (const [y, ypath] of options.external_Matrix_store) {
            if (y === x.matrix) {
                for (const f of ["OBJECT", "matrix.h5"]) {
                    await globals.copy(jsp.joinPath(ypath, f), jsp.joinPath(path, f));
                }
                return;
            }
        }
    }

    let fhandle = await globals.h5create(jsp.joinPath(path, "matrix.h5"))
    let handle_stack = [fhandle];
    let success = false;
    try {
        scran.writeSparseMatrixToHdf5(x.matrix, fhandle._path, "compressed_sparse_matrix", { format: "tenx_matrix", saveShape: false, overwrite: false });
        let ghandle = fhandle.open("compressed_sparse_matrix");
        handle_stack.push(ghandle);
        ghandle.writeDataSet("shape", "Uint32", [2], [x.matrix.numberOfRows(), x.matrix.numberOfColumns()]);
        ghandle.writeAttribute("layout", "String", [], ["CSC"]);
        sucess = true;
    } finally {
        for (const handle of handle_stack.reverse()) {
            handle.close();
        }
        await globals.h5finish(handle, !success);
    }

    await globals.write(jsp.joinPath(path, "OBJECT"), JSON.stringify({
        type: "compressed_sparse_matrix",
        compressed_sparse_matrix: { version: "1.0" }
    }));

    if ("external_Matrix_store" in options) {
        options.external_Matrix_store.push([x, path]);
    }
}

export async function saveNormalizedMatrix(x, path, globals, options) {
    let fhandle = await globals.h5create(jsp.joinPath(path, "array.h5"));
    let handle_stack = [fhandle];
    let success = false;
    try {
        // Saving the division by log(2).
        let dhandle = fhandle.createGroup("logcounts");
        handle_stack.push(dhandle);
        dhandle.writeAttribute("delayed_type", "String", [], "operation");
        dhandle.writeAttribute("delayed_operation", "String", [], "unary arithmetic");
        dhandle.writeDataSet("value", "Float64", [], { data: [Math.log(2)] });
        dhandle.writeDataSet("method", "String", [], { data: ["/"] });
        dhandle.writeDataSet("side", "String", [], { data: ["right"] });

        // Saving the log-transformation.
        let l1phandle = dhandle.createGroup("seed");
        handle_stack.push(l1handle);
        l1phandle.writeAttribute("delayed_type", "String", [], { data: ["operation"] });
        l1phandle.writeAttribute("delayed_operation", "String", [], { data: ["unary math"] });
        l1phandle.writeDataSet("method", "String", [], { data: ["log1p"] });

        // Saving the division by the size factors.
        let sfhandle = l1phandle.createGroup("seed");
        handle_stack.push(sfhandle);
        sfhandle.writeAttribute("delayed_type", "String", [], { data: ["operation"] });
        sfhandle.writeAttribute("delayed_operation", "String", [], { data: ["unary arithmetic"] });
        sfhandle.writeDataSet("value", "Float64", [x.sf.length], { data: x.sf });
        sfhandle.writeDataSet("method", "String", [], { data: ["/"] });
        sfhandle.writeDataSet("side", "String", [], { data: ["right"] });
        sfhandle.writeDataSet("along", "Int32", [], { data: [1] });

        // Saving the original seed as a custom array.
        let xhandle = sfhandle.createGroup("seed");
        handle_stack.push(xhandle);
        xhandle.writeAttribute("delayed_type", "String", [], { data: ["array"] });
        xhandle.writeAttribute("delayed_array", "String", [], { data: ["custom takane seed array"] });
        xhandle.writeDataSet("dimensions", "Int32", [2], { data: [mat.numberOfRows(), mat.numberOfColumns()] });
        xhandle.writeDataSet("type", "String", [], { data: ["FLOAT"] });
        xhandle.writeDataSet("index", "Uint8", [], { data: [0] });

        let seed_dir = jsp.joinPath(path, "seeds");
        await globals.mkdir(seed_dir);
        await saveSparseMatrix(x.matrix, jsp.joinPath(seed_dir, "0"), globals, options);
        success = true;

    } finally {
        for (const handle of handle_stack.reverse()) {
            handle.close();
        }
        await globals.h5finish(handle, !success);
    }

    await globals.write(jsp.joinPath(path, "OBJECT"), JSON.stringify({
        type: "delayed_array",
        delayed_array: { version: "1.0" }
    }));
}
