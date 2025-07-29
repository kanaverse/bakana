import * as scran from "scran.js";
import * as jsp from "jaspagate";

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

    await globals.mkdir(path);
    let fhandle = await globals.h5create(jsp.joinPath(path, "matrix.h5"))
    let handle_stack = [fhandle];
    let success = false;
    try {
        scran.writeSparseMatrixToHdf5(x.matrix, fhandle._path, "compressed_sparse_matrix", { format: "tenx_matrix", saveShape: false, overwrite: false });
        let kids = fhandle.children();
        kids["compressed_sparse_matrix"] = "Group"; // need to trick it a little to get open() to work.
        let ghandle = fhandle.open("compressed_sparse_matrix");
        handle_stack.push(ghandle);
        ghandle.createDataSet("shape", "Uint32", [2], { data: [x.matrix.numberOfRows(), x.matrix.numberOfColumns()] }).close();
        ghandle.writeAttribute("layout", "String", [], ["CSC"]);
        ghandle.writeAttribute("type", "String", [], ["integer"]); // assume that we're always dealing with count matrices.
        success = true;
    } finally {
        for (const handle of handle_stack.reverse()) {
            handle.close();
        }
        await globals.h5finish(fhandle, !success);
    }

    await globals.write(jsp.joinPath(path, "OBJECT"), JSON.stringify({
        type: "compressed_sparse_matrix",
        compressed_sparse_matrix: { version: "1.0" }
    }));

    if ("external_Matrix_store" in options) {
        options.external_Matrix_store.push([x.matrix, path]);
    }
}

export async function saveNormalizedMatrix(x, path, globals, options) {
    await globals.mkdir(path);
    let fhandle = await globals.h5create(jsp.joinPath(path, "array.h5"));
    let handle_stack = [fhandle];
    let success = false;
    try {
        // Saving the division by log(2).
        let dhandle = fhandle.createGroup("delayed_array");
        handle_stack.push(dhandle);
        dhandle.writeAttribute("delayed_type", "String", [], ["operation"]);
        dhandle.writeAttribute("delayed_operation", "String", [], ["unary arithmetic"]);
        dhandle.writeAttribute("delayed_version", "String", [], ["1.1"]);
        let vhandle = dhandle.createDataSet("value", "Float64", [], { data: [Math.log(2)] });
        vhandle.writeAttribute("type", "String", [], ["FLOAT"]);
        handle_stack.push(vhandle);
        dhandle.createDataSet("method", "String", [], { data: ["/"] }).close();
        dhandle.createDataSet("side", "String", [], { data: ["right"] }).close();

        // Saving the log-transformation.
        let l1phandle = dhandle.createGroup("seed");
        handle_stack.push(l1phandle);
        l1phandle.writeAttribute("delayed_type", "String", [], ["operation"]);
        l1phandle.writeAttribute("delayed_operation", "String", [], ["unary math"]);
        l1phandle.createDataSet("method", "String", [], { data: ["log1p"] }).close();

        // Saving the division by the size factors.
        let sfhandle = l1phandle.createGroup("seed");
        handle_stack.push(sfhandle);
        sfhandle.writeAttribute("delayed_type", "String", [], ["operation"]);
        sfhandle.writeAttribute("delayed_operation", "String", [], ["unary arithmetic"]);
        vhandle = sfhandle.createDataSet("value", "Float64", [x.sf.length], { data: x.sf });
        vhandle.writeAttribute("type", "String", [], ["FLOAT"]);
        handle_stack.push(vhandle);
        sfhandle.createDataSet("method", "String", [], { data: ["/"] }).close();
        sfhandle.createDataSet("side", "String", [], { data: ["right"] }).close();
        sfhandle.createDataSet("along", "Uint8", [], { data: [1] }).close();

        // Saving the original seed as a custom array.
        let xhandle = sfhandle.createGroup("seed");
        handle_stack.push(xhandle);
        xhandle.writeAttribute("delayed_type", "String", [], ["array"]);
        xhandle.writeAttribute("delayed_array", "String", [], ["custom takane seed array"]);
        xhandle.createDataSet("dimensions", "Uint64", [2], { data: [x.matrix.numberOfRows(), x.matrix.numberOfColumns()] }).close();
        xhandle.createDataSet("type", "String", [], { data: ["INTEGER"] }).close();
        xhandle.createDataSet("index", "Uint8", [], { data: [0] }).close();

        let seed_dir = jsp.joinPath(path, "seeds");
        await globals.mkdir(seed_dir);
        await saveSparseMatrix(new MockSparseMatrix(x.matrix), jsp.joinPath(seed_dir, "0"), globals, options);
        success = true;

    } finally {
        for (const handle of handle_stack.reverse()) {
            handle.close();
        }
        await globals.h5finish(fhandle, !success);
    }

    await globals.write(jsp.joinPath(path, "OBJECT"), JSON.stringify({
        type: "delayed_array",
        delayed_array: { version: "1.0" }
    }));
}
