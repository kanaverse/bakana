import * as jsp from "jaspagate";
import * as bioc from "bioconductor";

export class MockReducedDimensionMatrix {
    #values;
    #nr;
    #nc;

    constructor(nr, nc, values) {
        this.#nr = nr;
        this.#nc = nc;
        this.#values = values;
    }

    _bioconductor_NUMBER_OF_ROWS() {
        return this.#nr;
    }

    _bioconductor_NUMBER_OF_COLUMNS() {
        return this.#nc;
    }

    get values() {
        return this.#values;
    }
}

export async function saveReducedDimensionMatrix(x, path, globals, options) {
    await globals.mkdir(path);
    let fhandle = await globals.h5create(jsp.joinPath(path, "array.h5"));
    let handle_stack = [fhandle]
    let success = false;
    try {
        let dhandle = fhandle.createGroup("dense_array");
        handle_stack.push(dhandle);
        dhandle.createDataSet("data", "Float64", [bioc.NUMBER_OF_COLUMNS(x), bioc.NUMBER_OF_ROWS(x)], { data: x.values }); // [ncol, nrow] as HDF5 uses row-major storage.
        dhandle.writeAttribute("type", "String", [], ["number"]);
        dhandle.writeAttribute("transposed", "Int8", [], [1]);
        success = true;
    } finally {
        for (const handle of handle_stack.reverse()) {
            handle.close();
        }
        await globals.h5finish(fhandle, !success);
    }

    await globals.write(jsp.joinPath(path, "OBJECT"), JSON.stringify({
        type: "dense_array",
        dense_array: { version: "1.0" }
    }));
}
