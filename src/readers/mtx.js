import * as scran from "scran.js";
import * as rutils from "./utils/index.js";
import * as afile from "../abstract/file.js";

export function abbreviate(args) {
    var formatted = {
        "format": "MatrixMarket", 
        "mtx": rutils.formatFile(args.mtx, true)
    };

    if ("genes" in args) {
        formatted.genes = rutils.formatFile(args.genes, true);
    }

    if ("annotations" in args) {
        formatted.annotations = rutils.formatFile(args.annotations, true);
    }

    return formatted;
}

function extract_features(gene_file, numberOfRows) {
    const content = new Uint8Array(gene_file.content.buffer());
    var is_gz = gene_file.name.endsWith(".gz");
    let parsed = rutils.readTable(content, { compression: (is_gz ? "gz" : "none") });

    if (parsed.length == numberOfRows + 1) {
        // If it seems to have a header, we just use that directly.
        let output = {};
        let headers = parsed.shift();
        headers.forEach((x, i) => {
            output[x] = parsed.map(y => y[i]);
        });
        return output;

    } else {
        // Otherwise, we assume it's standard 10X CellRanger output, without a header.
        if (parsed.length !== numberOfRows) {
            throw new Error("number of matrix rows is not equal to the number of rows in '" + gene_file.name + "'");
        } 

        var ids = [], symb = [];
        parsed.forEach(x => {
            ids.push(x[0]);
            symb.push(x[1]);
        });
        let output = { "id": ids, "symbol": symb };

        if (parsed[0].length >= 3) {
            let types = [];
            parsed.forEach(x => { types.push(x[2]); });
            output.type = types;
        }

        return output;
    }
}

function extract_annotations(annotation_file, numberOfColumns, { summary = false, summaryLimit = 50 } = {}) {
    const content = new Uint8Array(annotation_file.content.buffer());
    var is_gz = annotation_file.name.endsWith(".gz");
    let parsed = rutils.readTable(content, { compression: (is_gz ? "gz" : "none") });

    // Check if a header is present or not. Standard 10X output doesn't have a 
    // header but we'd like to support some kind of customization.
    let headerFlag = true;
    let diff = numberOfColumns - parsed.length;
    if (diff === 0) {
        headerFlag = false;
    } else if (diff !== -1) {
        throw "number of matrix columns is not equal to the number of rows in '" + annotation_file.name + "'";
    }

    let headers;
    if (headerFlag) {
        headers = parsed.shift();
    } else {
        headers = parsed[0]; // whatever, just using the first row. Hope it's unique enough!
    }

    let annotations = {}
    headers.forEach((x, i) => {
        annotations[x] = parsed.map(y => y[i]);
    });
    for (const [k, v] of Object.entries(annotations)) {
        let conv = rutils.promoteToNumber(v);
        if (conv !== null) {
            annotations[k] = conv;
        }
        if (summary) {
            annotations[k] = rutils.summarizeArray(annotations[k], { limit: summaryLimit });
        }
    }

    return annotations;
}

export async function preflight(args) {
    let output = {};

    let dims = null;
    if ("genes" in args || "annotations" in args) {
        let mtx_data = rutils.formatFile(args.mtx, false);
        var is_gz = mtx_data.name.endsWith(".gz");
        let mm = afile.realizeMatrixMarket(mtx_data.content);
        let headers = scran.extractMatrixMarketDimensions(mm, { "compressed": is_gz });
        dims = [headers.rows, headers.columns];
    }

    if ("genes" in args) {
        let gene_data = rutils.formatFile(args.genes, false);
        let gene_info = extract_features(gene_data, dims[0]);

        let split_out = rutils.presplitByFeatureType(gene_info);
        if (split_out !== null) {
            output.genes = split_out.genes;
        } else {
            output.genes = { "RNA": gene_info };
        }
    } else {
        output.genes = null;
    }

    if ("annotations" in args) {
        let anno_data = rutils.formatFile(args.annotations, false);
        output.annotations = extract_annotations(anno_data, dims[1], { summary: true });
    } else {
        output.annotations = null;
    }

    return output;
}

export class Reader {
    #mtx;
    #genes;
    #annotations;

    constructor(args, formatted = false) {
        if (!formatted) {
            this.#mtx = rutils.formatFile(args.mtx, false);

            if ("genes" in args) {
                this.#genes = rutils.formatFile(args.genes, false);
            } else {
                this.#genes = null;
            }

            if ("annotations" in args) {
                this.#annotations = rutils.formatFile(args.annotations, false);
            } else {
                this.#annotations = null;
            }
        } else {
            this.#mtx = args.mtx;

            if ("genes" in args) {
                this.#genes = args.genes;
            } else {
                this.#genes = null;
            }

            if ("annotations" in args) {
                this.#annotations = args.annotations;
            } else {
                this.#annotations = null;
            }
        }
    }

    load() {
        var ext = this.#mtx.name.split('.').pop();
        var is_compressed = (ext == "gz");
        let mm = afile.realizeMatrixMarket(this.#mtx.content);

        let output = { matrix: new scran.MultiMatrix };
        try {
            let loaded = scran.initializeSparseMatrixFromMatrixMarket(mm, { "compressed": is_compressed });
            let out_mat = loaded.matrix;
            let out_ids = loaded.row_ids;
            output.matrix.add("RNA", out_mat);

            let raw_gene_info = null;
            if (this.#genes !== null) {
                raw_gene_info = extract_features(this.#genes, out_mat.numberOfRows());
            }
            let gene_info = rutils.reorganizeGenes(out_mat.numberOfRows(), out_ids, raw_gene_info);

            let split_out = rutils.splitByFeatureType(out_mat, out_ids, gene_info);
            if (split_out !== null) {
                scran.free(out_mat);
                output.matrix = split_out.matrices;
                output.row_ids = split_out.row_ids;
                output.genes = split_out.genes;
            } else {
                output.genes = { RNA: gene_info };
                output.row_ids = { RNA: out_ids };
            }

            if (this.#annotations !== null) {
                output.annotations = extract_annotations(this.#annotations, output.matrix.numberOfColumns());
            } else {
                output.annotations = null;
            }

        } catch (e) {
            scran.free(output.matrix);
            throw e;
        }

        return output;
    }

    format() {
        return "MatrixMarket";
    }

    async serialize(embeddedSaver) {
        let files = [await rutils.standardSerialize(this.#mtx, "mtx", embeddedSaver)];

        if (this.#genes !== null) {
            files.push(await rutils.standardSerialize(this.#genes, "genes", embeddedSaver));
        }

        if (this.#annotations !== null) {
            files.push(await rutils.standardSerialize(this.#annotations, "annotations", embeddedSaver));
        }

        return files;
    }
}

export async function unserialize(values, embeddedLoader) {
    let args = {};

    // This should contain 'mtx', and possibly 'genes' and/or 'annotations'.
    for (const x of values) {
        args[x.type] = await rutils.standardUnserialize(x, embeddedLoader);
    }

    return new Reader(args, true);
}
