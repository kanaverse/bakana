import * as scran from "scran.js";
import * as utils from "./../utils/general.js";
import * as rutils from "./../utils/reader.js";
import * as afile from "./../abstract/file.js";

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

function extract_features(gene_file, { numberOfRows = null, includeFeatureType = false } = {}) {
    const content = new Uint8Array(gene_file.content.buffer());
    var is_gz = gene_file.name.endsWith(".gz");
    let parsed = rutils.readTable(content, { compression: (is_gz ? "gz" : "none") });
    if (numberOfRows !== null && parsed.length !== numberOfRows) {
        throw new Error("number of matrix rows is not equal to the number of genes in '" + gene_file.name + "'");
    }

    var ids = [], symb = [];
    parsed.forEach(x => {
        ids.push(x[0]);
        symb.push(x[1]);
    });

    let output = { "id": ids, "symbol": symb };

    if (includeFeatureType && parsed[0].length >= 3) {
        let types = [];
        parsed.forEach(x => { types.push(x[2]); });
        output.type = types;
    }

    return output;
}

function extract_annotations(annotation_file, { numberOfColumns = null, namesOnly = false } = {}) {
    const content = new Uint8Array(annotation_file.content.buffer());
    var is_gz = annotation_file.name.endsWith(".gz");
    let parsed = rutils.readTable(content, { compression: (is_gz ? "gz" : "none"), firstOnly: namesOnly });

    // Check if a header is present or not
    let headerFlag = true;
    if (numberOfColumns !== null) {
        let diff = numberOfColumns - parsed.length;
        if (diff === 0) {
            headerFlag = false;
        } else if (diff !== -1) {
            throw "number of annotations rows is not equal to the number of cells in '" + annotation_file.name + "'";
        }
    }

    let headers;
    if (headerFlag) {
        headers = parsed.shift();
    } else {
        headers = parsed[0]; // whatever, just using the first row. Hope they're unique enough!
    }

    if (namesOnly) {
        return headers;
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
    }

    return annotations;
}

export function preflight(args) {
    let output = {};

    if ("genes" in args) {
        let gene_data = rutils.formatFile(args.genes, false);
        output.genes = extract_features(gene_data);
    } else {
        output.genes = null;
    }

    if ("annotations" in args) {
        let anno_data = rutils.formatFile(args.annotations, false);
        output.annotations = extract_annotations(anno_data, { namesOnly: true });
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
        var contents = new Uint8Array(this.#mtx.content.buffer());
        var ext = this.#mtx.name.split('.').pop();
        var is_compressed = (ext == "gz");

        let matrices = new rutils.MultiMatrix;
        let output;
        try {
            let out_mat = scran.initializeSparseMatrixFromMatrixMarketBuffer(contents, { "compressed": is_compressed });
            matrices.add("RNA", out_mat);

            let gene_info = null;
            if (this.#genes !== null) {
                gene_info = extract_features(this.#genes, { numberOfRows: out_mat.numberOfRows(), includeFeatureType: true });
            }
            let genes = { RNA: rutils.reorganizeGenes(out_mat, gene_info) };

            let split_out = rutils.splitByFeatureType(out_mat, genes.RNA);
            if (split_out !== null) {
                utils.freeCache(out_mat);
                matrices = split_out.matrices;
                genes = split_out.genes;
            }

            output = {
                matrix: matrices,
                genes: genes
            };

            if (this.#annotations !== null) {
                let first = matrices.get(matrices.available()[0]);
                output.annotations = extract_annotations(this.#annotations, { numberOfColumns: first.numberOfColumns() });
            }

        } catch (e) {
            utils.freeCache(matrices);
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
