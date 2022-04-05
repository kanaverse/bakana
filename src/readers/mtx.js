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

function extract_features(gene_file, { numberOfRows = null } = {}) {
    const content = new Uint8Array(gene_file.content.buffer());
    let parsed = rutils.readDSVFromBuffer(content, gene_file);
    if (numberOfRows !== null && parsed.length !== numberOfRows) {
        throw new Error("number of matrix rows is not equal to the number of genes in '" + gene_file.name + "'");
    }

    var ids = [], symb = [];
    parsed.forEach(x => {
        ids.push(x[0]);
        symb.push(x[1]);
    });

    return { "id": ids, "symbol": symb };
}

function extract_annotations(annotation_file, { numberOfColumns = null, namesOnly = false } = {}) {
    const content = new Uint8Array(annotation_file.content.buffer());
    let parsed = rutils.readDSVFromBuffer(content, annotation_file);

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
            this.#genes = args.genes;
            this.#annotations = args.annotations;
        }
    }

    load() {
        var contents = new Uint8Array(this.#mtx.content.buffer());
        var ext = this.#mtx.name.split('.').pop();
        var is_compressed = (ext == "gz");

        let output = {};
        try {
            output.matrix = scran.initializeSparseMatrixFromMatrixMarketBuffer(contents, { "compressed": is_compressed });

            if (this.#genes !== null) {
                output.genes = extract_features(this.#genes, { numberOfRows: output.matrix.numberOfRows() });
            } else {
                output.genes = null;
            }

            if (this.#annotations !== null) {
                output.annotations = extract_annotations(this.#annotations, { numberOfColumns: output.matrix.numberOfColumns() });
            } else {
                output.annotations = null;
            }
        } catch (e) {
            utils.freeCache(output.matrix);
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
