import * as scran from "scran.js";

/**
 * Reorganize an array collection based on the permutations to the identities in a ScranMatrix.
 *
 * @param {ScranMatrix} matrix - A ScranMatrix containing matrix data, possibly with reorganization of the rows. 
 * @param {?object} geneInfo - An object containing arrays of length equal to the number of rows in `matrix`.
 * Each array contains some kind of per-gene annotation, where each element corresponds to a row in the original matrix prior to any reorganization.
 *
 * @return If `geneInfo` is not `null`, it is modified in place by reordering the values in each array so that they match the identities of rows in `matrix`.
 * A reference to `geneInfo` is then returned.
 *
 * If `geneInfo` is `null`, an object is returned containing an `id` field with an array of placeholder identifiers for all rows.
 * Each placeholder is named after the position of the row in the original matrix.
 */
export function reorganizeGenes(matrix, geneInfo) {
    if (geneInfo === null) {
        let genes = [];
        if (matrix.isReorganized()) {
            let ids = matrix.identities();
            for (const i of ids) {
                genes.push(`Gene ${i + 1}`);
            }
        } else {
            for (let i = 0; i < matrix.numberOfRows(); i++) {
                genes.push(`Gene ${i + 1}`);
            }
        }
        geneInfo = { "id": genes };
    } else {
        if (matrix.isReorganized()) {
            scran.matchFeatureAnnotationToRowIdentities(matrix, geneInfo);
        }
    }
    return geneInfo;
}

export function splitByFeatureType(matrix, genes) { 
    if (!("type" in genes)) {
        return null;
    }

    let types0 = scran.splitByFactor(genes.type);
    if (Object.keys(types0).length == 1) {
        return null;
    }

    // Standardizing the names to something the rest of the pipeline
    // will recognize. By default, we check the 10X vocabulary here.
    let types = {};
    for (const [k, v] of Object.entries(types0)) {
        if (k.match(/gene expression/i)) {
            types["RNA"] = v;
        } else if (k.match(/antibody capture/i)) {
            types["ADT"] = v;
        }
    }

    let output = {};

    // Skipping 'type', as it's done its purpose now.
    let gene_deets = { ...genes };
    delete gene_deets.type;
    output.genes = scran.splitArrayCollection(gene_deets, types);

    if (matrix !== null) {
        // Allocating the split matrices. Note that this is skipped in the
        // 'null' case to support feature splitting for the preflight requests
        // (where the matrix is not loaded, obviously).
        let out_mats;
        try {
            out_mats = new scran.MultiMatrix({ store: scran.splitRows(matrix, types) });
            output.matrices = out_mats;
        } catch (e) {
            scran.safeFree(out_mats);
            throw e;
        }
    }

    return output;
}
