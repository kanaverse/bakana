import * as scran from "scran.js";

export function reorganizeGenes(nrows, ids, geneInfo) {
    if (geneInfo === null) {
        let genes = [];
        if (ids !== null) {
            for (const i of ids) {
                genes.push(`Gene ${i + 1}`);
            }
        } else {
            for (let i = 0; i < nrows; i++) {
                genes.push(`Gene ${i + 1}`);
            }
        }
        geneInfo = { "id": genes };
    } else {
        if (ids != null ){
            for (const [k, v] of Object.entries(geneInfo)) {
                geneInfo[k] = scran.quickSliceArray(ids, v);
            }
        }
    }
    return geneInfo;
}

function rawSplitByFeatureType(genes) { 
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

    let output = { types: types };

    // Skipping 'type', as it's done its purpose now.
    let gene_deets = { ...genes };
    delete gene_deets.type;
    output.genes = scran.splitArrayCollection(gene_deets, types);

    return output;
}

export function presplitByFeatureType(genes) { 
    let output = rawSplitByFeatureType(genes);
    if (output === null) {
        return null;
    }
    delete output.types;
    return output;
}

export function splitByFeatureType(matrix, row_ids, genes) { 
    let output = rawSplitByFeatureType(genes);
    if (output === null) {
        return null;
    }

    let types = output.types;
    delete output.types;

    // Allocating the split matrices. 
    let out_mats;
    try {
        out_mats = new scran.MultiMatrix({ store: scran.splitRows(matrix, types) });
        output.matrices = out_mats;
    } catch (e) {
        scran.free(out_mats);
        throw e;
    }

    // Also saving the fragmented ids.
    if (row_ids === null) {
        output.row_ids = types;
    } else {
        output.row_ids = {};
        for (const [k, v] of Object.entries(types)) {
            output.row_ids[k] = scran.quickSliceArray(v, row_ids);
        }
    }

    return output;
}
