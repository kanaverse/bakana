import * as bioc from "bioconductor";
import * as df from "./DataFrame.js";
import { MockSparseMatrix, MockNormalizedMatrix } from "./assays.js";
import { MockReducedDimensionMatrix } from "./reducedDimensions.js";

export async function formatSingleCellExperiment(state, path, { reportOneIndex = false, storeModalityColumnData = false } = {}) {
    let all_rowdata = state.inputs.fetchFeatureAnnotations();
    let modalities = Object.keys(row_info);
    let main = "RNA";
    if (!(main in row_info)) {
        main = modalities[0];
    }

    let all_metadata = {};
    for (const mod of modalities) {
        all_metadata[m] = new bioc.List;
    }

    let all_coldata = df.formatColumnData(
        state,
        modalities,
        main,
        all_metadata,
        storeModalityColumnData
    );

    // Setting up the SEs.
    let all_se = {}
    for (const m of modalities) {
        let mat = state.cell_filtering.fetchFilteredMatrix().get(m);
        let assays = { counts: new MockSparseMatrix(mat) };
        let reddim = {}; 
        let metadata;
        if (m in all_metadata) {
            metadata = all_metadata[m];
        } else {
            metadata = new bioc.List;
        }

        let step = null;
        switch (m) {
            case "RNA":
                step = state.rna_normalization;
                break;
            case "ADT":
                step = state.adt_normalization;
                break;
            case "CRISPR":
                step = state.crispr_normalization;
                break;
        }
        if (step !== null) {
            let sf = step.fetchSizeFactors();
            assays.logcounts = new MockNormalizedMatrix(mat, sf);
        }

        step = null;
        switch (m) {
            case "RNA":
                step = state.rna_pca;
                break;
            case "ADT":
                step = state.adt_pca;
                break;
            case "CRISPR":
                step = state.crispr_pca;
                break;
        }
        if (step !== null) {
            let pcs = step.fetchPCs();
            reddim.pca = new MockReducedDimensionMatrix(pcs.numberOfCells(), pcs.numberOfPCs(), pcs.principalComponents({ copy: "view" }));
            let tv = pcs.totalVariance();
            metadata.pca = { variance_explained: pcs.varianceExplained({ copy: "view" }).map(x => x / tv) };
        }

        all_se[m] = new bioc.SingleCellExperiment(
            assays,
            {
                rowData: all_rowdata[m],
                columnData: all_coldata[m],
                metadata: metadata,
                reducedDimensions: reddim
            }
        );
    }

    // Saving the dimensionality reduction results.
    for (const name of [ "tsne", "umap" ]) {
        let res = await state[name].fetchResults({ copy: false });
        let payload = new Float64Array(res.x.length * 2);
        payload.set(res.x);
        payload.set(res.y, res.x.length);
        all_se[main].setReducedDimension(name, new MockReducedDimensionMatrix(res.x.length, 2, payload), { inPlace: true });
    }

    // Saving extra metadata.
    {
        let customs = state.custom_selections.fetchSelections({ copy: true, force: "Int32Array" });
        if (reportOneIndex) {
            for (const [k, v] of Object.entries(customs)) { // incrementing by 1, if requested.
                v.forEach((x, i) => { v[i] = x + 1; });
            }
        }
        let meta = all_se[main].metadata();
        meta.set("custom_selections", customs, { inPlace: true });
        all_se[main].setMetadata(meta, { inPlace: true });
    }

    let main_se = all_se[main];
    for (const mod of modalities) {
        if (mod !== main) {
            main_se.setAlternativeExperiment(mod, all_se[mod], { inPlace: true });
        }
    }

    return main_se;
}
