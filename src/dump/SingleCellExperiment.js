import * as bioc from "bioconductor";
import * as df from "./DataFrame.js";

/************************************************
 ************************************************/

function formatColumnData(state, all_modalities, main_modality, all_other_metadata, store_per_modality) {
    let disc = state.cell_filtering.fetchDiscards();
    let keep = [];
    if (disc !== null) {
        disc.forEach((x, i) => {
            if (!x) { keep.push(i); }
        });
    }

    let all_coldata = {};
    {
        let full = state.inputs.fetchCellAnnotations();
        all_coldata[main_modality] = (disc === null ? full : bioc.SLICE(full, keep));

        let nrows = all_coldata[main_modality].numberOfRows();
        for (const k of all_modalities) {
            if (k !== main_modality) {
                all_coldata[k] = new bioc.DataFrame({}, { numberOfRows: nrows });
            }
        }
    }

    // Quality control.
    if (state.rna_quality_control.valid()) {
        let prefix = (store_per_modality ? "kana::quality_control" : "kana::RNA::quality_control");
        let target = (store_per_modality ? "RNA" : main_modality);

        let rdf = all_coldata[target];
        rdf = rdf.setColumn(prefix + "::sums", state.cell_filtering.applyFilter(state.rna_quality_control.fetchMetrics().sums({ copy: false })));
        rdf = rdf.setColumn(prefix + "::detected", state.cell_filtering.applyFilter(state.rna_quality_control.fetchMetrics().detected({ copy: false })));
        rdf = rdf.setColumn(prefix + "::proportions", state.cell_filtering.applyFilter(state.rna_quality_control.fetchMetrics().subsetProportions(0, { copy: false })));
        all_coldata[target] = rdf;

        all_other_metadata[target].set(
            prefix,
            { 
                "filters": {
                    "sums": state.rna_quality_control.fetchFilters().thresholdsSums(),
                    "detected": state.rna_quality_control.fetchFilters().thresholdsDetected(),
                    "proportions": state.rna_quality_control.fetchFilters().thresholdsSubsetProportions(0)
                }
            },
            { inPlace: true }
        );
    }

    if (state.adt_quality_control.valid()) {
        let prefix = (store_per_modality ? "kana::quality_control" : "kana::ADT::quality_control");
        let target = (store_per_modality ? "ADT" : main_modality);

        let adf = all_coldata[target];
        adf = adf.setColumn(prefix + "::sums", state.cell_filtering.applyFilter(state.adt_quality_control.fetchMetrics().sums({ copy: false })));
        adf = adf.setColumn(prefix + "::detected", state.cell_filtering.applyFilter(state.adt_quality_control.fetchMetrics().detected({ copy: false })));
        adf = adf.setColumn(prefix + "::igg_totals", state.cell_filtering.applyFilter(state.adt_quality_control.fetchMetrics().subsetTotals(0, { copy: false })));
        all_coldata[target] = adf;

        all_other_metadata[target].set(
            prefix,
            {
                "filters": {
                    "detected": state.adt_quality_control.fetchFilters().thresholdsDetected(),
                    "igg_totals": state.adt_quality_control.fetchFilters().thresholdsSubsetTotals(0)
                }
            },
            { inPlace: true }
        );
    }

    if (state.crispr_quality_control.valid()) {
        let prefix = (store_per_modality ? "kana::quality_control" : "kana::CRISPR::quality_control");
        let target = (store_per_modality ? "CRISPR" : main_modality);

        let cdf = all_coldata[target];
        cdf = cdf.setColumn(prefix + "::sums", state.cell_filtering.applyFilter(state.crispr_quality_control.fetchMetrics().sums({ copy: false })));
        cdf = cdf.setColumn(prefix + "::detected", state.cell_filtering.applyFilter(state.crispr_quality_control.fetchMetrics().detected({ copy: false })));
        cdf = cdf.setColumn(prefix + "::max_proportion", state.cell_filtering.applyFilter(state.crispr_quality_control.fetchMetrics().maxProportions({ copy: false })));
        cdf = cdf.setColumn(prefix + "::max_index", state.cell_filtering.applyFilter(state.crispr_quality_control.fetchMetrics().maxIndex({ copy: false })));
        all_coldata[target] = cdf;

        all_other_metadata[target].set(
            prefix,
            {
                "filters": {
                    "max_count": state.crispr_quality_control.fetchFilters().thresholdsMaxCount()
                }
            },
            { inPlace: true }
        );
    }

    if (disc !== null) {
        all_coldata[main_modality] = all_coldata[main_modality].setColumn("kana::quality_control::retained_indices", keep);
    }

    // Size Factors.
    if (state.rna_normalization.valid()) {
        let field = (store_per_modality ? "kana::size_factors" : "kana::RNA::size_factors");
        let target = (store_per_modality ? "RNA" : main_modality);
        all_coldata[target] = all_coldata[target].setColumn(field, state.rna_normalization.fetchSizeFactors());
    }

    if (state.adt_normalization.valid()) {
        let field = (store_per_modality ? "kana::size_factors" : "kana::ADT::size_factors");
        let target = (store_per_modality ? "ADT" : main_modality);
        all_coldata[target] = all_coldata[target].setColumn(field, state.adt_normalization.fetchSizeFactors());
    }

    if (state.crispr_normalization.valid()) {
        let field = (store_per_modality ? "kana::size_factors" : "kana::CRISPR::size_factors");
        let target = (store_per_modality ? "CRISPR" : main_modality);
        all_coldata[target] = all_coldata[target].setColumn(field, state.crispr_normalization.fetchSizeFactors());
    }

    {
        // Incrementing to avoid cluster names starting from 0. Note that there's
        // no need to respect the reportOneIndex setting, as cluster names are
        // not indices with respect to anything. The only thing they need to match
        // with is the marker table names, and we increment there (in markers.js) as well.
        let clusters = state.choose_clustering.fetchClusters().map(x => x + 1);
        all_coldata[main_modality] = all_coldata[main_modality].setColumn("kana::clusters", clusters);
    }

    {
        let block = state.cell_filtering.fetchFilteredBlock();
        if (block !== null) {
            let stringy = new Array(block.length);
            let levels = state.inputs.fetchBlockLevels();
            block.forEach((x, i) => { stringy[i] = levels[x]; }); 
            all_coldata[main_modality] = all_coldata[main_modality].setColumn("kana::block", stringy);
        }
    }

    // Custom selections, stored as boolean arrays.
    {
        let customs = state.custom_selections.fetchSelections({ copy: false });
        let nrows = all_coldata[main_modality].numberOfRows();
        for (const [v, k] of Object.entries(customs)) {
            let as_bool = new Uint8Array(nrows);
            as_bool.fill(0);
            k.forEach(index => { as_bool[index] = 1; });
            all_coldata[main_modality] = all_coldata[main_modality].setColumn("kana::custom_selections::" + v, as_bool);
        }
    }

    return all_coldata;
}

/************************************************
 ************************************************/

class MockSparseMatrix {
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

class MockNormalizedMatrix {
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

class MockReducedDimensionMatrix {
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

/************************************************
 ************************************************/

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

    let all_coldata = formatColumnData(
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
