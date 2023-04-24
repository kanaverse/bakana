import * as bioc from "bioconductor";
import * as df from "./DataFrame.js";
import * as assay from "./assays.js";
import * as reddim from "./reducedDimensions.js";
import * as list from "./List.js";

/************************************************
 ************************************************/

function dumpColumnData(state, modality_prefixes, main_modality, all_sce_metadata, all_other_metadata, all_files, forceBuffer, store_per_modality) {
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
        for (const k of Object.keys(modality_prefixes)) {
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

        all_other_metadata[target][prefix] = { 
            "filters": {
                "sums": state.rna_quality_control.fetchFilters().thresholdsSums(),
                "detected": state.rna_quality_control.fetchFilters().thresholdsDetected(),
                "proportions": state.rna_quality_control.fetchFilters().thresholdsSubsetProportions(0)
            }
        };
    }

    if (state.adt_quality_control.valid()) {
        let prefix = (store_per_modality ? "kana::quality_control" : "kana::ADT::quality_control");
        let target = (store_per_modality ? "ADT" : main_modality);

        let adf = all_coldata[target];
        adf = adf.setColumn(prefix + "::sums", state.cell_filtering.applyFilter(state.adt_quality_control.fetchMetrics().sums({ copy: false })));
        adf = adf.setColumn(prefix + "::detected", state.cell_filtering.applyFilter(state.adt_quality_control.fetchMetrics().detected({ copy: false })));
        adf = adf.setColumn(prefix + "::igg_totals", state.cell_filtering.applyFilter(state.adt_quality_control.fetchMetrics().subsetTotals(0, { copy: false })));
        all_coldata[target] = adf;

        all_other_metadata[target][prefix] = {
            "filters": {
                "detected": state.adt_quality_control.fetchFilters().thresholdsDetected(),
                "igg_totals": state.adt_quality_control.fetchFilters().thresholdsSubsetTotals(0)
            }
        };
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

        all_other_metadata[target][prefix] = {
            "filters": {
                "max_count": state.crispr_quality_control.fetchFilters().thresholdsMaxCount()
            }
        };
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

    // Dumping everything to file.
    for (const [name, prefix] of Object.entries(modality_prefixes)) {
        let collected = df.writeHdf5DataFrame(all_coldata[name], prefix + "coldata", { forceBuffer });
        collected.self.metadata.is_child = true;
        all_sce_metadata[name].summarized_experiment.column_data.resource.path = collected.self.metadata.path;
        all_files.push(collected.self);
        for (const x of collected.children) {
            all_files.push(x);
        }
    }

    return;
}

/************************************************
 ************************************************/

function mockSingleCellExperimentMetadata(p) {
    return {
        "$schema": "single_cell_experiment/v1.json",
        "path": p,
        "summarized_experiment": {
            "row_data": {
                "resource": {
                    "type": "local",
                    "path": null
                }
            },
            "column_data": {
                "resource": {
                    "type": "local",
                    "path": null
                }
            },
            "assays": [],
            "other_data": {
                "resource": {
                    "type": "local",
                    "path": null
                }
            }
        },
        "single_cell_experiment": {
            "alternative_experiments": [],
            "reduced_dimensions": []
        }
    };
}

/************************************************
 ************************************************/

export async function dumpSingleCellExperiment(state, path, { reportOneIndex = false, storeModalityColumnData = false, forceBuffer = false } = {}) {
    let row_info = state.inputs.fetchFeatureAnnotations();

    let modalities = Object.keys(row_info);
    let main = "RNA";
    if (!(main in row_info)) {
        main = modalities[0];
    } else {
        // Make sure it's always first.
        let replacement = [main];
        for (const m of modalities) {
            if (m !== main) {
                replacement.push(m);
            }
        }
    }

    let all_files = [];
    let all_top_meta = {};
    let all_prefixes = {};
    let all_metadata = {};

    for (const m of modalities) {
        let mprefix = path + "/";
        if (m != main) {
            mprefix += "altexp-" + m + "/";
        }
        all_prefixes[m] = mprefix;

        let sce_path = mprefix + "experiment.json";
        all_top_meta[m] = mockSingleCellExperimentMetadata(sce_path);

        let mat = state.cell_filtering.fetchFilteredMatrix().get(m);
        all_top_meta[m].summarized_experiment.dimensions = [mat.numberOfRows(), mat.numberOfColumns()];

        if (m == main) {
            all_top_meta[m].single_cell_experiment.main_experiment_name = m;
        } else {
            // As main modality is first, it's always guaranteed to exist by the time we want to push to it.
            all_top_meta[main].single_cell_experiment.alternative_experiments.push({ name: m, resource: { type: "local", path: sce_path } });
        }

        all_metadata[m] = {};
    }

    // Saving the column and row data.
    dumpColumnData(state, all_prefixes, main, all_top_meta, all_metadata, all_files, forceBuffer, storeModalityColumnData);

    {
        let row_info = state.inputs.fetchFeatureAnnotations();
        for (const [name, prefix] of Object.entries(all_prefixes)) {
            let collected = df.writeHdf5DataFrame(row_info[name], prefix + "rowdata", { forceBuffer });
            collected.self.metadata.is_child = true;
            all_top_meta[name].summarized_experiment.row_data.resource.path = collected.self.metadata.path;
            all_files.push(collected.self);
            for (const x of collected.children) {
                all_files.push(x);
            }
        }
    }

    // Saving the count assay.
    for (const [m, prefix] of Object.entries(all_prefixes)) {
        let mat = state.cell_filtering.fetchFilteredMatrix().get(m);
        let saved = assay.dumpCountMatrix(mat, prefix + "assay-counts", forceBuffer);
        saved.metadata.is_child = true;

        all_top_meta[m].summarized_experiment.assays.push({ name: "counts", resource: { type: "local", path: saved.metadata.path } });
        all_files.push(saved);
    }

    // Saving the log-normalized assay.
    for (const [m, prefix] of Object.entries(all_prefixes)) {
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
        if (step == null) {
            continue;
        }

        let mat = step.fetchNormalizedMatrix();
        let sf = step.fetchSizeFactors();
        let saved = assay.dumpNormalizedMatrix(mat, sf, prefix + "assay-logcounts", all_top_meta[m].summarized_experiment.assays[0].resource.path, forceBuffer);
        saved.metadata.is_child = true;

        all_top_meta[m].summarized_experiment.assays.push({ name: "logcounts", resource: { type: "local", path: saved.metadata.path } });
        all_files.push(saved);
    }

    // Saving the dimensionality reduction results.
    for (const [m, prefix] of Object.entries(all_prefixes)) {
        let step = null;
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
        if (step == null) {
            continue;
        }

        let pcs = step.fetchPCs();
        let saved = reddim.dumpPcaResultsToHdf5(pcs, prefix + "reddim-pca", forceBuffer);
        saved.metadata.is_child = true;
        all_files.push(saved);

        let tv = pcs.totalVariance();
        all_metadata[m].pca = { variance_explained: pcs.varianceExplained({ copy: "view" }).map(x => x / tv) };
        all_top_meta[m].single_cell_experiment.reduced_dimensions.push({ name: "PCA", resource: { type: "local", path: saved.metadata.path } });
    }

    for (const name of [ "tsne", "umap" ]) {
        let res = await state[name].fetchResults({ copy: false });
        let saved = reddim.dumpOtherReducedDimensionsToHdf5([ res.x, res.y ], all_prefixes[main] + "reddim-" + name, forceBuffer);
        saved.metadata.is_child = true;

        all_top_meta[main].single_cell_experiment.reduced_dimensions.push({ name: name.toUpperCase(), resource: { type: "local", path: saved.metadata.path } });
        all_files.push(saved);
    }

    // Saving extra metadata.
    all_metadata[main].cell_labelling = await state.cell_labelling.fetchResults();

    {
        let customs = state.custom_selections.fetchSelections({ copy: true, force: "Int32Array" });
        if (reportOneIndex) {
            for (const [k, v] of Object.entries(customs)) { // incrementing by 1, if requested.
                v.forEach((x, i) => { v[i] = x + 1; });
            }
        }
        all_metadata[main].custom_selections = customs;
    }

    for (const [m, prefix] of Object.entries(all_prefixes)) {
        let saved = list.dumpList(all_metadata[m], prefix + "other");
        all_top_meta[m].summarized_experiment.other_data.resource.path = saved.metadata.path;
        saved.metadata.is_child = true;
        all_files.push(saved);
    }

    for (const m of modalities) {
        if (m !== main) {
            all_files.push({ metadata: all_top_meta[m] });
        }
    }

    return { metadata: all_top_meta[main], files: all_files };
}
