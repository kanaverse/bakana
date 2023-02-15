import * as bioc from "bioconductor";
import * as df from "./DataFrame.js";
import * as assay from "./assays.js";
import * as markers from "./markers.js";
import * as reddim from "./reducedDimensions.js";
import * as list from "./List.js";

/************************************************
 ************************************************/

function dumpColumnData(state, modality_prefixes, main_modality, all_sce_metadata, all_other_metadata, all_files, forceBuffer) {
    let keep = [];
    state.cell_filtering.fetchDiscards().forEach((x, i) => {
        if (!x) { keep.push(i); }
    });
    let retained = keep.length;

    let all_coldata = {};
    for (const m of Object.keys(modality_prefixes)) {
        all_coldata[m] = new bioc.DataFrame({}, { numberOfRows: retained });
    }

    // Quality control.
    if (state.rna_quality_control.valid()) {
        let rdf = new bioc.DataFrame({}, { numberOfRows: retained });
        rdf.$setColumn("sums", state.cell_filtering.applyFilter(state.rna_quality_control.fetchMetrics().sums({ copy: false })));
        rdf.$setColumn("detected", state.cell_filtering.applyFilter(state.rna_quality_control.fetchMetrics().detected({ copy: false })));
        rdf.$setColumn("proportions", state.cell_filtering.applyFilter(state.rna_quality_control.fetchMetrics().subsetProportions(0, { copy: false })));
        all_coldata.RNA.$setColumn("rna_quality_control", rdf);

        all_other_metadata.RNA["rna_quality_control"] = { 
            "filters": {
                "sums": state.rna_quality_control.fetchFilters().thresholdsSums(),
                "detected": state.rna_quality_control.fetchFilters().thresholdsDetected(),
                "proportions": state.rna_quality_control.fetchFilters().thresholdsSubsetProportions(0)
            }
        };
    }

    if (state.adt_quality_control.valid()) {
        let adf = new bioc.DataFrame({}, { numberOfRows: retained });
        adf.$setColumn("sums", state.cell_filtering.applyFilter(state.adt_quality_control.fetchMetrics().sums({ copy: false })));
        adf.$setColumn("detected", state.cell_filtering.applyFilter(state.adt_quality_control.fetchMetrics().detected({ copy: false })));
        adf.$setColumn("igg_totals", state.cell_filtering.applyFilter(state.adt_quality_control.fetchMetrics().subsetTotals(0, { copy: false })));
        all_coldata.ADT.$setColumn("adt_quality_control", adf);

        all_other_metadata.ADT["adt_quality_control"] = {
            "filters": {
                "detected": state.adt_quality_control.fetchFilters().thresholdsDetected(),
                "igg_totals": state.adt_quality_control.fetchFilters().thresholdsSubsetTotals(0)
            }
        };
    }

    if (state.crispr_quality_control.valid()) {
        let cdf = new bioc.DataFrame({}, { numberOfRows: retained });
        cdf.$setColumn("sums", state.cell_filtering.applyFilter(state.crispr_quality_control.fetchMetrics().sums({ copy: false })));
        cdf.$setColumn("detected", state.cell_filtering.applyFilter(state.crispr_quality_control.fetchMetrics().detected({ copy: false })));
        cdf.$setColumn("max_proportion", state.cell_filtering.applyFilter(state.crispr_quality_control.fetchMetrics().maxProportion({ copy: false })));
        cdf.$setColumn("max_index", state.cell_filtering.applyFilter(state.crispr_quality_control.fetchMetrics().maxIndex({ copy: false })));
        all_coldata.CRISPR.$setColumn("crispr_quality_control", cdf);

        all_other_metadata.CRISPR["crispr_quality_control"] = {
            "filters": {
                "max_count": state.crispr_quality_control.fetchFilters().thresholdsMaxCount()
            }
        };
    }

    // Size Factors.
    if (state.adt_normalization.valid()) {
        all_coldata.ADT.$setColumn("sizeFactor", state.adt_normalization.fetchSizeFactors());
    }

    // Other bits and pieces.
    {
        let keepp1 = new Int32Array(keep.length);
        keep.forEach((x, i) => { keepp1[i] = x + 1; }); // 1-based indices.
        all_coldata[main_modality].$setColumn("retained", keepp1);
    }

    {
        let block = state.cell_filtering.fetchFilteredBlock();
        if (block !== null) {
            let realized = new Int32Array(block.length);
            block.forEach((x, i) => { realized[i] = x + 1 }); // 1-based indices.
            all_coldata[main_modality].$setColumn("block", realized);

            let stringy = new Array(block.length);
            let levels = state.inputs.fetchBlockLevels();
            block.forEach((x, i) => { stringy[i] = levels[x]; }); // 1-based indices.
            all_coldata[main_modality].$setColumn("named_block", stringy);

            all_other_metadata[main].block_levels = levels;
        }
    }

    all_coldata[main_modality].$setColumn("clusters", state.choose_clustering.fetchClusters());

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

function dumpRowData(state, modality_prefixes, main_modality, all_sce_metadata, all_other_metadata, all_files, forceBuffer) {
    let row_info = state.inputs.fetchFeatureAnnotations();

    let all_rowdata = {};
    for (const m of Object.keys(modality_prefixes)) {
        all_rowdata[m] = bioc.CLONE(row_info[m], { deepCopy: false });
    }

    if ("RNA" in modality_prefixes) {
        let res = state.feature_selection.fetchResults();
        let df = new bioc.DataFrame({
            mean: res.means({ copy: "view" }),
            variance: res.variances({ copy: "view" }),
            fitted: res.fitted({ copy: "view" }),
            residual: res.residuals({ copy: "view" })
        }, { 
            columnOrder: [ "mean", "variance", "fitted", "residual" ] 
        });
        all_rowdata.RNA.$setColumn("feature_selection", df);
    }

    markers.dumpMarkerDetectionResults(state, Object.keys(modality_prefixes), all_rowdata);

    markers.dumpCustomSelectionResults(state, Object.keys(modality_prefixes), main_modality, all_rowdata, all_other_metadata);

    for (const [name, prefix] of Object.entries(modality_prefixes)) {
        let collected = df.writeHdf5DataFrame(all_rowdata[name], prefix + "rowdata", { forceBuffer });
        collected.self.metadata.is_child = true;
        all_sce_metadata[name].summarized_experiment.row_data.resource.path = collected.self.metadata.path;
        all_files.push(collected.self);
        for (const x of collected.children) {
            all_files.push(x);
        }
    }
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

export async function dumpSingleCellExperiment(state, { forceBuffer = false } = {}) {
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
        let mprefix = (m == main ? "" : "altexp-" + m + "/");
        all_prefixes[m] = mprefix;

        let sce_path = mprefix + "experiment.json";
        all_top_meta[m] = mockSingleCellExperimentMetadata(sce_path);

        if (m == main) {
            all_top_meta[m].single_cell_experiment.main_experiment_name = m;
        } else {
            // As main modality is first, it's always guaranteed to exist by the time we want to push to it.
            all_top_meta[main].single_cell_experiment.alternative_experiments.push({ name: m, resource: { type: "local", path: sce_path } });
        }

        all_metadata[m] = {};
    }

    // Saving the column and row data.
    dumpColumnData(state, all_prefixes, main, all_top_meta, all_metadata, all_files, forceBuffer);

    dumpRowData(state, all_prefixes, main, all_top_meta, all_metadata, all_files, forceBuffer);

    // Saving the count assay.
    for (const [m, prefix] of Object.entries(all_prefixes)) {
        let mat = state.cell_filtering.fetchFilteredMatrix().get(m);
        let saved = assay.dumpCountMatrix(mat, prefix + "assay-counts", forceBuffer);
        saved.metadata.is_child = true;

        all_top_meta[m].summarized_experiment.assays.push({ name: "counts", resource: { type: "local", path: saved.metadata.path } });
        all_files.push(saved);
    }

    // Saving the dimensionality reduction results.
    for (const [m, prefix] of Object.entries(all_prefixes)) {
        let pcs = null;
        switch (m) {
            case "RNA":
                pcs = state.rna_pca.fetchPCs();
                break;
            case "ADT":
                pcs = state.adt_pca.fetchPCs();
                break;
            case "CRISPR":
                pcs = state.crispr_pca.fetchPCs();
                break;
        }

        if (pcs == null) {
            continue;
        }

        let saved = reddim.dumpPcaResultsToHdf5(pcs, prefix + "reddim-pca", forceBuffer);
        saved.metadata.is_child = true;
        all_files.push(saved);

        let tv = pcs.totalVariance();
        all_metadata[m].pca = { variance_explained: pcs.varianceExplained({ copy: "view" }).map(x => x / tv) };
        all_top_meta[m].single_cell_experiment.reduced_dimensions.push({ name: "PCA", resource: { type: "local", path: saved.metadata.path } });
    }

    for (const name of [ "tsne", "umap" ]) {
        let res = await state[name].fetchResults({ copy: false });
        let saved = reddim.dumpOtherReducedDimensionsToHdf5([ res.x, res.y ], "reddim-" + name, forceBuffer);
        saved.metadata.is_child = true;

        all_top_meta[main].single_cell_experiment.reduced_dimensions.push({ name: name.toUpperCase(), resource: { type: "local", path: saved.metadata.path } });
        all_files.push(saved);
    }

    // Saving extra metadata (including cluster labelling assignments).
    all_metadata[main].cell_labelling = state.cell_labelling.fetchResults();

    for (const [m, prefix] of Object.entries(all_prefixes)) {
        let saved = list.dumpList(all_metadata[m], prefix + "other");
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
