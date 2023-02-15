import * as bioc from "bioconductor";
import * as writedf from "./writeDataFrame.js";
import * as reddim from "./dumpReducedDimensions.js";

const translate_effects = { "lfc": "lfc", "delta_detected": "deltaDetected", "auc": "auc", "cohen": "cohen" };

export function dumpToSingleCellExperiment(state, { forceArrayBuffer = false } = {}) {
    let row_info = state.inputs.fetchFeatureAnnotations();
    let modalities = Object.keys(row_info);
    let main = "RNA";
    if (!(main in row_info)) {
        main = modalities[0];
    }

    let all_files = [];
    let all_meta = {};
    for (const m of modalities) {
        all_meta[m] = {
            "$schema": "single_cell_experiment/v1.json",
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
                "main_experiment_name": m,
                "alternative_experiments": [],
                "reduced_dimensions": []
            }
        };
    }

    let altpath = m => "altexp-" + m;
    for (const m of modalities) {
        if (m !== main) {
            all_meta[main].single_cell_experiment.alternative_experiments.push({ name: m, resource: { type: "local", path: altpath(m) } });
        }
    }

    let all_metadata = {};
    for (const m of modalities) {
        all_metadata[m] = {};
    }

    /**** Assembling the column data for each modality. ****/
    {
        let keep = [];
        state.cell_filtering.fetchDiscards().forEach((x, i) => {
            if (!x) { keep.push(i); }
        });

        let all_coldata = {};
        for (const m of modalities) {
            all_coldata[m] = new bioc.DataFrame({}, { numberOfRows: keep.length });
        }

        if ("RNA" in row_info) {
            all_coldata.RNA.$setColumn("qc_sums", state.cell_filtering.applyFilter(state.rna_quality_control.fetchMetrics().sums({ copy: false })));
            all_coldata.RNA.$setColumn("qc_detected", state.cell_filtering.applyFilter(state.rna_quality_control.fetchMetrics().detected({ copy: false })));
            all_coldata.RNA.$setColumn("qc_proportions", state.cell_filtering.applyFilter(state.rna_quality_control.fetchMetrics().subsetProportions(0, { copy: false })));

            all_metadata.RNA["qc_sums_filter"] = state.rna_quality_control.fetchFilters().thresholdsSums();
            all_metadata.RNA["qc_detected_filter"] = state.rna_quality_control.fetchFilters().thresholdsDetected();
            all_metadata.RNA["qc_proportions_filter"] = state.rna_quality_control.fetchFilters().thresholdsSubsetProportions(0);
        }

        if ("ADT" in row_info) {
            all_coldata.ADT.$setColumn("qc_sums", state.cell_filtering.applyFilter(state.adt_quality_control.fetchMetrics().sums({ copy: false })));
            all_coldata.ADT.$setColumn("qc_detected", state.cell_filtering.applyFilter(state.adt_quality_control.fetchMetrics().detected({ copy: false })));
            all_coldata.ADT.$setColumn("qc_igg_totals", state.cell_filtering.applyFilter(state.adt_quality_control.fetchMetrics().subsetTotals(0, { copy: false })));

            all_metadata.ADT["qc_detected_filter"] = state.adt_quality_control.fetchFilters().thresholdsDetected();
            all_metadata.ADT["qc_igg_totals_filter"] = state.adt_quality_control.fetchFilters().thresholdsSubsetTotals(0);
        }

        if ("CRISPR" in row_info) {
            all_coldata.CRISPR.$setColumn("qc_sums", state.cell_filtering.applyFilter(state.adt_quality_control.fetchMetrics().sums({ copy: false })));
            all_coldata.CRISPR.$setColumn("qc_detected", state.cell_filtering.applyFilter(state.adt_quality_control.fetchMetrics().detected({ copy: false })));
            all_coldata.CRISPR.$setColumn("qc_max_proportion", state.cell_filtering.applyFilter(state.adt_quality_control.fetchMetrics().maxProportion({ copy: false })));
            all_coldata.CRISPR.$setColumn("qc_max_index", state.cell_filtering.applyFilter(state.adt_quality_control.fetchMetrics().maxIndex({ copy: false })));

            all_metadata.CRISPR["qc_max_count_filter"] = state.adt_quality_control.fetchFilters().thresholdsMaxCount();
        }

        all_coldata[main].$setColumn("qc_retained", new Int32Array(keep));

        let block = state.cell_filtering.fetchFilteredBlock();
        if (block !== null) {
            let levels = state.inputs.fetchBlockLevels();
            let realized = new Array(block.length);
            block.forEach((x, i) => { realized[i] = levels[x] });
            all_coldata[main].$setColumn("block", realized);
        }

        if ("ADT" in row_info) {
            all_coldata.ADT.$setColumn("sizeFactor", state.adt_normalization.fetchSizeFactors());
        }

        all_coldata[main].$setColumn("clusters", state.choose_clustering.fetchClusters());

        for (const m of modalities) {
            let target = (m == main ? "" : altpath(m) + "/");
            let collected = writedf.writeHdf5DataFrame(all_coldata[m], target + "coldata", { forceArrayBuffer });
            all_meta[m].summarized_experiment.column_data.resource.path = collected.self.metadata.path;
            all_files.push(collected.self);
            for (const x of collected.children) {
                all_files.push(x);
            }
        }
    }

    /**** Assembling the row data for each modality. ****/
    {
        let all_rowdata = {};
        for (const m of modalities) {
            all_rowdata[m] = bioc.CLONE(row_info[m], { deepCopy: false });
        }

        if ("RNA" in row_info) {
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

        // Storing marker detection results.
        {
            const translate_summary = { "min": 0, "mean": 1, "min_rank": 4 };
            const do_auc = state.marker_detection.fetchParameters().compute_auc;

            for (const m of modalities) {
                let res = state.marker_detection.fetchResults()[m];
                let ngroups = res.numberOfGroups();
                let nfeatures = all_rowdata[m].numberOfRows();
                let rdf = new bioc.DataFrame({}, { numberOfRows: nfeatures });

                for (var group = 0; group < ngroups; group++) {
                    let mdf = new bioc.DataFrame({}, { numberOfRows: nfeatures });

                    for (const x of [ "means", "detected" ]) {
                        mdf.$setColumn(x, res[x](group, { copy: "view" }));
                    }

                    for (const [eff, trans_eff] of Object.entries(translate_effects)) {
                        if (eff == "auc" && !do_auc) {
                            continue;
                        }
                        for (const [summ, trans_summ] of Object.entries(translate_summary)) {
                            mdf.$setColumn(eff + "-" + summ, res[trans_eff](group, { summary: trans_summ }));
                        }
                    }

                    rdf.$setColumn(String(group), mdf);
                }

                all_rowdata[m].$setColumn("marker_detection", rdf);
            }
        }

        // Storing custom selection marker results.
        {
            const translate_summary = { "min": 0, "mean": 1, "min_rank": 4 };
            const do_auc = state.custom_selections.fetchParameters().compute_auc;

            let all_sel = state.custom_selections.fetchSelections();
            all_metadata[main] = { custom_selections: all_sel };

            for (const m of modalities) {
                let nfeatures = all_rowdata[m].numberOfRows();
                let rdf = new bioc.DataFrame({}, { numberOfRows: nfeatures });

                for (const sel of Object.keys(all_sel)) {
                    let res = state.custom_selections.fetchResults(sel, m);
                    let mdf = new bioc.DataFrame({}, { numberOfRows: nfeatures });

                    for (const x of [ "means", "detected" ]) {
                        mdf.$setColumn(x, res[x](1, { copy: "view" }));
                    }

                    for (const [eff, trans_eff] of Object.entries(translate_effects)) {
                        if (eff == "auc" && !do_auc) {
                            continue;
                        }
                        mdf.$setColumn(eff, res[trans_eff](1, { copy: "view" }));
                    }

                    rdf.$setColumn(String(group), mdf);
                }
                all_rowdata[m].$setColumn("custom_selections", rdf);
            }
        }

        for (const m of modalities) {
            let target = (m == main ? "" : altpath(m) + "/");
            let collected = writedf.writeHdf5DataFrame(all_rowdata[m], target + "rowdata", { forceArrayBuffer });
            all_meta[m].summarized_experiment.row_data.resource.path = collected.self.metadata.path;
            all_files.push(collected.self);
            console.log(collected.children);
            for (const x of collected.children) {
                all_files.push(x);
            }
        }
    }

    // Saving the count assay.
    {
        for (const m of modalities) {
            let mat = state.cell_filtering.fetchFilteredMatrix().get(m);
            let target = (m == main ? "" : altpath(m) + "/");
            let meta = {
                "$schema": "hdf5_sparse_matrix/v1.json",
                "path": target + "assay-counts/matrix.h5",
                "array": {
                    "dimensions": [mat.numberOfRows(), mat.numberOfColumns()]
                },
                "hdf5_sparse_matrix": {
                    "group": "matrix",
                    "format": "tenx_matrix"
                }
            };

            let temppath = scran.chooseTemporaryPath({ extension: ".h5" });
            let contents = temppath;
            try {
                scran.writeSparseMatrixToHdf5(mat, temppath, "matrix", { format: "tenx_matrix" });
                if (forceArrayBuffer) {
                    contents = scran.readFile(temppath);
                    scran.removeFile(temppath);
                }
            } catch (e) {
                scran.removeFile(temppath);
                throw e;
            }

            all_meta[m].summarized_experiment.assays.push({ name: "counts", resource: { type: "local", path: meta.path } });
            all_files.push({ metadata: meta, contents: contents });
        }
    }

    // Saving the dimensionality reduction results.
    {
        let dump_pca = (m, pcs) => {
            let saved = reddim.dump_pca_to_hdf5(pcs, forceArrayBuffer);
            saved.metadata.path = altpath(m) + "reddim-pca/matrix.h5"
            all_files.push(saved);

            let tv = pcs.totalVariance();
            all_metadata[m].pca_variance_explained = pcs.varianceExplained({ copy: "view" }).map(x => x / tv);
            all_meta[m].single_cell_experiment.reduced_dimensions.push({ name: "PCA", { resource: { type: "local", path: saved.metadata.path } } });
        };

        if ("RNA" in row_info) {
            dump_pca("RNA", state.rna_pca.fetchPCs());
        }

        if ("ADT" in row_info) {
            dump_pca("ADT", state.adt_pca.fetchPCs());
        }

        if ("CRISPR" in row_info) {
            dump_pca("CRISPR", state.crispr_pca.fetchPCs());
        }

        {
            let res = await state.tsne.fetchResults({ copy: false });
            let saved = reddim.dumpOtherReducedDimensionsToHdf5([ res.x, res.y ], forceArrayBuffer);
            saved.metadata.path = "reddim-tsne/matrix.h5"
            all_meta[main].single_cell_experiment.reduced_dimensions.push({ name: "TSNE", { resource: { type: "local", path: saved.metadata.path } } });
        }

        {
            let res = await state.umap.fetchResults({ copy: false });
            let saved = reddim.dumpOtherReducedDimensionsToHdf5([ res.x, res.y ], forceArrayBuffer);
            saved.metadata.path = "reddim-umap/matrix.h5"
            all_meta[main].single_cell_experiment.reduced_dimensions.push({ name: "UMAP", { resource: { type: "local", path: saved.metadata.path } } });
        }
    }

    // Saving the metadata (including cluster labelling assignments).
    {
    }

    for (const m of modalities) {
        if (m !== main) {
            all_files.push({ metadata: all_meta[m] });
        }
    }

    return { metadata: all_meta[main], files: all_files };
}
