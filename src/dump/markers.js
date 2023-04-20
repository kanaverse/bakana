import * as bioc from "bioconductor";
import * as df from "./DataFrame.js";

const translate_effects = { "lfc": "lfc", "delta_detected": "deltaDetected", "auc": "auc", "cohen": "cohen" };

export function dumpMarkerDetectionResults(state, modality_names, prefix, all_files, forceBuffer) {
    const translate_summary = { "min": 0, "mean": 1, "min_rank": 4 };
    const do_auc = state.marker_detection.fetchParameters().compute_auc;
    let all_rowdata = state.inputs.fetchFeatureAnnotations();

    for (const [m, rn] of Object.entries(modality_names)) {
        let res = state.marker_detection.fetchResults()[m];
        let ngroups = res.numberOfGroups();
        let nfeatures = all_rowdata[m].numberOfRows();

        for (var group = 0; group < ngroups; group++) {
            let mdf = new bioc.DataFrame({}, { numberOfRows: nfeatures, rowNames: rn });

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

            let new_name = group + 1; // incrementing to avoid cluster names starting from 0.
            let collected = df.writeHdf5DataFrame(mdf, prefix + "/marker_detection/" + m + "/" + String(new_name), { forceBuffer });
            all_files.push(collected.self);
        }
    }

    return;
}

export function dumpCustomSelectionResults(state, modality_names, prefix, all_files, forceBuffer) {
    const do_auc = state.custom_selections.fetchParameters().compute_auc;
    let all_sel = state.custom_selections.fetchSelections();
    let all_rowdata = state.inputs.fetchFeatureAnnotations();

    for (const [m, rn] of Object.entries(modality_names)) {
        let nfeatures = all_rowdata[m].numberOfRows();

        for (const sel of Object.keys(all_sel)) {
            let res = state.custom_selections.fetchResults(sel)[m];
            let mdf = new bioc.DataFrame({}, { numberOfRows: nfeatures, rowNames: rn });

            for (const x of [ "means", "detected" ]) {
                mdf.$setColumn(x, res[x](1, { copy: "view" }));
            }

            for (const [eff, trans_eff] of Object.entries(translate_effects)) {
                if (eff == "auc" && !do_auc) {
                    continue;
                }
                mdf.$setColumn(eff, res[trans_eff](1, { copy: "view" }));
            }

            let collected = df.writeHdf5DataFrame(mdf, prefix + "/custom_selections/" + m + "/" + sel, { forceBuffer });
            all_files.push(collected.self);
        }
    }

    return;
}

export function dumpFeatureSelectionResults(state, rna_names, prefix, all_files, forceBuffer) {
    let res = state.feature_selection.fetchResults();
    let fdf = new bioc.DataFrame(
        {
            mean: res.means({ copy: "view" }),
            variance: res.variances({ copy: "view" }),
            fitted: res.fitted({ copy: "view" }),
            residual: res.residuals({ copy: "view" })
        }, 
        { 
            columnOrder: [ "mean", "variance", "fitted", "residual" ],
            rowNames: rna_names
        }
    );

    let collected = df.writeHdf5DataFrame(fdf, prefix + "/feature_selection", { forceBuffer });
    all_files.push(collected.self);
    return;
}
