import * as bioc from "bioconductor";
import * as df from "./DataFrame.js";

const translate_effects = { "lfc": "lfc", "delta_detected": "deltaDetected", "auc": "auc", "cohen": "cohen" };

export function formatMarkerDetectionResults(state, modality_names) {
    const translate_summary = { "min": 0, "mean": 1, "min_rank": 4 };
    const do_auc = state.marker_detection.fetchParameters().compute_auc;
    let all_rowdata = state.inputs.fetchFeatureAnnotations();
    let all_output = {};

    for (const [m, rn] of Object.entries(modality_names)) {
        let res = state.marker_detection.fetchResults()[m];
        let ngroups = res.numberOfGroups();
        let nfeatures = all_rowdata[m].numberOfRows();

        for (var group = 0; group < ngroups; group++) {
            let mdf = new bioc.DataFrame({}, { numberOfRows: nfeatures, rowNames: rn });

            for (const x of [ "means", "detected" ]) {
                mdf.setColumn(x, res[x](group, { copy: "view" }), { inPlace: true });
            }

            for (const [eff, trans_eff] of Object.entries(translate_effects)) {
                if (eff == "auc" && !do_auc) {
                    continue;
                }
                for (const [summ, trans_summ] of Object.entries(translate_summary)) {
                    mdf.setColumn(eff + "-" + summ, res[trans_eff](group, { summary: trans_summ }), { inPlace: true });
                }
            }

            let new_name = group + 1; // incrementing to avoid cluster names starting from 0.
            all_output[m + "/" + String(new_name)] = mdf;
        }
    }

    return all_output;
}

export function dumpCustomSelectionResults(state, modality_names) {
    const do_auc = state.custom_selections.fetchParameters().compute_auc;
    let all_sel = state.custom_selections.fetchSelections();
    let all_rowdata = state.inputs.fetchFeatureAnnotations();
    let all_output = {};

    for (const [m, rn] of Object.entries(modality_names)) {
        let nfeatures = all_rowdata[m].numberOfRows();

        for (const sel of Object.keys(all_sel)) {
            let res = state.custom_selections.fetchResults(sel)[m];
            let mdf = new bioc.DataFrame({}, { numberOfRows: nfeatures, rowNames: rn });

            for (const x of [ "means", "detected" ]) {
                mdf.setColumn(x, res[x](1, { copy: "view" }), { inPlace: true });
            }

            for (const [eff, trans_eff] of Object.entries(translate_effects)) {
                if (eff == "auc" && !do_auc) {
                    continue;
                }
                mdf.setColumn(eff, res[trans_eff](1, { copy: "view" }), { inPlace: true });
            }

            all_output[m + "/" + self] = mdf;
        }
    }

    return all_output;
}

export function formatFeatureSelectionResults(state, rna_names) {
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

    return fdf;
}
