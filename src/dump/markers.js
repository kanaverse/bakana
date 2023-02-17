import * as bioc from "bioconductor";

const translate_effects = { "lfc": "lfc", "delta_detected": "deltaDetected", "auc": "auc", "cohen": "cohen" };

export function dumpMarkerDetectionResults(state, modalities, all_rowdata) {
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

    return;
}

export function dumpCustomSelectionResults(state, modalities, main, all_rowdata, all_other_metadata) {
    const do_auc = state.custom_selections.fetchParameters().compute_auc;

    let all_sel = state.custom_selections.fetchSelections();
    all_other_metadata[main] = { custom_selections: all_sel };

    for (const m of modalities) {
        let nfeatures = all_rowdata[m].numberOfRows();
        let rdf = new bioc.DataFrame({}, { numberOfRows: nfeatures });

        for (const sel of Object.keys(all_sel)) {
            let res = state.custom_selections.fetchResults(sel)[m];
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

            rdf.$setColumn(sel, mdf);
        }
        all_rowdata[m].$setColumn("custom_selections", rdf);
    }
}
