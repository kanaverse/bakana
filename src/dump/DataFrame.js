import * as wa from "wasmarrays.js";
import * as bioc from "bioconductor";

// Monkey-patching these methods so that we can use these WasmArrays as columns in a bioc.DataFrame.
wa.Uint8WasmArray.prototype._bioconductor_LENGTH = function() { return this.length; };
wa.Int32WasmArray.prototype._bioconductor_LENGTH = function() { return this.length; };
wa.Float64WasmArray.prototype._bioconductor_LENGTH = function() { return this.length; };

export function formatColumnData(state, all_modalities, main_modality, all_other_metadata, store_per_modality) {
    let keep_raw = state.cell_filtering.fetchKeep();
    let keep = [];
    if (keep_raw !== null) {
        keep_raw.forEach((x, i) => {
            if (x) { keep.push(i); }
        });
    }

    let all_coldata = {};
    {
        let full = state.inputs.fetchCellAnnotations();
        all_coldata[main_modality] = (keep_raw === null ? full : bioc.SLICE(full, keep));

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
        rdf = rdf.setColumn(prefix + "::sums", state.cell_filtering.applyFilter(state.rna_quality_control.fetchMetrics().sum({ copy: false })));
        rdf = rdf.setColumn(prefix + "::detected", state.cell_filtering.applyFilter(state.rna_quality_control.fetchMetrics().detected({ copy: false })));
        rdf = rdf.setColumn(prefix + "::proportions", state.cell_filtering.applyFilter(state.rna_quality_control.fetchMetrics().subsetProportion(0, { copy: false })));
        all_coldata[target] = rdf;

        all_other_metadata[target].set(
            prefix,
            { 
                "filters": {
                    "sums": state.rna_quality_control.fetchFilters().sum(),
                    "detected": state.rna_quality_control.fetchFilters().detected(),
                    "proportions": state.rna_quality_control.fetchFilters().subsetProportion(0)
                }
            },
            { inPlace: true }
        );
    }

    if (state.adt_quality_control.valid()) {
        let prefix = (store_per_modality ? "kana::quality_control" : "kana::ADT::quality_control");
        let target = (store_per_modality ? "ADT" : main_modality);

        let adf = all_coldata[target];
        adf = adf.setColumn(prefix + "::sums", state.cell_filtering.applyFilter(state.adt_quality_control.fetchMetrics().sum({ copy: false })));
        adf = adf.setColumn(prefix + "::detected", state.cell_filtering.applyFilter(state.adt_quality_control.fetchMetrics().detected({ copy: false })));
        adf = adf.setColumn(prefix + "::igg_totals", state.cell_filtering.applyFilter(state.adt_quality_control.fetchMetrics().subsetSum(0, { copy: false })));
        all_coldata[target] = adf;

        all_other_metadata[target].set(
            prefix,
            {
                "filters": {
                    "detected": state.adt_quality_control.fetchFilters().detected(),
                    "igg_totals": state.adt_quality_control.fetchFilters().subsetSum(0)
                }
            },
            { inPlace: true }
        );
    }

    if (state.crispr_quality_control.valid()) {
        let prefix = (store_per_modality ? "kana::quality_control" : "kana::CRISPR::quality_control");
        let target = (store_per_modality ? "CRISPR" : main_modality);

        let cdf = all_coldata[target];
        cdf = cdf.setColumn(prefix + "::sums", state.cell_filtering.applyFilter(state.crispr_quality_control.fetchMetrics().sum({ copy: false })));
        cdf = cdf.setColumn(prefix + "::detected", state.cell_filtering.applyFilter(state.crispr_quality_control.fetchMetrics().detected({ copy: false })));
        cdf = cdf.setColumn(prefix + "::max_proportion", state.cell_filtering.applyFilter(state.crispr_quality_control.fetchMetrics().maxProportion({ copy: false })));
        cdf = cdf.setColumn(prefix + "::max_index", state.cell_filtering.applyFilter(state.crispr_quality_control.fetchMetrics().maxIndex({ copy: false })));
        all_coldata[target] = cdf;

        all_other_metadata[target].set(
            prefix,
            {
                "filters": {
                    "max_count": state.crispr_quality_control.fetchFilters().maxValue()
                }
            },
            { inPlace: true }
        );
    }

    if (keep_raw !== null) {
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
