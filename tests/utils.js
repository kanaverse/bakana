import * as path from "path";
import * as fs from "fs";
import * as bakana from "../src/index.js";
import * as valkana from "valkana";
import * as scran from "scran.js";
import * as wa from "wasmarrays.js";

export async function initializeAll() {
    await bakana.initialize({ localFile: true });
    await valkana.initialize({ localFile: true });
}

export function validateState(path, embedded = true) {
    valkana.validateState(path, embedded, bakana.kanaFormatVersion);
}

export function baseParams() {
    let output = bakana.analysisDefaults();

    // Cut down on the work.
    output.rna_pca.num_pcs = 10;

    // Avoid getting held up by pointless iterations.
    output.tsne.iterations = 10;
    output.umap.num_epochs = 10;

    // Actually do something.
    output.cell_labelling = {
        mouse_references: [ "ImmGen" ],
        human_references: [ "BlueprintEncode" ]
    };
    return output;
}

bakana.setCellLabellingDownload(url => {
    let fpath = path.basename(url);
    return fs.readFileSync("files/references/" + fpath).slice(); // Mimic a buffer from fetch().
});

export function mockOffsets(paths) {
    let offsets = {};
    let sofar = 0;
    for (const p of paths) {
        offsets[sofar] = p;
        sofar += fs.statSync(p).size;
    }
    return offsets;
}

function is_same(left, right) {
    if (left.length != right.length) {
        throw new Error("mismatch in array lengths");
    }
    for (var i = 0; i < left.length; i++) {
        if (left[i] != right[i]) {
            return false;
        }
    }
    return true;
}

export function checkReorganization(matrix, ids, names, loadedMatrix, loadedIds, loadedNames, { mustDiffer = true, referenceSubset = false } = {}) {
    if (!referenceSubset && ids !== null) { 
        throw new Error("reference matrix should not be reorganized");
    } else if (referenceSubset && ids === null) {
        throw new Error("subsetted reference matrix should be reorganized");
    }
    if (loadedIds === null) {
        throw new Error("loaded matrix should be reorganized");
    }

    let NR = matrix.numberOfRows();
    let NC = matrix.numberOfColumns();
    if (loadedMatrix.numberOfRows() != NR || NC != loadedMatrix.numberOfColumns()) {
        throw new Error("loaded and reference matrix have different dimensions");
    }

    if (mustDiffer) {
        let same = true;
        for (var i = 0; i < loadedIds.length; i++) {
            if (loadedIds[i] != i) {
                same = false;
                break;
            }
        }
        if (same) {
            throw new Error("identities should not be trivial after reorganization");
        }
    }

    let sorted_ids = loadedIds.slice().sort((a, b) => a - b);
    let permuter;
    if (referenceSubset) {
        // Assume reference IDs are already sorted when referenceSubset = true.
        if (!is_same(ids, sorted_ids)) {
            throw new Error("reference and loaded identities should have identical elements");
        }

        // Creating a mapping to permute the reference to the loaded order.
        let mapping = {};
        ids.forEach((x, i) => {
            mapping[x] = i;
        });
        permuter = Array.from(loadedIds).map(x => mapping[x]);
    } else {
        // Should contain everything at least once.
        for (var i = 0; i < sorted_ids.length; i++) {
            if (sorted_ids[i] != i) {
                throw new error("loaded identities should contain all consecutive integers");
            }
        }
        permuter = Array.from(loadedIds);
    }

    // Checking that the reorganization matches up with the reference for every 100th column.
    let at_least_one_difference = false;

    for (var c = 0; c < NC; c += 100) {
        let loaded_first = loadedMatrix.column(c);
        let reference_first = matrix.column(c);
        if (mustDiffer && !is_same(loaded_first, reference_first)) {
            at_least_one_difference = true;
        }

        let converted = new reference_first.constructor(NR);
        permuter.forEach((x, i) => {
            converted[i] = reference_first[x];
        });
        if (!is_same(loaded_first, converted)) {
            throw new Error("loaded matrix's reorganization is not consistent"); 
        }
    }

    if (mustDiffer && !at_least_one_difference) {
        throw new Error("loaded and reference columns should have some kind of difference");
    }

    // Now checking the identifiers.
    if (mustDiffer && is_same(names, loadedNames)) {
        throw new Error("loaded and reference names should not have been the same");
    }

    let converted = new names.constructor(NR);
    permuter.forEach((x, i) => {
        converted[i] = names[x];
    });

    if (!is_same(converted, loadedNames)) {
        throw new Error("loaded matrix's name reorganization is not consistent");
    }
}

export async function checkStateResultsMinimal(state) {
    // Inputs:
    {
        let counts = state.inputs.fetchCountMatrix();
        expect(counts instanceof scran.MultiMatrix).toBe(true);
        let feat_anno = state.inputs.fetchFeatureAnnotations();
        let row_ids = state.inputs.fetchRowIds();

        for (const a of counts.available()){ 
            let NR = counts.get(a).numberOfRows();
            expect(a in feat_anno).toBe(true);
            expect(feat_anno[a].numberOfRows()).toBe(NR);
            expect(a in row_ids).toBe(true);
            expect(row_ids[a].length).toBe(NR);
        }
    }

    let nfiltered = state.cell_filtering.fetchFilteredMatrix().numberOfColumns();

    // Cell filtering:
    {
        let counts = state.inputs.fetchCountMatrix();
        let discards = state.cell_filtering.fetchDiscards();
        if (discards !== null) {
            expect(discards.length).toBe(counts.numberOfColumns());
            expect(discards instanceof wa.Uint8WasmArray).toBe(true);
        }

        let last_filtered = nfiltered - 1;
        let idx = [0, last_filtered];
        state.cell_filtering.undoFilter(idx);
        expect(idx[1]).toBeGreaterThan(last_filtered);

        let filtered = state.cell_filtering.fetchFilteredMatrix();
        expect(counts.available()).toEqual(filtered.available());

        for (const a of counts.available()) {
            let acounts = counts.get(a);
            let afiltered = filtered.get(a);
            expect(acounts.numberOfRows()).toEqual(afiltered.numberOfRows());
            expect(acounts.column(idx[0])).toEqual(afiltered.column(0));
            expect(acounts.column(idx[1])).toEqual(afiltered.column(last_filtered));
        }
    }

    // Combine embeddings:
    {
        let nc = state.combine_embeddings.fetchNumberOfCells();
        expect(nc).toEqual(nfiltered);

        let nd = state.combine_embeddings.fetchNumberOfDimensions();
        let com = state.combine_embeddings.fetchCombined();
        expect(com.length).toEqual(nc * nd);
    }

    // Batch correction:
    {
        let nc = state.batch_correction.fetchNumberOfCells();
        expect(nc).toEqual(state.combine_embeddings.fetchNumberOfCells());

        let nd = state.batch_correction.fetchNumberOfDimensions();
        expect(nd).toEqual(state.combine_embeddings.fetchNumberOfDimensions());

        let com = state.batch_correction.fetchCorrected();
        expect(com.length).toEqual(nc * nd);
    }

    let nclusters = (new Set(state.choose_clustering.fetchClusters().array())).size;

    // Clustering:
    {
        let clusters = state.choose_clustering.fetchClusters();
        expect(clusters.array() instanceof Int32Array).toBe(true);
        expect(clusters.length).toBe(nfiltered);
        expect(nclusters).toBeGreaterThan(1);
    }

    // t-SNE:
    {
        let coords = await state.tsne.fetchResults();
        expect(coords.x.length).toEqual(nfiltered);
        expect(coords.y.length).toEqual(nfiltered);
        expect(coords.iterations).toBeGreaterThan(0);
    }

    // UMAP:
    {
        let coords = await state.umap.fetchResults();
        expect(coords.x.length).toEqual(nfiltered);
        expect(coords.y.length).toEqual(nfiltered);
        expect(coords.iterations).toBeGreaterThan(0);
    }

    // Markers:
    {
        let counts = state.inputs.fetchCountMatrix();
        let res = state.marker_detection.fetchResults();
        expect(counts.available().sort()).toEqual(Object.keys(res).sort());

        for (const a of counts.available()) {
            let ngenes = counts.get(a).numberOfRows();
            let ares = res[a];
            expect(ares.numberOfGroups()).toEqual(nclusters);

            let res0 = ares.cohen(0);
            expect(res0 instanceof Float64Array).toBe(true);
            expect(res0.length).toEqual(ngenes);

            let resN = ares.auc(nclusters - 1);
            expect(resN instanceof Float64Array).toBe(true);
            expect(resN.length).toEqual(ngenes);
        }
    }

    return;
}

export async function checkStateResultsRna(state, { exclusive = false } = {}) {
    // Inputs:
    {
        let counts = state.inputs.fetchCountMatrix();
        expect(counts.has("RNA")).toBe(true);

        if (exclusive) {
            expect(state.inputs.fetchCountMatrix().has("ADT")).toBe(false);
            expect("ADT" in state.inputs.fetchFeatureAnnotations()).toBe(false);
            expect("ADT" in state.inputs.fetchRowIds()).toBe(false);
        }
    }

    // Quality control:
    {
        let metres = state.rna_quality_control.fetchMetrics();
        expect(metres instanceof scran.PerCellRnaQcMetricsResults).toBe(true);

        let sumvec = metres.sums();
        expect(sumvec instanceof Float64Array).toBe(true);
        let ncells = state.inputs.fetchCountMatrix().numberOfColumns();
        expect(sumvec.length).toBe(ncells);

        expect(state.rna_quality_control.fetchDiscards().length).toEqual(ncells);

        let filtres = state.rna_quality_control.fetchFilters();
        expect(filtres instanceof scran.SuggestRnaQcFiltersResults).toBe(true);
    }

    let nfiltered = state.cell_filtering.fetchFilteredMatrix().numberOfColumns();

    // Cell filtering:
    {
        let detvec = state.rna_quality_control.fetchMetrics().detected();
        expect(detvec.length).toBeGreaterThan(nfiltered);
        let refiltered = state.cell_filtering.applyFilter(detvec);
        expect(refiltered.length).toEqual(nfiltered);

        if (exclusive) {
            expect(state.cell_filtering.fetchDiscards().owner).toBe(state.rna_quality_control.fetchDiscards()); // i.e. a view
        }

        let last_filtered = nfiltered - 1;
        let idx = [0, last_filtered];
        state.cell_filtering.undoFilter(idx);
        expect(idx[1]).toBeGreaterThan(last_filtered);

        let counts = state.inputs.fetchCountMatrix().get("RNA");
        let filtered = state.cell_filtering.fetchFilteredMatrix().get("RNA");
        expect(counts.column(idx[0])).toEqual(filtered.column(0));
        expect(counts.column(idx[1])).toEqual(filtered.column(last_filtered));
    }

    let ngenes = state.inputs.fetchCountMatrix().get("RNA").numberOfRows();

    // Normalization:
    {
        let normed = state.rna_normalization.fetchNormalizedMatrix();
        expect(normed instanceof scran.ScranMatrix).toBe(true);
        expect(normed.numberOfColumns()).toBe(nfiltered);
        expect(normed.numberOfRows()).toBe(ngenes);

        let sf = state.rna_normalization.fetchSizeFactors();
        expect(sf instanceof wa.Float64WasmArray).toBe(true);
        expect(sf.length).toBe(nfiltered);
    }

    // Feature selection:
    {
        let feats = state.feature_selection.fetchResults();
        expect(feats instanceof scran.ModelGeneVarResults).toBe(true);

        let resids = feats.residuals();
        expect(resids.length).toBe(ngenes);

        let sresids = state.feature_selection.fetchSortedResiduals();
        resids.sort((a, b) => a - b);
        expect(sresids).toEqual(resids);
    }

    // PCA:
    {
        let pcs = state.rna_pca.fetchPCs();
        expect(pcs instanceof scran.RunPCAResults).toBe(true);
        expect(pcs.numberOfCells()).toBe(nfiltered);
        expect(pcs.numberOfPCs()).toBe(state.rna_pca.fetchParameters().num_pcs);
    }

    // Combine embeddings:
    {
        let nd = state.combine_embeddings.fetchNumberOfDimensions();

        if (exclusive) {
            expect(nd).toEqual(state.rna_pca.fetchPCs().numberOfPCs());
            let com = state.combine_embeddings.fetchCombined();
            expect(com.owner).toEqual({}); // a.k.a. it's a view.
        } else {
            expect(nd).toBeGreaterThanOrEqual(state.rna_pca.fetchPCs().numberOfPCs());
        }
    }

    // ADTs are no-ops.
    if (exclusive) {
        expect(state.adt_quality_control.fetchMetrics()).toBeUndefined();
        expect(state.adt_normalization.fetchSizeFactors()).toBeUndefined();
        expect(state.adt_pca.fetchPCs()).toBeUndefined();

        expect(state.rna_quality_control.valid()).toBe(true);
        expect(state.rna_normalization.valid()).toBe(true);
        expect(state.rna_pca.valid()).toBe(true);

        expect(state.adt_quality_control.valid()).toBe(false);
        expect(state.adt_normalization.valid()).toBe(false);
        expect(state.adt_pca.valid()).toBe(false);
    }

    return;
}

export async function checkStateResultsUnblocked(state) {
    // Inputs.
    {
        expect(state.inputs.fetchBlock()).toBeNull();
        expect(state.inputs.fetchBlockLevels()).toBeNull();
    }

    // Quality control.
    if (state.rna_quality_control.valid()) {
        let filtres = state.rna_quality_control.fetchFilters();
        expect(filtres.thresholdsSums().length).toBe(1);
    }
    if (state.adt_quality_control.valid()) {
        let filtres = state.adt_quality_control.fetchFilters();
        expect(filtres.thresholdsDetected().length).toBe(1);
    }

    // Batch correction (no-op as it's a view).
    {
        let com = state.batch_correction.fetchCorrected();
        expect(com.owner).toEqual({});
    }

    // Cell filtering.
    {
        expect(state.cell_filtering.fetchFilteredBlock()).toBeNull();
    }

    // Marker detection has a single block.
    {
        let res = state.marker_detection.fetchResults();
        for (const v of Object.values(res)) {
            expect(v.numberOfBlocks()).toBe(1);
        }
    }
}

export async function checkStateResultsBlocked(state) {
    // Checking the non-NULL blocking factors.
    let nlevels = state.inputs.fetchBlockLevels().length;
    {
        expect(nlevels).toBeGreaterThan(1);

        let ids = state.inputs.fetchBlock();
        let counts = {};
        ids.forEach((x, i) => {
            if (x in counts) {
                counts[x] = 1;
            } else {
                counts[x]++;
            }
        });
        expect(Object.keys(counts).length).toBe(nlevels);
    }

    // Check that multiple QC thresholds exist.
    if (state.rna_quality_control.valid()) {
        let res = state.rna_quality_control.fetchFilters();
        let props = res.thresholdsSubsetProportions(0);
        expect(props.length).toEqual(nlevels);
    }
    if (state.adt_quality_control.valid()) {
        let res = state.adt_quality_control.fetchFilters();
        let props = res.thresholdsSubsetTotals(0);
        expect(props.length).toEqual(nlevels);
    }
    if (state.crispr_quality_control.valid()) {
        let res = state.crispr_quality_control.fetchFilters();
        let props = res.thresholdsMaxCount();
        expect(props.length).toEqual(nlevels);
    }

    // Check that some correction actually occured.
    {
        let corrected_pcs = state.batch_correction.fetchCorrected();
        let original_pcs = state.combine_embeddings.fetchCombined();
        expect(corrected_pcs.length).toEqual(original_pcs.length);
        expect(corrected_pcs.array()).not.toEqual(original_pcs.array());
    }

    // Checking that the marker results show up with multiple blocks.
    {
        let res = state.marker_detection.fetchResults();
        for (const v of Object.values(res)) {
            expect(v.numberOfBlocks()).toEqual(nlevels);
        }
    }

    return;
}

export function checkStateResultsAdt(state, { exclusive = false } = {}) {
    // Checking that ADTs exist.
    {
        expect(state.inputs.fetchCountMatrix().has("ADT")).toBe(true);
        let adtmat = state.inputs.fetchCountMatrix().get("ADT");
        expect(state.inputs.fetchFeatureAnnotations()["ADT"].numberOfRows()).toEqual(adtmat.numberOfRows());
        expect(state.inputs.fetchRowIds()["ADT"].length).toEqual(adtmat.numberOfRows());
    }

    // Checking the QCs.
    {
        let amet = state.adt_quality_control.fetchMetrics();
        expect(amet instanceof scran.PerCellAdtQcMetricsResults).toBe(true);

        let positive_total = 0;
        amet.subsetTotals(0, { copy: false }).forEach(x => { positive_total += (x > 0); });
        expect(positive_total).toBeGreaterThan(0);

        expect(state.adt_quality_control.fetchDiscards().length).toEqual(state.inputs.fetchCountMatrix().numberOfColumns());

        let afilt = state.adt_quality_control.fetchFilters();
        expect(afilt instanceof scran.SuggestAdtQcFiltersResults).toBe(true);
        expect(afilt.thresholdsDetected()[0]).toBeGreaterThan(0);
        expect(afilt.thresholdsSubsetTotals(0)[0]).toBeGreaterThan(0);
    }

    let nfiltered = state.cell_filtering.fetchFilteredMatrix().numberOfColumns();

    // Cell filtering:
    {
        let detvec = state.adt_quality_control.fetchMetrics().detected();
        expect(detvec.length).toBeGreaterThan(nfiltered);
        let refiltered = state.cell_filtering.applyFilter(detvec);
        expect(refiltered.length).toEqual(nfiltered);

        if (exclusive) {
            expect(state.cell_filtering.fetchDiscards().owner).toBe(state.adt_quality_control.fetchDiscards()); // i.e. a view
        }

        let last_filtered = nfiltered - 1;
        let idx = [0, last_filtered];
        state.cell_filtering.undoFilter(idx);
        expect(idx[1]).toBeGreaterThan(last_filtered);

        let counts = state.inputs.fetchCountMatrix().get("ADT");
        let filtered = state.cell_filtering.fetchFilteredMatrix().get("ADT");
        expect(counts.column(idx[0])).toEqual(filtered.column(0));
        expect(counts.column(idx[1])).toEqual(filtered.column(last_filtered));
    }

    // Normalization.
    {
        let sf = state.adt_normalization.fetchSizeFactors();
        expect(sf.length).toEqual(state.cell_filtering.fetchFilteredMatrix().numberOfColumns());

        let positive_total = 0;
        sf.forEach(x => { positive_total += (x > 0); });
        expect(positive_total).toBeGreaterThan(0);

        let norm = state.adt_normalization.fetchNormalizedMatrix();
        expect(norm.numberOfColumns()).toEqual(sf.length);
        expect(norm.numberOfRows()).toEqual(state.inputs.fetchCountMatrix().get("ADT").numberOfRows());
    }

    // PCA.
    {
        let pcs = state.adt_pca.fetchPCs();
        expect(pcs.numberOfCells()).toEqual(nfiltered);
        expect(pcs.numberOfPCs()).toBeGreaterThan(0);
        expect(pcs.numberOfPCs()).toBeLessThan(state.adt_pca.fetchParameters().num_pcs);
    }

    // Combined embeddings.
    {
        let nd = state.combine_embeddings.fetchNumberOfDimensions();

        if (exclusive) {
            expect(nd).toEqual(state.adt_pca.fetchPCs().numberOfPCs());
            let com = state.combine_embeddings.fetchCombined();
            expect(com.owner).toEqual({}); // a.k.a. it's a view.
        } else {
            expect(nd).toBeGreaterThanOrEqual(state.adt_pca.fetchPCs().numberOfPCs());
        }
    }

    // RNA are no-ops.
    if (exclusive) {
        expect(state.rna_quality_control.fetchMetrics()).toBeUndefined();
        expect(state.rna_normalization.fetchSizeFactors()).toBeUndefined();
        expect(state.rna_pca.fetchPCs()).toBeUndefined();

        expect(state.rna_quality_control.valid()).toBe(false);
        expect(state.rna_normalization.valid()).toBe(false);
        expect(state.rna_pca.valid()).toBe(false);

        expect(state.adt_quality_control.valid()).toBe(true);
        expect(state.adt_normalization.valid()).toBe(true);
        expect(state.adt_pca.valid()).toBe(true);
    }

    return;
}

export function checkStateResultsCrispr(state, { exclusive = false, use_embeddings = true } = {}) {
    // Checking that CRISPR exists.
    {
        expect(state.inputs.fetchCountMatrix().has("CRISPR")).toBe(true);
        let crisprmat = state.inputs.fetchCountMatrix().get("CRISPR");
        expect(state.inputs.fetchFeatureAnnotations()["CRISPR"].numberOfRows()).toEqual(crisprmat.numberOfRows());
        expect(state.inputs.fetchRowIds()["CRISPR"].length).toEqual(crisprmat.numberOfRows());
    }

    // Checking the QCs.
    {
        let cmet = state.crispr_quality_control.fetchMetrics();
        expect(cmet instanceof scran.PerCellCrisprQcMetricsResults).toBe(true);

        expect(state.crispr_quality_control.fetchDiscards().length).toEqual(state.inputs.fetchCountMatrix().numberOfColumns());

        let cfilt = state.crispr_quality_control.fetchFilters();
        expect(cfilt instanceof scran.SuggestCrisprQcFiltersResults).toBe(true);
        expect(cfilt.thresholdsMaxCount()[0]).toBeGreaterThan(0);
    }

    let nfiltered = state.cell_filtering.fetchFilteredMatrix().numberOfColumns();

    // Cell filtering:
    {
        let detvec = state.crispr_quality_control.fetchMetrics().detected();
        expect(detvec.length).toBeGreaterThan(nfiltered);
        let refiltered = state.cell_filtering.applyFilter(detvec);
        expect(refiltered.length).toEqual(nfiltered);

        if (exclusive) {
            expect(state.cell_filtering.fetchDiscards().owner).toBe(state.crispr_quality_control.fetchDiscards()); // i.e. a view
        }

        let last_filtered = nfiltered - 1;
        let idx = [0, last_filtered];
        state.cell_filtering.undoFilter(idx);
        expect(idx[1]).toBeGreaterThan(last_filtered);

        let counts = state.inputs.fetchCountMatrix().get("CRISPR");
        let filtered = state.cell_filtering.fetchFilteredMatrix().get("CRISPR");
        expect(counts.column(idx[0])).toEqual(filtered.column(0));
        expect(counts.column(idx[1])).toEqual(filtered.column(last_filtered));
    }

    // Normalization.
    {
        let sf = state.crispr_normalization.fetchSizeFactors();
        expect(sf.length).toEqual(state.cell_filtering.fetchFilteredMatrix().numberOfColumns());

        let positive_total = 0;
        sf.forEach(x => { positive_total += (x > 0); });
        expect(positive_total).toBeGreaterThan(0);

        let norm = state.crispr_normalization.fetchNormalizedMatrix();
        expect(norm.numberOfColumns()).toEqual(sf.length);
        expect(norm.numberOfRows()).toEqual(state.inputs.fetchCountMatrix().get("CRISPR").numberOfRows());
    }

    // PCA.
    {
        let pcs = state.crispr_pca.fetchPCs();
        expect(pcs.numberOfCells()).toEqual(nfiltered);
        expect(pcs.numberOfPCs()).toBeGreaterThan(0);
        expect(pcs.numberOfPCs()).toEqual(state.crispr_pca.fetchParameters().num_pcs);
    }

    // Combined embeddings.
    if (use_embeddings) {
        let nd = state.combine_embeddings.fetchNumberOfDimensions();

        if (exclusive) {
            expect(nd).toEqual(state.crispr_pca.fetchPCs().numberOfPCs());
            let com = state.combine_embeddings.fetchCombined();
            expect(com.owner).toEqual({}); // a.k.a. it's a view.
        } else {
            expect(nd).toBeGreaterThanOrEqual(state.crispr_pca.fetchPCs().numberOfPCs());
        }
    }

    // RNA are no-ops.
    if (exclusive) {
        expect(state.rna_quality_control.fetchMetrics()).toBeUndefined();
        expect(state.rna_normalization.fetchSizeFactors()).toBeUndefined();
        expect(state.rna_pca.fetchPCs()).toBeUndefined();

        expect(state.rna_quality_control.valid()).toBe(false);
        expect(state.rna_normalization.valid()).toBe(false);
        expect(state.rna_pca.valid()).toBe(false);

        expect(state.crispr_quality_control.valid()).toBe(true);
        expect(state.crispr_normalization.valid()).toBe(true);
        expect(state.crispr_pca.valid()).toBe(true);
    }

    return;
}

export function checkStateResultsRnaPlusAdt(state) {
    // Cell filtering responds to both modalities.
    let nfiltered = state.cell_filtering.fetchFilteredMatrix().numberOfColumns();
    {
        expect(state.cell_filtering.fetchDiscards().owner).toBeNull(); // i.e. not a view.

        let rna_only = 0;
        state.rna_quality_control.fetchDiscards().forEach(x => { rna_only += (x > 0); });

        let adt_only = 0;
        state.adt_quality_control.fetchDiscards().forEach(x => { adt_only += (x > 0); });

        expect(nfiltered).toBeGreaterThan(rna_only);
        expect(nfiltered).toBeGreaterThan(adt_only);
    }

    // Combined embeddings.
    {
        let rna_dims = state.rna_pca.fetchPCs().numberOfPCs();
        let adt_dims = state.adt_pca.fetchPCs().numberOfPCs();
        expect(state.combine_embeddings.fetchNumberOfDimensions()).toEqual(rna_dims + adt_dims);
        expect(state.combine_embeddings.fetchNumberOfCells()).toEqual(nfiltered);
    }

    return;
}

export function checkStateResultsRnaPlusCrispr(state, { use_embeddings = true } = {}) {
    // Cell filtering responds to both modalities.
    let nfiltered = state.cell_filtering.fetchFilteredMatrix().numberOfColumns();
    {
        expect(state.cell_filtering.fetchDiscards().owner).toBeNull(); // i.e. not a view.

        let rna_only = 0;
        state.rna_quality_control.fetchDiscards().forEach(x => { rna_only += (x > 0); });

        let crispr_only = 0;
        state.crispr_quality_control.fetchDiscards().forEach(x => { crispr_only += (x > 0); });

        expect(nfiltered).toBeGreaterThan(rna_only);
        expect(nfiltered).toBeGreaterThan(crispr_only);
    }

    // Combined embeddings.
    if (use_embeddings) {
        let rna_dims = state.rna_pca.fetchPCs().numberOfPCs();
        let crispr_dims = state.crispr_pca.fetchPCs().numberOfPCs();
        expect(state.combine_embeddings.fetchNumberOfDimensions()).toEqual(rna_dims + crispr_dims);
        expect(state.combine_embeddings.fetchNumberOfCells()).toEqual(nfiltered);
    }

    return;
}

export async function overlordCheckStandard(state) {
    await checkStateResultsMinimal(state);
    await checkStateResultsRna(state, { exclusive: true });
    await checkStateResultsUnblocked(state);
}

export async function overlordCheckBlocked(state) {
    await checkStateResultsMinimal(state);
    await checkStateResultsRna(state, { exclusive: true });
    await checkStateResultsBlocked(state);
}

export async function compareStates(left, right, { checkRna = true, checkAdt = false, checkCrispr = false } = {}) {
    // Parameter checks.
    for (const step of Object.keys(left)) {
        let params = right[step].fetchParameters();
        let ref = left[step].fetchParameters();
        expect(ref).toEqual(params);
    }

    // Inputs:
    {
        let lcounts = left.inputs.fetchCountMatrix();
        let rcounts = right.inputs.fetchCountMatrix();

        let lavailable = lcounts.available();
        expect(lavailable).toEqual(rcounts.available());
        expect(lcounts.numberOfColumns()).toEqual(rcounts.numberOfColumns());

        for (const a of lavailable) {
            let lmat = lcounts.get(a);
            let rmat = rcounts.get(a);
            let NR = lmat.numberOfRows();
            expect(NR).toEqual(rmat.numberOfRows());
            expect(lmat.row(0)).toEqual(rmat.row(0));
            expect(lmat.row(NR-1)).toEqual(rmat.row(NR-1));
        }

        // Checking that the permutation is unchanged on reload.
        for (const a of lavailable) {
            let old_ids = left.inputs.fetchRowIds()[a];
            let new_ids = right.inputs.fetchRowIds()[a];
            expect(old_ids.length).toBeGreaterThan(0);
            expect(old_ids).toEqual(new_ids);

            let old_ids2 = left.inputs.fetchFeatureAnnotations()[a];
            let new_ids2 = right.inputs.fetchFeatureAnnotations()[a];
            expect(old_ids2.columnNames()).toEqual(new_ids2.columnNames());
            for (const col of old_ids2.columnNames()) {
                expect(old_ids2.column(col)).toEqual(new_ids2.column(col));
            }
        }
    }

    // Quality control:
    if (checkRna) {
        let lmetrics = left.rna_quality_control.fetchMetrics();
        let rmetrics = right.rna_quality_control.fetchMetrics();
        expect(lmetrics.sums()).toEqual(rmetrics.sums());
        expect(lmetrics.detected()).toEqual(rmetrics.detected());
        expect(lmetrics.subsetProportions(0)).toEqual(rmetrics.subsetProportions(0));

        let lfilters = left.rna_quality_control.fetchFilters();
        let rfilters = right.rna_quality_control.fetchFilters();
        expect(lfilters.thresholdsSums()).toEqual(rfilters.thresholdsSums());
        expect(lfilters.thresholdsSubsetProportions(0)).toEqual(rfilters.thresholdsSubsetProportions(0));

        let ldiscards = left.rna_quality_control.fetchDiscards().array();
        let rdiscards = right.rna_quality_control.fetchDiscards().array();
        expect(ldiscards).toEqual(rdiscards);
    }

    if (checkAdt) {
        let lmetrics = left.adt_quality_control.fetchMetrics();
        let rmetrics = right.adt_quality_control.fetchMetrics();
        expect(lmetrics.sums()).toEqual(rmetrics.sums());
        expect(lmetrics.detected()).toEqual(rmetrics.detected());
        expect(lmetrics.subsetTotals(0)).toEqual(rmetrics.subsetTotals(0));

        let lfilters = left.adt_quality_control.fetchFilters();
        let rfilters = right.adt_quality_control.fetchFilters();
        expect(lfilters.thresholdsDetected()).toEqual(rfilters.thresholdsDetected());
        expect(lfilters.thresholdsSubsetTotals(0)).toEqual(rfilters.thresholdsSubsetTotals(0));

        let ldiscards = left.adt_quality_control.fetchDiscards().array();
        let rdiscards = right.adt_quality_control.fetchDiscards().array();
        expect(ldiscards).toEqual(rdiscards);
    }

    if (checkCrispr) {
        let lmetrics = left.crispr_quality_control.fetchMetrics();
        let rmetrics = right.crispr_quality_control.fetchMetrics();
        expect(lmetrics.sums()).toEqual(rmetrics.sums());
        expect(lmetrics.detected()).toEqual(rmetrics.detected());
        expect(lmetrics.maxProportions()).toEqual(rmetrics.maxProportions());

        let lfilters = left.crispr_quality_control.fetchFilters();
        let rfilters = right.crispr_quality_control.fetchFilters();
        expect(lfilters.thresholdsMaxCount()).toEqual(rfilters.thresholdsMaxCount());

        let ldiscards = left.crispr_quality_control.fetchDiscards().array();
        let rdiscards = right.crispr_quality_control.fetchDiscards().array();
        expect(ldiscards).toEqual(rdiscards);
    }

    // Cell filtering:
    {
        let ldiscard = left.cell_filtering.fetchDiscards();
        let rdiscard = right.cell_filtering.fetchDiscards();
        if (ldiscard == null) {
            expect(rdiscard).toBeNull();
        } else {
            expect(ldiscard.array()).toEqual(rdiscard.array());
        }

        let lfiltered = left.cell_filtering.fetchFilteredMatrix();
        let rfiltered = right.cell_filtering.fetchFilteredMatrix();
        expect(lfiltered.available()).toEqual(rfiltered.available());

        for (const a of lfiltered.available()) {
            let lmat = lfiltered.get(a);
            let rmat = rfiltered.get(a);
            let NR = lmat.numberOfRows();
            expect(NR).toEqual(rmat.numberOfRows());
            expect(lmat.row(0)).toEqual(rmat.row(0));
            expect(lmat.row(NR-1)).toEqual(rmat.row(NR-1));
        }
    }

    // Normalization:
    if (checkRna) {
        let lmat = left.rna_normalization.fetchNormalizedMatrix();
        let rmat = right.rna_normalization.fetchNormalizedMatrix();

        let NR = lmat.numberOfRows();
        expect(NR).toEqual(rmat.numberOfRows());
        expect(lmat.row(0)).toEqual(rmat.row(0));
        expect(lmat.row(NR-1)).toEqual(rmat.row(NR-1));

        let lsf = left.rna_normalization.fetchSizeFactors();
        let rsf = right.rna_normalization.fetchSizeFactors();
        expect(lsf.array()).toEqual(rsf.array());
    }

    if (checkAdt) {
        let lmat = left.adt_normalization.fetchNormalizedMatrix();
        let rmat = right.adt_normalization.fetchNormalizedMatrix();

        let NR = lmat.numberOfRows();
        expect(NR).toEqual(rmat.numberOfRows());
        expect(lmat.row(0)).toEqual(rmat.row(0));
        expect(lmat.row(NR-1)).toEqual(rmat.row(NR-1));

        let lsf = left.adt_normalization.fetchSizeFactors();
        let rsf = right.adt_normalization.fetchSizeFactors();
        expect(lsf.array()).toEqual(rsf.array());
    }

    // Feature selection:
    if (checkRna) {
        let old_fres = left.feature_selection.fetchResults();
        let new_fres = right.feature_selection.fetchResults();
        for (const x of [ "means", "residuals" ]) {
            expect(old_fres[x]()).toEqual(new_fres[x]());
        }

        let old_sresids = left.feature_selection.fetchSortedResiduals();
        let new_sresids = right.feature_selection.fetchSortedResiduals();
        expect(old_sresids).toEqual(new_sresids);
    }

    // PCA:
    if (checkRna) {
        let lpcs = left.rna_pca.fetchPCs();
        let rpcs = right.rna_pca.fetchPCs();
        expect(lpcs.principalComponents()).toEqual(rpcs.principalComponents());

        let lvp = lpcs.varianceExplained();
        lvp.forEach((x, i) => { lvp[i] /= lpcs.totalVariance(); });
        let rvp = rpcs.varianceExplained();
        rvp.forEach((x, i) => { rvp[i] /= rpcs.totalVariance(); });

        expect(lvp).toEqual(rvp);
    }

    if (checkAdt) {
        let lpcs = left.adt_pca.fetchPCs();
        let rpcs = right.adt_pca.fetchPCs();
        expect(lpcs.principalComponents()).toEqual(rpcs.principalComponents());

        let lvp = lpcs.varianceExplained();
        lvp.forEach((x, i) => { lvp[i] /= lpcs.totalVariance(); });
        let rvp = rpcs.varianceExplained();
        rvp.forEach((x, i) => { rvp[i] /= rpcs.totalVariance(); });

        expect(lvp).toEqual(rvp);
    }

    if (checkCrispr) {
        let lpcs = left.crispr_pca.fetchPCs();
        let rpcs = right.crispr_pca.fetchPCs();
        expect(lpcs.principalComponents()).toEqual(rpcs.principalComponents());

        let lvp = lpcs.varianceExplained();
        lvp.forEach((x, i) => { lvp[i] /= lpcs.totalVariance(); });
        let rvp = rpcs.varianceExplained();
        rvp.forEach((x, i) => { rvp[i] /= rpcs.totalVariance(); });

        expect(lvp).toEqual(rvp);
    }

    // Combine embeddings:
    {
        expect(left.combine_embeddings.fetchCombined().array()).toEqual(right.combine_embeddings.fetchCombined().array());
        expect(left.combine_embeddings.fetchNumberOfCells()).toEqual(right.combine_embeddings.fetchNumberOfCells());
        expect(left.combine_embeddings.fetchNumberOfDimensions()).toEqual(right.combine_embeddings.fetchNumberOfDimensions());
    }

    // Batch correction:
    {
        expect(left.batch_correction.fetchCorrected().array()).toEqual(right.batch_correction.fetchCorrected().array());
        expect(left.batch_correction.fetchNumberOfCells()).toEqual(right.batch_correction.fetchNumberOfCells());
        expect(left.batch_correction.fetchNumberOfDimensions()).toEqual(right.batch_correction.fetchNumberOfDimensions());
    }

    // Clustering:
    {
        expect(left.choose_clustering.fetchClusters().array()).toEqual(right.choose_clustering.fetchClusters().array());
    }

    // t-SNE:
    {
        let lres = await left.tsne.fetchResults();
        let rres = await right.tsne.fetchResults();
        expect(rres.x.length).toEqual(lres.x.length); // direct comparison fails for reasons I don't understand.
        expect(rres.x.buffer).toEqual(lres.x.buffer);
        expect(rres.y.length).toEqual(lres.y.length);
        expect(rres.y.buffer).toEqual(lres.y.buffer);
    }

    // UMAP:
    {
        let lres = await left.umap.fetchResults();
        let rres = await right.umap.fetchResults();
        expect(rres.x.length).toEqual(lres.x.length);
        expect(rres.x.buffer).toEqual(lres.x.buffer);
        expect(rres.y.length).toEqual(lres.y.length);
        expect(rres.y.buffer).toEqual(lres.y.buffer);
    }

    // Markers:
    {
        let lres = left.marker_detection.fetchResults();
        let rres = right.marker_detection.fetchResults();
        let available = Object.keys(lres).sort();
        expect(Object.keys(rres).sort()).toEqual(available);

        for (const a of available) {
            let ng = lres[a].numberOfGroups();
            expect(ng).toEqual(rres[a].numberOfGroups());

            for (var g = 0; g < ng; g++) {
                expect(lres[a].cohen(g)).toEqual(rres[a].cohen(g));
                expect(lres[a].auc(g)).toEqual(rres[a].auc(g));
                expect(lres[a].means(g)).toEqual(rres[a].means(g));
                expect(lres[a].detected(g)).toEqual(rres[a].detected(g));
            }
        }
    }

    return;
}

export function checkClusterVersusMode(state) {
    let nclusters = (new Set(state.choose_clustering.fetchClusters().array())).size;

    let vres = state.marker_detection.computeVersus(nclusters - 1, 0);
    expect(vres.results.RNA instanceof scran.ScoreMarkersResults).toBe(true);
    let lfcs = vres.results.RNA.lfc(vres.left);

    let vres2 = state.marker_detection.computeVersus(0, nclusters - 1);
    expect(vres2.results.RNA instanceof scran.ScoreMarkersResults).toBe(true);
    expect(vres2.left).toBe(vres.right);

    let lfcs2 = vres2.results.RNA.lfc(vres2.left);
    lfcs2.forEach((x, i) => {
        lfcs2[i] *= -1;
    });
    expect(lfcs).toEqual(lfcs2);

    return vres;
}

export function launchCustomSelections(state) {
    let ncells = state.cell_filtering.fetchFilteredMatrix().numberOfColumns();
    state.custom_selections.addSelection("first", [0,1,2,3,4]);
    state.custom_selections.addSelection("last", [ncells - 5, ncells - 4, ncells - 3, ncells - 2, ncells - 1]);

    return {
        first: state.custom_selections.fetchResults("first"),
        last: state.custom_selections.fetchResults("last", "cohen", "RNA"),
        versus: state.custom_selections.computeVersus("last", "first")
    };
}

export async function triggerAnimation(state) {
    let collected = { tsne: [], umap: [] }
    let fun = bakana.setVisualizationAnimate((type, x, y, iterations) => {
        collected[type].push(iterations);
    });

    let p = [state.tsne.animate(), state.umap.animate()];
    await Promise.all(p);
    bakana.setVisualizationAnimate(null)

    expect(collected.tsne.length).toBeGreaterThan(0);
    expect(collected.umap.length).toBeGreaterThan(0);
    return;
}
