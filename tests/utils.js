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
    output.pca.num_pcs = 10;

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

export async function checkStateResultsBase(state, { mutable = true } = {}) {
    // Inputs:
    {
        let counts = state.inputs.fetchCountMatrix();
        expect(counts instanceof scran.MultiMatrix).toBe(true);
        expect(counts.has("RNA")).toBe(true);

        let feat_anno = state.inputs.fetchFeatureAnnotations();
        expect("RNA" in feat_anno).toBe(true);
        expect(feat_anno.RNA.numberOfRows()).toBe(counts.get("RNA").numberOfRows());

        let row_ids = state.inputs.fetchRowIds();
        expect("RNA" in row_ids).toBe(true);
        expect(row_ids.RNA.length).toBe(counts.get("RNA").numberOfRows());
    }

    // Quality control:
    {
        let metres = state.quality_control.fetchMetrics();
        expect(metres instanceof scran.PerCellQCMetricsResults).toBe(true);

        let sumvec = metres.sums();
        expect(sumvec instanceof Float64Array).toBe(true);
        expect(sumvec.length).toBe(state.inputs.fetchCountMatrix().numberOfColumns());

        let filtres = state.quality_control.fetchFilters();
        expect(filtres instanceof scran.PerCellQCFiltersResults).toBe(true);
    }

    let nfiltered = state.cell_filtering.fetchFilteredMatrix().numberOfColumns();

    // Cell filtering:
    {
        let detvec = state.quality_control.fetchMetrics().detected();
        expect(detvec.length).toBeGreaterThan(nfiltered);
        let refiltered = state.cell_filtering.applyFilter(detvec);
        expect(refiltered.length).toEqual(nfiltered);

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
        let normed = state.normalization.fetchNormalizedMatrix();
        // expect(normed instanceof scran.ScranMatrix).toBe(true);
        expect(normed.numberOfColumns()).toBe(nfiltered);
        expect(normed.numberOfRows()).toBe(ngenes);

        let sf = state.normalization.fetchSizeFactors();
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
        let pcs = state.pca.fetchPCs();
        expect(pcs instanceof scran.RunPCAResults).toBe(true);
        expect(pcs.numberOfCells()).toBe(nfiltered);
        expect(pcs.numberOfPCs()).toBe(state.pca.fetchParameters().num_pcs);
    }

    // Combine embeddings:
    {
        let nc = state.combine_embeddings.fetchNumberOfCells();
        expect(nc).toEqual(nfiltered);

        let nd = state.combine_embeddings.fetchNumberOfDimensions();
        expect(nd).toBeGreaterThanOrEqual(state.pca.fetchPCs().numberOfPCs());

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

    // Animation catcher works correctly.
    {
        let collected = { tsne: [], umap: [] }
        let fun = bakana.setVisualizationAnimate((type, x, y, iterations) => {
            collected[type].push(iterations);
        });

        let p = [state.tsne.animate(), state.umap.animate()];
        await Promise.all(p);
        bakana.setVisualizationAnimate(null)

        expect(collected.tsne.length).toBeGreaterThan(0);
        expect(collected.umap.length).toBeGreaterThan(0);
    }

    // Markers:
    {
        let res = state.marker_detection.fetchResults()["RNA"];
        expect(res.numberOfGroups()).toEqual(nclusters);

        let res0 = res.cohen(0);
        expect(res0 instanceof Float64Array).toBe(true);
        expect(res0.length).toEqual(ngenes);

        let resN = res.auc(nclusters - 1);
        expect(resN instanceof Float64Array).toBe(true);
        expect(resN.length).toEqual(ngenes);
    }

    // Markers in versus mode:
    if (mutable) {
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
    }

    return;
}

export async function checkStateResultsSimple(state, { skipBasic = false } = {}) {
    if (!skipBasic) {
        await checkStateResultsBase(state);
    }

    // Inputs.
    {
        expect(state.inputs.fetchBlock()).toBeNull();
        expect(state.inputs.fetchBlockLevels()).toBeNull();
    }

    // Quality control.
    {
        let filtres = state.quality_control.fetchFilters();
        expect(filtres.thresholdsSums().length).toBe(1);
    }

    // Cell filtering.
    {
        expect(state.cell_filtering.fetchFilteredBlock()).toBeNull();
    }

    // Combined embeddings (no-op).
    {
        let nd = state.combine_embeddings.fetchNumberOfDimensions();
        expect(nd).toEqual(state.pca.fetchPCs().numberOfPCs());
        let com = state.combine_embeddings.fetchCombined();
        expect(com.owner).toEqual({});
    }

    // Batch correction (no-op).
    {
        let com = state.batch_correction.fetchCorrected();
        expect(com.owner).toEqual({});
    }

    // Marker detection has a single block.
    {
        let res = state.marker_detection.fetchResults();
        expect(res.RNA.numberOfBlocks()).toBe(1);
    }

    // ADTs are no-ops.
    {
        expect(state.adt_quality_control.fetchMetrics()).toBeUndefined();
        expect(state.adt_normalization.fetchSizeFactors()).toBeNull();
        expect(state.adt_pca.fetchPCs()).toBeNull();
    }

    return;
}

export async function checkStateResultsBatched(state, { skipBasic = false } = {}) {
    if (!skipBasic) {
        await checkStateResultsBase(state);
    }

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
    {
        let res = state.quality_control.fetchFilters();
        let props = res.thresholdsSubsetProportions(0);
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
        expect(res["RNA"].numberOfBlocks()).toEqual(nlevels);
    }
}

export async function checkStateResultsAdt(state, { skipBasic = false } = {}) {
    if (!skipBasic) {
        await checkStateResultsBase(state);
    }

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
        let positive_total = 0;
        amet.subsetTotals(0, { copy: false }).forEach(x => { positive_total += (x > 0); });
        expect(positive_total).toBeGreaterThan(0);

        let afilt = state.adt_quality_control.fetchFilters();
        expect(afilt.thresholdsDetected()[0]).toBeGreaterThan(0);
        expect(afilt.thresholdsSubsetTotals(0)[0]).toBeGreaterThan(0);
    }

    // Cell filtering responds to both modalities.
    let nfiltered = state.cell_filtering.fetchFilteredMatrix().numberOfColumns();
    {
        let rna_only = 0;
        state.quality_control.fetchDiscards().forEach(x => { rna_only += (x > 0); });

        let adt_only = 0;
        state.adt_quality_control.fetchDiscards().forEach(x => { adt_only += (x > 0); });

        expect(nfiltered).toBeGreaterThan(rna_only);
        expect(nfiltered).toBeGreaterThan(adt_only);
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
    }

    // Combined embeddings.
    {
        let rna_dims = state.pca.fetchPCs().numberOfPCs();
        let adt_dims = state.adt_pca.fetchPCs().numberOfPCs();
        expect(state.combine_embeddings.fetchNumberOfDimensions()).toEqual(rna_dims + adt_dims);
        expect(state.combine_embeddings.fetchNumberOfCells()).toEqual(nfiltered);
    }

    // Markers.
    {
        let res = state.marker_detection.fetchResults();
        expect(res["ADT"].numberOfGroups()).toEqual(res["RNA"].numberOfGroups());
        let d = res["ADT"].cohen(0);
        expect(d.length).toBe(state.inputs.fetchCountMatrix().get("ADT").numberOfRows());
    }

    return;
}

export async function compareStates(left, right) {
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
    {
        let lmetrics = left.quality_control.fetchMetrics();
        let rmetrics = right.quality_control.fetchMetrics();
        expect(lmetrics.sums()).toEqual(rmetrics.sums());
        expect(lmetrics.detected()).toEqual(rmetrics.detected());
        expect(lmetrics.subsetProportions(0)).toEqual(rmetrics.subsetProportions(0));

        let lfilters = left.quality_control.fetchFilters();
        let rfilters = right.quality_control.fetchFilters();
        expect(lfilters.discardOverall()).toEqual(rfilters.discardOverall());
        expect(lfilters.thresholdsSums()).toEqual(rfilters.thresholdsSums());
    }

    // Cell filtering:
    {
        let ldiscard = left.cell_filtering.fetchDiscards();
        let rdiscard = right.cell_filtering.fetchDiscards();
        expect(ldiscard).toEqual(rdiscard);

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
    {
        let lmat = left.normalization.fetchNormalizedMatrix();
        let rmat = right.normalization.fetchNormalizedMatrix();

        let NR = lmat.numberOfRows();
        expect(NR).toEqual(rmat.numberOfRows());
        expect(lmat.row(0)).toEqual(rmat.row(0));
        expect(lmat.row(NR-1)).toEqual(rmat.row(NR-1));

        let lsf = left.normalization.fetchSizeFactors();
        let rsf = right.normalization.fetchSizeFactors();
        expect(lsf.array()).toEqual(rsf.array());
    }

    // Feature selection:
    {
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
    {
        let lpcs = left.pca.fetchPCs();
        let rpcs = right.pca.fetchPCs();
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

        let ng = lres.RNA.numberOfGroups();
        expect(ng).toEqual(rres.RNA.numberOfGroups());

        for (const a of available) {
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
