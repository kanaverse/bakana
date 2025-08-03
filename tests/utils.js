import * as path from "path";
import * as fs from "fs";
import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as gesel from "gesel";
import * as wa from "wasmarrays.js";
import * as bioc from "bioconductor";

export async function initializeAll() {
    await bakana.initialize({ localFile: true });
}

export function purgeDirectory(path) {
    if (fs.existsSync(path)) {
        fs.rmSync(path, { recursive: true, force: true });
    }
}

export function baseParams() {
    let output = bakana.analysisDefaults();

    // Cut down on the work.
    output.rna_pca.num_pcs = 10;

    // Avoid getting held up by pointless iterations.
    output.tsne.iterations = 10;
    output.umap.num_epochs = 10;

    // Avoid loading in the feature sets.
    output.feature_set_enrichment.skip = true;

    // Avoid loading in the reference data.
    output.cell_labelling.references = [];

    return output;
}

/****************************
 **** Download overrides ****
 ****************************/

bakana.CellLabellingState.setDownload(url => {
    let fpath = path.basename(url);
    return fs.readFileSync("files/references/" + fpath).slice(); // Mimic a buffer from fetch().
});

bakana.RnaQualityControlState.setDownload(url => {
    let fpath = path.basename(url);
    return fs.readFileSync("files/mito-lists/" + fpath).slice(); // Mimic a buffer from fetch().
});

async function retrieve(file, old) {
    let cache_path = path.join("files", "feature-sets", file);
    let contents = fs.readFileSync(cache_path);
    let buffer = (new Uint8Array(contents)).buffer;
    return { ok: true, arrayBuffer: () => buffer }; // mimic Response object.
}

gesel.referenceDownload(retrieve);
gesel.geneDownload(retrieve);

/***********************************/

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

export function checkMatrixContents(matrix, names, loadedMatrix, loadedNames, { mustDiffer = true } = {}) {
    let NR = matrix.numberOfRows();
    let NC = matrix.numberOfColumns();
    if (loadedMatrix.numberOfRows() != NR || NC != loadedMatrix.numberOfColumns()) {
        throw new Error("loaded and reference matrix have different dimensions");
    }

    if (names.length != NR) {
        throw new Error("length of 'names' is different from number of matrix rows"); 
    }
    if (!is_same(names, loadedNames)) {
        throw new Error("loaded and reference names should be the same");
    }

    // Checking that the contents match up with the reference for every 100th column.
    for (var c = 0; c < NC; c += 100) {
        let loaded_first = loadedMatrix.column(c);
        let reference_first = matrix.column(c);
        if (!is_same(loaded_first, reference_first)) {
            throw new Error("loaded and reference matrix contents are not the same"); 
        }
    }
}

/***********************************/

function isDataFrameWithSimpleColumns(x) {
    expect(x instanceof bioc.DataFrame).toBe(true);
    for (const cn of x.columnNames()) {
        let col = x.column(cn);
        expect(col instanceof Array || (ArrayBuffer.isView(col) && !(col instanceof DataView))).toBe(true);
    }
}

function isArrayOfUniqueNames(x) {
    expect(x instanceof Array).toBe(true);
    for (const y of x) {
        expect(typeof y).toEqual("string") 
    }
    expect((new Set(x)).size).toEqual(x.length);
}

export async function checkDatasetGeneral(dataset) {
    expect(typeof dataset.constructor.format()).toEqual("string");
    expect(dataset.abbreviate().constructor).toBe(Object);
}

export async function checkDatasetSerialize(dataset) {
    let serialized = await dataset.serialize();
    expect(serialized.constructor).toBe(Object); 
    for (const fentry of serialized.files) {
        expect(typeof fentry.type).toEqual("string");
        expect(fentry.file instanceof bakana.SimpleFile).toEqual(true);
    }
    expect(serialized.options.constructor).toBe(Object); 

    let copy = await dataset.constructor.unserialize(serialized.files, serialized.options);
    expect(copy instanceof dataset.constructor).toBe(true);
    expect(copy.constructor.format()).toEqual(dataset.constructor.format());
    expect(copy.abbreviate()).toEqual(dataset.abbreviate());
    return copy;
}

export async function checkDatasetSummary(dataset) {
    let summ = await dataset.summary();
    if ("all_features" in summ) {
        expect("modality_features" in summ).toBe(false);
        isDataFrameWithSimpleColumns(summ.all_features);
        expect(summ.all_features.numberOfRows()).toBeGreaterThan(0);
    } else {
        expect("modality_features" in summ).toBe(true);
        for (const [mod, df] of Object.entries(summ.modality_features)) {
            isDataFrameWithSimpleColumns(df);
            expect(df.numberOfRows()).toBeGreaterThan(0);
        }
    }

    isDataFrameWithSimpleColumns(summ.cells);

    if ("all_assay_names" in summ) {
        expect("modality_assay_names" in summ).toBe(false);
        isArrayOfUniqueNames(summ.all_assay_names);
    } else {
        expect("all_assay_names" in summ).toBe(false);
        for (const [mod, nms] of Object.entries(summ.modality_assay_names)) {
            isArrayOfUniqueNames(nms);
        }
    }

    let uncached = await dataset.summary({ cache: false });
    expect(Object.keys(uncached)).toEqual(Object.keys(summ));
    return summ;
}

function sameDataFrame(left, right) {
    expect(left.numberOfRows()).toEqual(right.numberOfRows());
    expect(left.columnNames()).toEqual(right.columnNames());
}

export function sameDatasetSummary(left, right) {
    expect(Object.keys(left)).toEqual(Object.keys(right));

    if ("all_features" in left) {
        sameDataFrame(left.all_features, right.all_features);
    } else {
        expect(Object.keys(left.modality_features)).toEqual(Object.keys(right.modality_features));
        for (const [mod, df] of Object.entries(left.modality_features)) {
            sameDataFrame(df, right.modality_features[mod]);
        }
    }

    sameDataFrame(left.cells, right.cells);

    if ("all_assay_names" in left) {
        expect(left.all_assay_names).toEqual(right.all_assay_names);
    } else {
        expect(left.modality_assay_names).toEqual(right.modality_assay_names);
    }
}

export async function checkDatasetLoad(dataset) {
    let loaded = await dataset.load();

    expect(loaded.matrix.numberOfColumns()).toBeGreaterThan(0);
    let available_modalities = loaded.matrix.available();
    expect(available_modalities.length).toBeGreaterThan(0);
    for (const mod of available_modalities) {
        expect(["RNA", "ADT", "CRISPR"].indexOf(mod)).toBeGreaterThanOrEqual(0);
        let mat = loaded.matrix.get(mod);
        expect(mat.numberOfRows()).toBeGreaterThan(0);
        expect(mat.isSparse()).toEqual(true); // we force all count matrices to be sparse.
    }

    isDataFrameWithSimpleColumns(loaded.cells);
    expect(loaded.cells.numberOfRows()).toEqual(loaded.matrix.numberOfColumns());

    expect(available_modalities).toEqual(Object.keys(loaded.features));
    expect(available_modalities).toEqual(Object.keys(loaded.primary_ids));
    for (const mod of available_modalities) {
        const df = loaded.features[mod];
        expect(df.numberOfRows()).toEqual(loaded.matrix.get(mod).numberOfRows());
        isDataFrameWithSimpleColumns(df);
        const pid = loaded.primary_ids[mod];
        expect(df.numberOfRows()).toEqual(pid.length);
        expect(pid.every(y => typeof y === "string")).toBe(true);
    }

    let preview = await dataset.previewPrimaryIds();
    expect(preview).toEqual(loaded.primary_ids);

    let uncached = await dataset.load({ cache: false });
    expect(Object.keys(uncached)).toEqual(Object.keys(loaded));
    return loaded;
}

export function sameDatasetLoad(left, right) {
    let available_modalities = left.matrix.available();
    expect(available_modalities).toEqual(right.matrix.available());
    for (const mod of available_modalities) {
        expect(left.matrix.get(mod).numberOfRows()).toEqual(right.matrix.get(mod).numberOfRows());
        expect(left.matrix.get(mod).numberOfColumns()).toEqual(right.matrix.get(mod).numberOfColumns());
        expect(left.matrix.get(mod).row(0)).toEqual(right.matrix.get(mod).row(0));
        expect(left.matrix.get(mod).column(0)).toEqual(right.matrix.get(mod).column(0));
    }

    expect(Object.keys(left.features)).toEqual(available_modalities);
    expect(Object.keys(right.features)).toEqual(available_modalities);
    for (const mod of available_modalities) {
        sameDataFrame(left.features[mod], right.features[mod]);
    }

    sameDataFrame(left.cells, right.cells);
    expect(left.primary_ids).toEqual(right.primary_ids);
}

export async function checkResultSummary(result) {
    let summ = await result.summary();
    if ("all_features" in summ) {
        expect("modality_features" in summ).toBe(false);
        isDataFrameWithSimpleColumns(summ.all_features);
    } else {
        expect("modality_features" in summ).toBe(true);
        for (const [mod, df] of Object.entries(summ.modality_features)) {
            isDataFrameWithSimpleColumns(df);
        }
    }

    isDataFrameWithSimpleColumns(summ.cells);

    if ("all_assay_names" in summ) {
        expect("modality_assay_names" in summ).toBe(false);
        isArrayOfUniqueNames(summ.all_assay_names);
    } else {
        expect("all_assay_names" in summ).toBe(false);
        for (const [mod, nms] of Object.entries(summ.modality_assay_names)) {
            isArrayOfUniqueNames(nms);
        }
    }

    if ("reduced_dimension_names" in summ) {
        isArrayOfUniqueNames(summ.reduced_dimension_names);
    }
    if ("other_metadata" in summ) {
        expect(summ.other_metadata.constructor).toBe(Object);
    }

    let uncached = await result.summary({ cache: false });
    expect(Object.keys(uncached)).toEqual(Object.keys(summ));
    return summ;
}

export async function checkResultLoad(result) {
    let loaded = await result.load();

    expect(loaded.matrix.numberOfColumns()).toBeGreaterThan(0);
    let available_modalities = loaded.matrix.available();
    expect(available_modalities.length).toBeGreaterThan(0);
    for (const mod of available_modalities) {
        expect(loaded.matrix.get(mod).numberOfRows()).toBeGreaterThan(0);
    }

    isDataFrameWithSimpleColumns(loaded.cells);
    expect(loaded.cells.numberOfRows()).toEqual(loaded.matrix.numberOfColumns());

    expect(available_modalities).toEqual(Object.keys(loaded.features));
    for (const mod of available_modalities) {
        const df = loaded.features[mod];
        isDataFrameWithSimpleColumns(df);
        expect(df.numberOfRows()).toEqual(loaded.matrix.get(mod).numberOfRows());
    }

    if ("reduced_dimension_names" in loaded) {
        for (const [name, vals] of Object.entries(loaded)) {
            for (const dim of vals) {
                expect(dim instanceof Float64Array).toBe(true);
                expect(dim.length).toEqual(loaded.cells.numberOfRows());
            }
        }
    }
    if ("other_metadata" in loaded) {
        expect(loaded.other_metadata.constructor).toBe(Object);
    }

    let uncached = await result.load({ cache: false });
    expect(Object.keys(uncached)).toEqual(Object.keys(loaded));
    return loaded;
}

/***********************************/

export async function checkStateResultsMinimal(state, { mustFilter = true } = {}) {
    // Inputs:
    {
        let counts = state.inputs.fetchCountMatrix();
        expect(counts instanceof scran.MultiMatrix).toBe(true);
        let feat_anno = state.inputs.fetchFeatureAnnotations();

        for (const a of counts.available()){ 
            let NR = counts.get(a).numberOfRows();
            expect(a in feat_anno).toBe(true);
            expect(feat_anno[a].numberOfRows()).toBe(NR);
        }
    }

    let nfiltered = state.cell_filtering.fetchFilteredMatrix().numberOfColumns();

    // Cell filtering:
    {
        let counts = state.inputs.fetchCountMatrix();
        let keep = state.cell_filtering.fetchKeep();
        if (keep !== null) {
            expect(keep.length).toBe(counts.numberOfColumns());
            expect(keep instanceof wa.Uint8WasmArray).toBe(true);
        }

        let last_filtered = nfiltered - 1;
        let idx = [0, last_filtered];
        state.cell_filtering.undoFilter(idx);
        if (mustFilter) {
            expect(idx[1]).toBeGreaterThan(last_filtered);
        }

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

            let res0 = ares.cohensD(0);
            expect(res0 instanceof Float64Array).toBe(true);
            expect(res0.length).toEqual(ngenes);

            let resN = ares.auc(nclusters - 1);
            expect(resN instanceof Float64Array).toBe(true);
            expect(resN.length).toEqual(ngenes);
        }
    }

    return;
}

function checkSizeFactors(raw, norm, sf) {
    expect(sf instanceof wa.Float64WasmArray).toBe(true);
    expect(sf.length).toBe(norm.numberOfColumns());
    sf = sf.slice();

    // Checking every 200'th cell for consistent normalization.
    for (var c = 0; c < raw.numberOfColumns(); c += 200) {
        let expected = norm.column(c);
        let observed = raw.column(c);
        observed.forEach((x, i) => { observed[i] = Math.log2(observed[i] / sf[c] + 1); });

        // Cap precision for floating-point comparisons.
        observed.forEach((x, i) => { observed[i] = Math.round(x * 1000000) / 1000000; });
        expected.forEach((x, i) => { expected[i] = Math.round(x * 1000000) / 1000000; });
        expect(observed).toEqual(expected);
    }
}

export async function checkStateResultsRna(state, { exclusive = false, mustFilter = true, hasMito = true } = {}) {
    // Inputs:
    {
        let counts = state.inputs.fetchCountMatrix();
        expect(counts.has("RNA")).toBe(true);

        if (exclusive) {
            expect(state.inputs.fetchCountMatrix().has("ADT")).toBe(false);
            expect("ADT" in state.inputs.fetchFeatureAnnotations()).toBe(false);
        }
    }

    // Quality control:
    {
        let metres = state.rna_quality_control.fetchMetrics();
        expect(metres instanceof scran.PerCellRnaQcMetricsResults).toBe(true);

        let sumvec = metres.sum();
        expect(sumvec instanceof Float64Array).toBe(true);
        let ncells = state.inputs.fetchCountMatrix().numberOfColumns();
        expect(sumvec.length).toBe(ncells);

        // Check that the gene ID guessers find the mitochondrial genes.
        let mitovec = metres.subsetProportion(0);
        let sum = mitovec.reduce((a, b) => a+b);
        if (hasMito) {
            expect(sum).toBeGreaterThan(0);
        } else {
            expect(sum).toEqual(0);
        }

        expect(state.rna_quality_control.fetchKeep().length).toEqual(ncells);

        let filtres = state.rna_quality_control.fetchFilters();
        expect(filtres instanceof scran.SuggestRnaQcFiltersResults).toBe(true);
    }

    let nfiltered = state.cell_filtering.fetchFilteredMatrix().numberOfColumns();

    // Cell filtering:
    {
        let detvec = state.rna_quality_control.fetchMetrics().detected();
        if (mustFilter) {
            expect(detvec.length).toBeGreaterThan(nfiltered);
        }

        let refiltered = state.cell_filtering.applyFilter(detvec);
        expect(refiltered.length).toEqual(nfiltered);

        if (exclusive) {
            expect(state.cell_filtering.fetchKeep().owner).toBe(state.rna_quality_control.fetchKeep()); // i.e. a view
        }
    }

    let ngenes = state.inputs.fetchCountMatrix().get("RNA").numberOfRows();

    // Normalization:
    {
        let normed = state.rna_normalization.fetchNormalizedMatrix();
        expect(normed instanceof scran.ScranMatrix).toBe(true);
        expect(normed.numberOfColumns()).toBe(nfiltered);
        expect(normed.numberOfRows()).toBe(ngenes);

        let sf = state.rna_normalization.fetchSizeFactors();
        checkSizeFactors(state.cell_filtering.fetchFilteredMatrix().get("RNA"), normed, sf);
    }

    // Feature selection:
    {
        let feats = state.feature_selection.fetchResults();
        expect(feats instanceof scran.ModelGeneVariancesResults).toBe(true);

        let resids = feats.residuals();
        expect(resids.length).toBe(ngenes);

        let sresids = state.feature_selection.fetchSortedResiduals();
        resids.sort((a, b) => a - b);
        expect(sresids).toEqual(resids);
    }

    // PCA:
    {
        let pcs = state.rna_pca.fetchPCs();
        expect(pcs instanceof scran.RunPcaResults).toBe(true);
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
        expect(filtres.sum().length).toBe(1);
    }
    if (state.adt_quality_control.valid()) {
        let filtres = state.adt_quality_control.fetchFilters();
        expect(filtres.detected().length).toBe(1);
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
        let props = res.subsetProportion(0);
        expect(props.length).toEqual(nlevels);
    }
    if (state.adt_quality_control.valid()) {
        let res = state.adt_quality_control.fetchFilters();
        let props = res.subsetSum(0);
        expect(props.length).toEqual(nlevels);
    }
    if (state.crispr_quality_control.valid()) {
        let res = state.crispr_quality_control.fetchFilters();
        let props = res.maxValue();
        expect(props.length).toEqual(nlevels);
    }

    // Check that some correction actually occured.
    {
        let corrected_pcs = state.batch_correction.fetchCorrected();
        let original_pcs = state.combine_embeddings.fetchCombined();
        expect(corrected_pcs.length).toEqual(original_pcs.length);
        expect(corrected_pcs.array()).not.toEqual(original_pcs.array());
    }

    return;
}

export function checkStateResultsAdt(state, { exclusive = false } = {}) {
    // Checking that ADTs exist.
    {
        expect(state.inputs.fetchCountMatrix().has("ADT")).toBe(true);
        let adtmat = state.inputs.fetchCountMatrix().get("ADT");
        expect(state.inputs.fetchFeatureAnnotations()["ADT"].numberOfRows()).toEqual(adtmat.numberOfRows());
    }

    // Checking the QCs.
    {
        let amet = state.adt_quality_control.fetchMetrics();
        expect(amet instanceof scran.PerCellAdtQcMetricsResults).toBe(true);

        let positive_total = 0;
        amet.subsetSum(0, { copy: false }).forEach(x => { positive_total += (x > 0); });
        expect(positive_total).toBeGreaterThan(0);

        // Check that the feature ID guessers find the IgGs.
        let igvec = amet.subsetSum(0);
        let sum = igvec.reduce((a, b) => a+b);
        expect(sum).toBeGreaterThan(0);

        expect(state.adt_quality_control.fetchKeep().length).toEqual(state.inputs.fetchCountMatrix().numberOfColumns());

        let afilt = state.adt_quality_control.fetchFilters();
        expect(afilt instanceof scran.SuggestAdtQcFiltersResults).toBe(true);
        expect(afilt.detected()[0]).toBeGreaterThan(0);
        expect(afilt.subsetSum(0)[0]).toBeGreaterThan(0);
    }

    let nfiltered = state.cell_filtering.fetchFilteredMatrix().numberOfColumns();

    // Cell filtering:
    {
        let detvec = state.adt_quality_control.fetchMetrics().detected();
        expect(detvec.length).toBeGreaterThan(nfiltered);
        let refiltered = state.cell_filtering.applyFilter(detvec);
        expect(refiltered.length).toEqual(nfiltered);

        if (exclusive) {
            expect(state.cell_filtering.fetchKeep().owner).toBe(state.adt_quality_control.fetchKeep()); // i.e. a view
        }
    }

    // Normalization.
    {
        let sf = state.adt_normalization.fetchSizeFactors();
        let positive_total = 0;
        sf.forEach(x => { positive_total += (x > 0); });
        expect(positive_total).toBeGreaterThan(0);

        let norm = state.adt_normalization.fetchNormalizedMatrix();
        expect(norm.numberOfColumns()).toEqual(sf.length);
        expect(norm.numberOfRows()).toEqual(state.inputs.fetchCountMatrix().get("ADT").numberOfRows());

        checkSizeFactors(state.cell_filtering.fetchFilteredMatrix().get("ADT"), norm, sf);
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
    }

    // Checking the QCs.
    {
        let cmet = state.crispr_quality_control.fetchMetrics();
        expect(cmet instanceof scran.PerCellCrisprQcMetricsResults).toBe(true);

        expect(state.crispr_quality_control.fetchKeep().length).toEqual(state.inputs.fetchCountMatrix().numberOfColumns());

        let cfilt = state.crispr_quality_control.fetchFilters();
        expect(cfilt instanceof scran.SuggestCrisprQcFiltersResults).toBe(true);
        expect(cfilt.maxValue()[0]).toBeGreaterThan(0);
    }

    let nfiltered = state.cell_filtering.fetchFilteredMatrix().numberOfColumns();

    // Cell filtering:
    {
        let detvec = state.crispr_quality_control.fetchMetrics().detected();
        expect(detvec.length).toBeGreaterThan(nfiltered);
        let refiltered = state.cell_filtering.applyFilter(detvec);
        expect(refiltered.length).toEqual(nfiltered);

        if (exclusive) {
            expect(state.cell_filtering.fetchKeep().owner).toBe(state.crispr_quality_control.fetchKeep()); // i.e. a view
        }
    }

    // Normalization.
    {
        let sf = state.crispr_normalization.fetchSizeFactors();
        let positive_total = 0;
        sf.forEach(x => { positive_total += (x > 0); });
        expect(positive_total).toBeGreaterThan(0);

        let norm = state.crispr_normalization.fetchNormalizedMatrix();
        expect(norm.numberOfColumns()).toEqual(sf.length);
        expect(norm.numberOfRows()).toEqual(state.inputs.fetchCountMatrix().get("CRISPR").numberOfRows());

        checkSizeFactors(state.cell_filtering.fetchFilteredMatrix().get("CRISPR"), norm, sf);
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
        expect(state.cell_filtering.fetchKeep().owner).toBeNull(); // i.e. not a view.

        let rna_only = 0;
        state.rna_quality_control.fetchKeep().forEach(x => { rna_only += (x > 0); });

        let adt_only = 0;
        state.adt_quality_control.fetchKeep().forEach(x => { adt_only += (x > 0); });

        expect(nfiltered).toBeLessThan(rna_only);
        expect(nfiltered).toBeLessThan(adt_only);
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
        expect(state.cell_filtering.fetchKeep().owner).toBeNull(); // i.e. not a view.

        let rna_only = 0;
        state.rna_quality_control.fetchKeep().forEach(x => { rna_only += (x > 0); });

        let crispr_only = 0;
        state.crispr_quality_control.fetchKeep().forEach(x => { crispr_only += (x > 0); });

        expect(nfiltered).toBeLessThan(rna_only);
        expect(nfiltered).toBeLessThan(crispr_only);
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
        expect(lmetrics.sum()).toEqual(rmetrics.sum());
        expect(lmetrics.detected()).toEqual(rmetrics.detected());
        expect(lmetrics.subsetProportion(0)).toEqual(rmetrics.subsetProportion(0));

        let lfilters = left.rna_quality_control.fetchFilters();
        let rfilters = right.rna_quality_control.fetchFilters();
        expect(lfilters.sum()).toEqual(rfilters.sum());
        expect(lfilters.subsetProportion(0)).toEqual(rfilters.subsetProportion(0));

        let ldiscards = left.rna_quality_control.fetchKeep().array();
        let rdiscards = right.rna_quality_control.fetchKeep().array();
        expect(ldiscards).toEqual(rdiscards);
    }

    if (checkAdt) {
        let lmetrics = left.adt_quality_control.fetchMetrics();
        let rmetrics = right.adt_quality_control.fetchMetrics();
        expect(lmetrics.sum()).toEqual(rmetrics.sum());
        expect(lmetrics.detected()).toEqual(rmetrics.detected());
        expect(lmetrics.subsetSum(0)).toEqual(rmetrics.subsetSum(0));

        let lfilters = left.adt_quality_control.fetchFilters();
        let rfilters = right.adt_quality_control.fetchFilters();
        expect(lfilters.detected()).toEqual(rfilters.detected());
        expect(lfilters.subsetSum(0)).toEqual(rfilters.subsetSum(0));

        let ldiscards = left.adt_quality_control.fetchKeep().array();
        let rdiscards = right.adt_quality_control.fetchKeep().array();
        expect(ldiscards).toEqual(rdiscards);
    }

    if (checkCrispr) {
        let lmetrics = left.crispr_quality_control.fetchMetrics();
        let rmetrics = right.crispr_quality_control.fetchMetrics();
        expect(lmetrics.sum()).toEqual(rmetrics.sum());
        expect(lmetrics.detected()).toEqual(rmetrics.detected());
        expect(lmetrics.maxProportion()).toEqual(rmetrics.maxProportion());

        let lfilters = left.crispr_quality_control.fetchFilters();
        let rfilters = right.crispr_quality_control.fetchFilters();
        expect(lfilters.maxValue()).toEqual(rfilters.maxValue());

        let ldiscards = left.crispr_quality_control.fetchKeep().array();
        let rdiscards = right.crispr_quality_control.fetchKeep().array();
        expect(ldiscards).toEqual(rdiscards);
    }

    // Cell filtering:
    {
        let ldiscard = left.cell_filtering.fetchKeep();
        let rdiscard = right.cell_filtering.fetchKeep();
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
                expect(lres[a].cohensD(g)).toEqual(rres[a].cohensD(g));
                expect(lres[a].auc(g)).toEqual(rres[a].auc(g));
                expect(lres[a].mean(g)).toEqual(rres[a].mean(g));
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
    let lfcs = vres.results.RNA.deltaMean(vres.left);

    let vres2 = state.marker_detection.computeVersus(0, nclusters - 1);
    expect(vres2.results.RNA instanceof scran.ScoreMarkersResults).toBe(true);
    expect(vres2.left).toBe(vres.right);

    let lfcs2 = vres2.results.RNA.deltaMean(vres2.left);
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

export function mockBlocks(input, output, nblocks) {
    // Mocking up a blocking file with pretend batches.
    let f = fs.readFileSync(input);
    let buff = f.buffer.slice(f.byteOffset, f.byteOffset + f.byteLength);
    let stuff = bakana.readTable(new Uint8Array(buff));

    let ncells = stuff.length;
    let per_block = Math.ceil(ncells / nblocks);
    let blocks = new Array(ncells);
    for (var c = 0; c < ncells; c++) {
        blocks[c] = 'A' + String(Math.floor(c / per_block));
    }

    fs.writeFileSync(output, blocks.join("\n"));
    return "A0";
}

