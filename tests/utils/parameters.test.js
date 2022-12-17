import * as params from "../../src/defaults.js";

test("configureBatchCorrection works as expected", () => {
    let stuff = params.analysisDefaults();

    params.configureBatchCorrection(stuff, "none");
    expect(stuff.batch_correction.method).toEqual('none');
    expect(stuff.rna_pca.block_method).toEqual('none');
    expect(stuff.adt_pca.block_method).toEqual('none');
    expect(params.guessBatchCorrectionConfig(stuff)).toBe("none");

    params.configureBatchCorrection(stuff, "regress");
    expect(stuff.batch_correction.method).toEqual('none');
    expect(stuff.rna_pca.block_method).toEqual('regress');
    expect(stuff.adt_pca.block_method).toEqual('regress');
    expect(params.guessBatchCorrectionConfig(stuff)).toBe("regress");

    params.configureBatchCorrection(stuff, "mnn");
    expect(stuff.batch_correction.method).toEqual('mnn');
    expect(stuff.rna_pca.block_method).toEqual('weight');
    expect(stuff.adt_pca.block_method).toEqual('weight');
    expect(params.guessBatchCorrectionConfig(stuff)).toBe("mnn");

    expect(() => params.configureBatchCorrection(stuff, "foo")).toThrow("unknown");
});

test("guessBatchCorrection handles ambiguity properly", () => {
    let stuff = params.analysisDefaults();

    stuff.batch_correction.method = "mnn";
    expect(params.guessBatchCorrectionConfig(stuff)).toBe("mnn");
    expect(params.guessBatchCorrectionConfig(stuff, { strict: true })).toBeNull();

    stuff.batch_correction.method = "none";
    stuff.rna_pca.block_method = "regress";
    expect(params.guessBatchCorrectionConfig(stuff)).toBe("regress");
    expect(params.guessBatchCorrectionConfig(stuff, { strict: true })).toBeNull();

    stuff.rna_pca.block_method = "none";
    stuff.adt_pca.block_method = "weight";
    expect(params.guessBatchCorrectionConfig(stuff)).toBe("none");
    expect(params.guessBatchCorrectionConfig(stuff, { strict: true })).toBeNull();

    stuff.rna_pca.block_method = "weight";
    expect(params.guessBatchCorrectionConfig(stuff)).toBe("none");
    expect(params.guessBatchCorrectionConfig(stuff, { strict: true })).toBeNull();
});

test("configureApproximateNeighbors works as expected", () => {
    let stuff = params.analysisDefaults();

    params.configureApproximateNeighbors(stuff, false);
    expect(stuff.batch_correction.approximate).toBe(false);
    expect(stuff.combine_embeddings.approximate).toBe(false);
    expect(stuff.neighbor_index.approximate).toBe(false);
    expect(params.guessApproximateNeighborsConfig(stuff)).toBe(false);

    params.configureApproximateNeighbors(stuff, true);
    expect(stuff.batch_correction.approximate).toBe(true);
    expect(stuff.combine_embeddings.approximate).toBe(true);
    expect(stuff.neighbor_index.approximate).toBe(true);
    expect(params.guessApproximateNeighborsConfig(stuff)).toBe(true);

    // Handles ambiguity properly.
    stuff.batch_correction.approximate = false;
    expect(params.guessApproximateNeighborsConfig(stuff)).toBe(true);
    expect(params.guessApproximateNeighborsConfig(stuff, { strict: true })).toBeNull();
});
