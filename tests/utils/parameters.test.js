import * as params from "../../src/defaults.js";

test("configureBatchCorrection works as expected", () => {
    let stuff = params.analysisDefaults();

    params.configureBatchCorrection(stuff, "none");
    expect(stuff.batch_correction.method).toEqual('none');
    expect(stuff.pca.block_method).toEqual('none');
    expect(stuff.adt_pca.block_method).toEqual('none');

    params.configureBatchCorrection(stuff, "regress");
    expect(stuff.batch_correction.method).toEqual('none');
    expect(stuff.pca.block_method).toEqual('regress');
    expect(stuff.adt_pca.block_method).toEqual('regress');

    params.configureBatchCorrection(stuff, "mnn");
    expect(stuff.batch_correction.method).toEqual('mnn');
    expect(stuff.pca.block_method).toEqual('weight');
    expect(stuff.adt_pca.block_method).toEqual('weight');

    expect(() => params.configureBatchCorrection(stuff, "foo")).toThrow("unknown");
});

test("configureApproximateNeighbors works as expected", () => {
    let stuff = params.analysisDefaults();

    params.configureApproximateNeighbors(stuff, false);
    expect(stuff.batch_correction.approximate).toBe(false);
    expect(stuff.combine_embeddings.approximate).toBe(false);
    expect(stuff.neighbor_index.approximate).toBe(false);

    params.configureApproximateNeighbors(stuff, true);
    expect(stuff.batch_correction.approximate).toBe(true);
    expect(stuff.combine_embeddings.approximate).toBe(true);
    expect(stuff.neighbor_index.approximate).toBe(true);
});
