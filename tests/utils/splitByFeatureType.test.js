import * as bakana from "./../../src/index.js";

beforeAll(async () => await bakana.initialize({ localFile: true }));
afterAll(async () => await bakana.terminate());

const rna = "Gene Expression";
const adt = "Antibody Capture";

test("splitting by feature type works as expected with ADTs", () => {
    let reader = new bakana.TenxMatrixMarketDataset(
        "files/datasets/immune_3.0.0-matrix.mtx.gz",
        "files/datasets/immune_3.0.0-features.tsv.gz",
        "files/datasets/immune_3.0.0-barcodes.tsv.gz"
    );

    let deets = reader.load();
    expect(deets.matrix.get(rna).numberOfRows()).toBe(deets.features[rna].numberOfRows());
    expect(deets.matrix.get(adt).numberOfRows()).toBe(deets.features[adt].numberOfRows());

    // Everything remaining in the main matrix matches Ensembl.
    let guessed = bakana.callScran(scran => scran.guessFeatures(deets.features[rna].column("id")));
    expect(guessed.species).toBe("human");
    expect(guessed.type).toBe("ensembl");
    expect(guessed.confidence).toBe(1);

    // The ADT matrix just contains ADTs, which do look like gene symbols.
    let alt_guessed = bakana.callScran(scran => scran.guessFeatures(deets.features[adt].column("id")));
    expect(alt_guessed.species).toBe("human");
    expect(alt_guessed.type).toBe("symbol");
})

test("splitting by feature type works as expected when there are no ADTs", () => {
    let reader = new bakana.TenxMatrixMarketDataset(
        "files/datasets/pbmc3k-matrix.mtx.gz",
        "files/datasets/pbmc3k-features.tsv.gz",
        "files/datasets/pbmc3k-barcodes.tsv.gz"
    );

    let deets = reader.load();
    expect(deets.matrix.get(rna).numberOfRows()).toBe(deets.features[rna].numberOfRows());
    expect(deets.matrix.has(adt)).toBe(false);
})
