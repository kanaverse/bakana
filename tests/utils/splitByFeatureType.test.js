import * as bakana from "./../../src/index.js";
import * as mtx from "./../../src/readers/mtx.js";
import * as rutils from "./../../src/utils/reader.js";

beforeAll(async () => await bakana.initialize({ localFile: true }));
afterAll(async () => await bakana.terminate());

test("splitting by feature type works as expected with ADTs", () => {
    let reader = new mtx.Reader({
        format: "MatrixMarket",
        mtx: "files/datasets/immune_3.0.0-matrix.mtx.gz",
        genes: "files/datasets/immune_3.0.0-features.tsv.gz",
        annotations: "files/datasets/immune_3.0.0-barcodes.tsv.gz"
    });

    let deets = reader.load();
    expect(deets.matrix.get("RNA").numberOfRows()).toBe(deets.genes.RNA.id.length);
    expect(deets.matrix.get("ADT").numberOfRows()).toBe(deets.genes.ADT.id.length);

    // Everything remaining in the main matrix matches Ensembl.
    let guessed = bakana.callScran(scran => scran.guessFeatures(deets.genes.RNA.id));
    expect(guessed.species).toBe("human");
    expect(guessed.type).toBe("ensembl");
    expect(guessed.confidence).toBe(1);

    // The ADT matrix just contains ADTs, which do look like gene symbols.
    let alt_guessed = bakana.callScran(scran => scran.guessFeatures(deets.genes.ADT.id));
    expect(alt_guessed.species).toBe("human");
    expect(alt_guessed.type).toBe("symbol");
})

test("splitting by feature type works as expected when there are no ADTs", () => {
    let reader = new mtx.Reader({
        format: "MatrixMarket",
        mtx: "files/datasets/pbmc3k-matrix.mtx.gz",
        genes: "files/datasets/pbmc3k-features.tsv.gz",
        annotations: "files/datasets/pbmc3k-barcodes.tsv.gz"
    });

    let deets = reader.load();
    expect(deets.matrix.numberOfRows()).toBe(deets.genes.RNA.id.length);
    expect(deets.matrix.get("RNA").numberOfRows()).toBe(deets.genes.RNA.id.length);
    expect(deets.matrix.has("ADT")).toBe(false);
})
