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
    expect(deets.matrix.numberOfRows()).toBe(deets.genes.type.length);
    expect("Antibody Capture" in deets.alternatives).toBe(true);

    // Everything remaining in the main matrix matches Ensembl.
    let guessed = bakana.callScran(scran => scran.guessFeatures(deets.genes.id));
    expect(guessed.species).toBe("human");
    expect(guessed.type).toBe("ensembl");
    expect(guessed.confidence).toBe(1);

    let found = Array.from(new Set(deets.genes.type));
    expect(found.length).toBe(1);
    expect(found[0]).toBe("Gene Expression");

    // The ADT matrix just contains ADTs, which do look like gene symbols.
    let alt = deets.alternatives["Antibody Capture"];
    expect(alt.matrix.numberOfRows()).toBe(alt.genes.id.length);

    let alt_guessed = bakana.callScran(scran => scran.guessFeatures(alt.genes.id));
    expect(alt_guessed.species).toBe("human");
    expect(alt_guessed.type).toBe("symbol");

    let alt_found = Array.from(new Set(alt.genes.type));
    expect(alt_found.length).toBe(1);
    expect(alt_found[0]).toBe("Antibody Capture");
})

test("splitting by feature type works as expected when there are no ADTs", () => {
    let reader = new mtx.Reader({
        format: "MatrixMarket",
        mtx: "files/datasets/pbmc3k-matrix.mtx.gz",
        genes: "files/datasets/pbmc3k-features.tsv.gz",
        annotations: "files/datasets/pbmc3k-barcodes.tsv.gz"
    });

    let deets = reader.load();
    expect(deets.matrix.numberOfRows()).toBe(deets.genes.type.length);
    expect("alternative" in deets).toBe(false);
})
