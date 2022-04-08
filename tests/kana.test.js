import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";
import * as fs from "fs";

beforeAll(async () => await bakana.initialize({ localFile: true }));
afterAll(async () => await bakana.terminate());

let files = { 
    default: {
        format: "MatrixMarket",
        mtx: "files/datasets/pbmc3k-matrix.mtx.gz",
        genes: "files/datasets/pbmc3k-features.tsv.gz",
        annotations: "files/datasets/pbmc3k-barcodes.tsv.gz"
    }
}

test("saving to and loading from a kana file works correctly (embedded)", async () => {
    let params = utils.baseParams();
    let state = await bakana.createAnalysis();
    await bakana.runAnalysis(state, files, params);
    let ref_pca = state.pca.summary();

    // Saving to a kana file.
    const path = "TEST_kana_state.h5";
    let collected = await bakana.saveAnalysis(state, path);
    bakana.freeAnalysis(state);

    const kpath = "TEST_state.kana";
    let res = await bakana.createKanaFile(path, collected.collected, { outputPath: kpath });
    expect(collected.total > 0).toBe(true);
    expect(fs.statSync(res).size).toBe(24 + fs.statSync(path).size + collected.total);

    // Alright - trying to unpack everything.
    let tmp = "TEST_kana_loader";
    if (fs.existsSync(tmp)) {
        fs.rmSync(tmp, { recursive: true, force: true });
    }
    fs.mkdirSync(tmp);
    const round_trip = "TEST_kana_state_again.h5";

    let loader = await bakana.parseKanaFile(kpath, round_trip, { stageDir: tmp });
    let reloaded = await bakana.loadAnalysis(round_trip, loader);

    // Checking that we got something that was reasonable.
    expect(reloaded.state.pca.summary()).toEqual(ref_pca);
})

test("saving to and loading from a kana file works with links", async () => {
    let params = utils.baseParams();
    let state = await bakana.createAnalysis();
    await bakana.runAnalysis(state, files, params);
    let ref_pca = state.pca.summary();

    // Links just re-use the file path for our Node tests, which is unique enough!
    let old_create = bakana.setCreateLink(x => x);
    let old_resolve = bakana.setResolveLink(x => x)

    // Saving to a kana file.
    const path = "TEST_kana_state2.h5";
    await bakana.saveAnalysis(state, path, { embedded: false });
    bakana.freeAnalysis(state);

    const kpath = "TEST_state2.kana";
    let res = await bakana.createKanaFile(path, null, { outputPath: kpath });
    expect(fs.statSync(res).size).toBe(24 + fs.statSync(path).size);

    // Alright - trying to unpack everything.
    const round_trip = "TEST_kana_state_again2.h5";
    let loader = await bakana.parseKanaFile(kpath, round_trip);
    let reloaded = await bakana.loadAnalysis(round_trip, loader);

    // Checking that we got something that was reasonable.
    expect(reloaded.state.pca.summary()).toEqual(ref_pca);

    // Reverting the linkers.
    bakana.setCreateLink(old_create);
    bakana.setResolveLink(old_resolve);
})
