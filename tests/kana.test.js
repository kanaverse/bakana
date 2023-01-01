import * as bakana from "../src/index.js";
import * as scran from "scran.js";
import * as utils from "./utils.js";
import * as fs from "fs";

beforeAll(utils.initializeAll);
afterAll(async () => await bakana.terminate());

let files = {
    default: new bakana.TenxMatrixMarketDataset(
            "files/datasets/pbmc3k-matrix.mtx.gz",
            "files/datasets/pbmc3k-features.tsv.gz",
            "files/datasets/pbmc3k-barcodes.tsv.gz"
        )
};

test("saving to and loading from a kana file works correctly (embedded)", async () => {
    let params = utils.baseParams();
    let state = await bakana.createAnalysis();
    await bakana.runAnalysis(state, files, params);
    await utils.overlordCheckStandard(state);

    const path = "TEST_kana_state.h5";
    let collected = await bakana.saveAnalysis(state, path);
    expect(collected.total > 0).toBe(true);

    const kpath = "TEST_state.kana";
    const round_trip = "TEST_kana_state_again.h5";

    // Saving to a kana file.
    {
        let res = await bakana.createKanaFile(path, collected.collected, { outputPath: kpath });
        expect(fs.statSync(res).size).toBe(24 + fs.statSync(path).size + collected.total);

        // Alright - trying to unpack everything.
        let tmp = "TEST_kana_loader";
        if (fs.existsSync(tmp)) {
            fs.rmSync(tmp, { recursive: true, force: true });
        }
        fs.mkdirSync(tmp);

        let loader = await bakana.parseKanaFile(kpath, round_trip, { stageDir: tmp });
        expect(loader.version).toBe(bakana.kanaFormatVersion);
        expect(loader.embedded).toBe(true);

        utils.validateState(round_trip);
        let reloaded = await bakana.loadAnalysis(round_trip, loader.load);
        await utils.compareStates(state, reloaded);

        // Cleaning up.
        await bakana.freeAnalysis(reloaded);
        scran.removeFile(kpath);
        scran.removeFile(round_trip);
    }

    // Round 2: converting all collected files into buffers.
    {
        let collected2 = collected.collected.map(x => (new bakana.SimpleFile(x)).buffer());
        let res = await bakana.createKanaFile(path, collected2, { outputPath: kpath });
        expect(fs.statSync(res).size).toBe(24 + fs.statSync(path).size + collected.total);

        let contents = fs.readFileSync(kpath); // Buffer is subclass of Uint8Array already!
        let loader = await bakana.parseKanaFile(contents, round_trip);
        utils.validateState(round_trip);
        let reloaded = await bakana.loadAnalysis(round_trip, loader.load);
        await utils.compareStates(state, reloaded);

        // Cleaning up.
        await bakana.freeAnalysis(reloaded);
        scran.removeFile(kpath);
        scran.removeFile(round_trip);
    }

    scran.removeFile(path);
    await bakana.freeAnalysis(state);
})

test("saving to and loading from a kana file works with links", async () => {
    let params = utils.baseParams();
    let state = await bakana.createAnalysis();
    await bakana.runAnalysis(state, files, params);
    await utils.overlordCheckStandard(state);

    // Links just re-use the file path for our Node tests, which is unique enough!
    let old_create = bakana.setCreateLink((format, file) => format + ":" + file.content());
    let old_resolve = bakana.setResolveLink(id => id.split(":")[1])

    // Saving to a kana file.
    const path = "TEST_kana_state2.h5";
    await bakana.saveAnalysis(state, path, { embedded: false });

    const kpath = "TEST_state2.kana";
    let res = await bakana.createKanaFile(path, null, { outputPath: kpath });
    expect(fs.statSync(res).size).toBe(24 + fs.statSync(path).size);

    scran.removeFile(path);
    expect(fs.existsSync(path)).toBe(false); // properly removed.

    // Alright - trying to unpack everything.
    const round_trip = "TEST_kana_state_again2.h5";
    let loader = await bakana.parseKanaFile(kpath, round_trip);
    expect(loader.version).toBe(bakana.kanaFormatVersion);
    expect(loader.embedded).toBe(false);

    utils.validateState(round_trip, false);
    let reloaded = await bakana.loadAnalysis(round_trip, loader.load);
    await utils.compareStates(state, reloaded);

    // Reverting the linkers.
    bakana.setCreateLink(old_create);
    bakana.setResolveLink(old_resolve);

    // Releasing all memory.
    await bakana.freeAnalysis(state);
    await bakana.freeAnalysis(reloaded);

    // Deleting the files.
    scran.removeFile(path);
    scran.removeFile(round_trip);
})
