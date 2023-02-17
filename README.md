# Backend for kana

## Overview

**bakana** provides the compute backend for the [**kana**](https://github.com/jkanche/kana) application.
It provides a pipeline for a routine single-cell RNA-seq analysis, starting from the count matrix and finishing with the usual results (markers, clusters, t-SNE and so on).
Datasets involving multiple samples in one or multiple matrices can also be analyzed with blocking and batch correction.
The pipeline can be executed both in the browser and on Node.js.
It supports in-memory caching of the analysis state for fast iterative re-analysis,
as well as serialization of the state for storage and distribution to other machines.

## Getting started

Install the package from [NPM](https://npmjs.com/package/bakana) using the usual method:

```sh
npm install bakana
```

See the [reference documentation](https://ltla.github.io/bakana/) for details on available functions.

## Running analyses

We perform an analysis with the following commands:

```js
import * as bakana from "bakana";
await bakana.initialize({ numberOfThreads: 8 }); 

let state = await bakana.createAnalysis();
let params = bakana.analysisDefaults();

await bakana.runAnalysis(state, 
    // Specify files using paths (Node.js) or File objects (browser).
    { my_data: new bakana.TenxHdf5Dataset("/some/file/path.h5") },
    params
);
```

Each step is represented by a `*State` class instance as a property of `state`, containing the analysis results.
We can extract results from each state for further inspection.

```js
state.rna_normalization.fetchNormalizedMatrix(); 
// ScranMatrix {}

state.rna_quality_control.fetchMetrics(); 
// PerCellRnaQcMetricsResults {}

state.rna_pca.fetchPCs(); 
// RunPCAResults {}
```

We can also supply a callback that processes summaries from each step as soon as it finishes.
This can be useful for, e.g., posting diagnostics to another system once they become available.

```js
function finishCallback(step) {
    if (state[step].changed) {
        console.log("Running " + step);
        // Do other stuff with the state's results...
    }
}

await bakana.runAnalysis(state, 
    { my_data: new bakana.TenxHdf5Dataset("/some/file/path.h5") }
    params,
    { finishFun: finishCallback }
);
```

If the analysis is re-run with different parameters, **bakana** will only re-run the affected steps.
This includes all steps downstream of any step with changed parameters.

```js
params.rna_pca.num_pcs = 15;
await bakana.runAnalysis(state, 
    { my_data: new bakana.TenxHdf5Dataset("/some/file/path.h5") }
    params,
    { finishFun: finishCallback }
);
// Running rna_pca
// Running combine_embeddings
// Running batch_correction
// Running neighbor_index
// Running kmeans_cluster
// Running snn_graph_cluster
// ...
```

## Saving results

Given an analysis state, we can dump its contents into a `SingleCellExperiment` for further examination:

```js
await bakana.saveSingleCellExperiment(state, "sce", { directory: "output" });
```

This stores the data and results into various fields of the `SingleCellExperiment`:

- The assays contain the sparse (QC-filtered) count matrix and its corresponding log-transformed normalized matrix as a `DelayedArray`.
- The row data contains gene identifiers along with variance modelling and marker detection results in nested `DataFrame`s.
- The column data contains quality control metrics in nested `DataFrame`s, along with other pieces like the clustering and blocking.
- The reduced dimensions contains PCA, t-SNE and UMAP results.
- The alternative experiments contains further nested `SingleCellExperiment`s for other modalities.

We use the [**alabaster**](https://github.com/ArtifactDB/alabaster.base) representation of a `SingleCellExperiment` to provide multi-language access to the results.
For example, we can load the `SingleCellExperiment` back into an R session via the `loadObject()` function:

```r
library(alabaster.base)
info <- acquireMetadata("output", "sce")
sce <- loadObject(info, "output")
```

## Saving configurations

Given an analysis state, we can save its configuration via the `serializeAnalysis()` function.
This returns an object that contains the analysis parameters, which can then be converted to JSON and saved to file.

```js
let saved = [];
let saveFileHandler = (k, f, file) => {
    saved.push(file.buffer());
    return String(saved.length);
};

let config = await bakana.serializeConfiguration(state, saveFileHandler);
```

Applications are responsible for deciding how to handle the input data files.
In the example above, we just store the file contents in a `saved` array of Uint8Arrays, e.g., for inclusion in a tarball with the configuration JSON.
More complex applications may create a staging directory on the file system in which to store the files (e.g., for Node.js),
or may register the file contents in a database for later extraction.

It is similarly easy to use a configuration to create a new analysis state via the `unserializeConfiguration()` function.
This will extract the parameters/data files and rerun the entire analysis via `runAnalysis()`,
allowing us to recover the same analysis state that went into `serializeConfiguration()`.
The example below uses a loading handler that just undoes the effect of `saveFileHandler`.

```js
let loadFileHandler = id => saved[Number(id) - 1];
let reloaded = await bakana.unserializeConfiguration(config, loadFileHandler);
```

## Terminating analyses

Once a particular analysis is finished, we should free the resources of its state.

```js
bakana.freeAnalysis(state);
bakana.freeAnalysis(reloaded);
bakana.freeAnalysis(reloaded2);
```

If all analyses are complete, we can terminate the entire session.
This is necessary on Node.js to allow the runtime to exit properly.

```js
bakana.terminate();
```

## Developer notes

See [here](docs/related/custom_readers.md) for instructions on adding custom dataset readers.
This allows us to use **bakana**'s analysis pipeline and serialization capabilities on datasets from other sources such as in-house databases.

Testing can be done with `npm run test` with Node 16+.
For older versions of Node, it requires some combination of the options below:

```sh
node --experimental-vm-modules \
    --experimental-wasm-threads \
    --experimental-wasm-bulk-memory \
    --experimental-wasm-bigint \
    node_modules/jest/bin/jest.js \
    --testTimeout=100000000 \
    --runInBand
```
