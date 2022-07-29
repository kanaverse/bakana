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
    { my_data: { format: "10X", h5: "/some/file/path" } },
    params
);
```

Each step is represented by a `*State` class instance as a property of `state`, containing the analysis results.
We can extract summaries for diagnostics or display on a UI:

```js
state.quality_control.summary();
// {
//   data: {
//     default: {
//       sums: [Float64Array],
//       detected: [Int32Array],
//       proportion: [Float64Array]
//     }
//   },
//   ...
// }

state.pca.summary():
// {
//   var_exp: Float64Array(20) [...]
// }

// Some summary() are promises, see the relevant documentation.
await state.tsne.summary();
// {
//   x: [Float64Array],
//   y: [Float64Array],
//   iterations: 500
// }
```

Alternatively, we can supply a callback that processes summaries from each step as soon as it finishes.
This can be useful for, e.g., posting diagnostics to another system once they become available.

```js
function finisher(step, results) {
    console.log(step);
    // Do other stuff with 'results' here.
}

await bakana.runAnalysis(state, 
    // Specify files using paths (Node.js) or File objects (browser).
    { my_data: { format: "10X", h5: "/some/file/path" } },
    params,
    { finishFun: finisher }
);
// inputs
// quality_control
// normalizaton
// feature_selection
// pca
// neighbor_index
// kmeans_cluster
// ...
```

If the analysis is re-run with different parameters, **bakana** will only re-run the affected steps.
This includes all steps downstream of any step with changed parameters.

```js
params.pca.num_pcs = 15;

await bakana.runAnalysis(state, 
    // Specify files using paths (Node.js) or File objects (browser).
    { my_data: { format: "10X", h5: "/some/file/path" } },
    params,
    { finishFun: finisher }
);
// pca
// neighbor_index
// ...
```

## Saving analyses

Given an analysis state, we can save the parameters and results to a HDF5 file.

```js
let collected = await bakana.saveAnalysis(state, "whee.h5");
```

`saveAnalysis` will return a promise that resolves to an array of paths (Node.js) or `File` objects for the original data files.
These can be used to assemble a `*.kana`-format file following the [**kanaval** specification](https://ltla.github.io/kanaval).
By default, this assumes that we are embedding the data files into the `*.kana` file.

```js
let res = await bakana.createKanaFile("whee.h5", collected.collected);
```

Advanced users may prefer to store links to the data files rather than embedding them.
This can be achieved using the `setCreateLink()` function to define a mechanism for creating application-specific links.
For example, the **kana** application uses IndexedDB to cache each file for later use in the browser.

```js
bakana.setCreateLink(save_to_some_db);
bakana.setResolveLink(load_from_some_db);
```

## Loading analyses

Given a `*.kana` file, it is straightforward to extract the various state and data files.

```js
let loader = await bakana.parseKanaFile("something.kana", "foo.h5");
```

Once the `*.kana` file is parsed, we can reload the analysis state into memory.
We can also report the parameters used at each step.

```js
let reloaded = await bakana.loadAnalysis("foo.h5", loader);
let params = bakana.retrieveParameters(reloaded);
```

These can be used to perform a new analysis, possibly after modifying some parameters.

```js
params.pca.num_hvgs = 3000;
await bakana.runAnalysis(new_state, null, new_params);
```

Again, we can supply a callback that runs on the results of each step as soon as they are loaded.

```js
let reloaded = await bakana.loadAnalysis("whee.h5", loader, { finishFun: finisher });
```

If links are present, it is assumed that the application will use `setResolveLink()` to specify a mechanism to resolve each link to a data file.

## Terminating analyses

Once a particular analysis is finished, we should free the resources of its state.

```js
bakana.freeAnalysis(state);
```

If all analyses are complete, we can terminate the entire session.
This is necessary on Node.js to allow the runtime to exit properly.

```js
bakana.terminate();
```

## Developer notes

See [here](docs/related/custom_readers.md) for instructions on adding custom dataset readers.
This allows us to use **bakana**'s analysis pipeline and serialization capabilities on datasets from other sources such as in-house databases.

Testing requires some combination of the options below,
depending on the version of Node.js available.

```sh
node --experimental-vm-modules \
    --experimental-wasm-threads \
    --experimental-wasm-bulk-memory \
    --experimental-wasm-bigint \
    node_modules/jest/bin/jest.js \
    --testTimeout=100000000 \
    --runInBand
```
