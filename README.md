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
This will also report the parameters used at each step.

```js
let reloaded = await bakana.loadAnalysis("foo.h5", loader);
let new_state = reloaded.state;
let new_params = reloaded.parameters;
```

These can be used to perform a new analysis, possibly after modifying some parameters.

```js
new_params.pca.num_hvgs = 3000;
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

## Adding custom readers

Developers can add new data readers to supplement the defaults in `src/readers`.
This is most useful for pulling data from external databases for entry into the **bakana** workflow.
A new reader should be implemented as an ES6 module with the required exports below.

**`abbreviate(args)`:** generate a summary of the matrix to check whether the inputs have changed on re-analysis.
- This should accept an object `args` containing information about a single count matrix in the reader's desired format.
  `args` may contain an arbitrary number of properties - their specification is left to the reader.
  By convention, properties referring to files should be `File` objects in the browser and strings containing file paths for Node.js.
- This should quickly return an object that summarizes the information about the count matrix.
  The returned object should be easily stringified for quick comparison,
  and should be unique enough to distinguish between most other values of `args`.
  We generally recommend generating a summary based on the file name and size, which is fast and unique enough for most comparisons.

**`preflight(args)`:** fetch annotations from the matrix before the analysis begins.
- This should accept an object `args` containing information about a single count matrix in the reader's desired format.
  The expected structure of `args` is identical to that described for `abbreviate()`.
- This should return an object containing the `genes` and `annotations` properties.
  `genes` should be an object where each key is a gene annotation field name and each value is an array of per-gene information, usually strings containing Ensembl identifiers or symbols.
  `annotations` should be an array of strings containing the names of the per-cell annotation fields.
- Alternatively, this method may return a promise that resolves to such an object.

**`Reader(args)`:** the `Reader` class, where each instance contains all information about a single count matrix.
- This constructor should accept an object `args` containing information about a single count matrix in the reader's desired format.
  The expected structure of `args` is identical to that described for `abbreviate()`.
- The structure of the class is left as an implementation detail for each reader.

**`Reader.format()`:** specify the format of the count matrix.
- This should return a string specifying the format of the count matrix corresponding to this reader.
  The value of the string is defined by the reader's developer.

**`Reader.load()`:** load the count matrix into memory.
- This should return an object containing `matrix`, a `ScranMatrix` object.
  The object may also contain `genes`, an object where each key is a gene annotation field name and each value is an array of per-gene information;
  and/or `annotations`, an object where each key is a cell annotation field name and each value is an array of per-cell information.
- Alternatively, this method may return a promise that resolves to such an object.

**`Reader.serialize(embeddedSaver)`:** register the count matrix in the state file. 
- This should return an array of objects where each object represents a data file used to create the count matrix.
  Each object should contain `type`, a string specifying the type of file; and `name`, the name of the file.
  Alternatively, this method may return a promise that resolves to this array.
- If `embeddedSaver` is a function, it should be called on each file, in the same order as they occur in the returned array.
  `embeddedSaver` will accept two arguments - a file and its size - and will return a Promise that resolves to an object containing `offset` and `size`.
  (Here, a "file" is defined as an `ArrayBuffer` containing the contents of the file in the browser, or a path to a file in Node.js.)
  The reported `offset` and `size` should then be included in the properties of the corresponding object in the returned array.
- If `embeddedSaver = null`, each object in the returned array should instead contain `id`.
  This should be a unique string, typically containing a link to some external database system that holds the corresponding file.
  The interpretation of `id` is left to the reader, as long as it can be used to meaningfully recover the same file in `unserialize()` below.

**`unserialize(values, embeddedLoader)`:** load the count matrix from the state file.
- `values` should be an array of objects where object represents a data file used to create the count matrix.
  Each object should contain `type` and `name`, as described in the `Reader.serialize()` method.
- If `embeddedLoader` is a function, each object in `values` will additionally contain `offset` and `size` integers.
  `embeddedLoader` will accept two arguments - the offset and size - and will return the corresponding file.
  (Again, a "file" is defined as an `ArrayBuffer` containing the contents of the file in the browser, or a path to a file in Node.js.)
  `embeddedLoader` should be called on each entry of `values` in order - note, this should be done even if not all files are used.
- If `embeddedLoader = null`, each object in `values` will additionally contain an `id` string.
  This usually contains a link to some external database system as defined by `Reader.serialize()`.
- This method should return an instance of the `Reader` class containing all information related to the count matrix represented by `values`.
  Presumably this is done using the file contents returned by `embeddedLoader` or the linked files from the `id` values.
  Alternatively, this method may return a promise that resolves to a `Reader` instance.

Applications can then register custom readers by assigning to the `bakana.availableReaders` object.
The property name should be the name of the reader format (the same as that returned by `Reader.format()`) and the value should be the ES6 module object implementing the reader.

## Developer notes

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

