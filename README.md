# Backend for kana

## Overview

**kanako** provides the compute backend for the [**kana**](https://github.com/jkanche/kana) application.
This can be used either in the browser or with Node.js.

```sh
node --experimental-vm-modules --experimental-wasm-threads --experimental-wasm-bulk-memory --experimental-wasm-bigint node_modules/jest/bin/jest.js 
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

**`Reader(args)`:** the `Reader` class, where each instance contains all information about a single count matrix.
- This constructor should accept an object `args` containing information about a single count matrix in the reader's desired format.
  The expected structure of `args` is identical to that described for `abbreviate()`.
- The structure of the class is left as an implementation detail for each reader.

**`Reader.format()`:** specify the format of the count matrix.
- This should return a string specifying the format of the count matrix corresponding to this reader.
  The value of the string is defined by the reader's developer.

**`Reader.load()`:** 
- This should return an object containing `matrix`, a `ScranMatrix` object.
  The object may also contain `genes`, an object where each key is a gene annotation field name and each value is an array of per-gene information;
  and/or `annotations`, an object where each key is a cell annotation field name and each value is an array of per-cell information.

**`Reader.serialize(embeddedSaver)`:**
- This should return an array of objects where each object represents a data file used to create the count matrix.
  Each object should contain `type`, a string specifying the type of file; and `name`, the name of the file.
- If `embeddedSaver` is a function, it should be called on each file, in the same order as they occur in the returned array.
  `embeddedSaver` will accept two arguments - a file and its size - and will return a Promise that resolves to an object containing `offset` and `size`.
  (Here, a "file" is defined as an `ArrayBuffer` containing the contents of the file in the browser, or a path to a file in Node.js.)
  The reported `offset` and `size` should then be included in the properties of the corresponding object in the returned array.
- If `embeddedSaver = null`, each object in the returned array should instead contain `id`.
  This should be a unique string, typically containing a link to some external database system that holds the corresponding file.
  The interpretation of `id` is left to the reader, as long as it can be used to meaningfully recover the same file in `unserialize()` below.
- This method should be `async`.

**`unserialize(values, embeddedLoader)`:**
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
- This method should be `async`.
