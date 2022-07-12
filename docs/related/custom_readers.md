# Adding custom readers

## Overview

Developers can add new data readers to supplement the defaults in `src/readers`.
This is most useful for pulling data from external databases for entry into the **bakana** workflow.
A new reader should be implemented as an ES6 module with the required exports below.

## Exports

```js
abbreviate(args) 
```

**Generate a summary of the matrix to check whether the inputs have changed on re-analysis.**

This function should accept an object `args` containing information about a single count matrix in the reader's desired format.
`args` may contain an arbitrary number of properties - their specification is left to the reader.
By convention, properties referring to files should be `File` objects in the browser and strings containing file paths for Node.js.

This function should _quickly_ return an object that summarizes the information about the count matrix.
The returned object should be easily stringified for quick comparison,
while still being unique enough to distinguish between most other values of `args`.
We generally recommend generating a summary based on the file name and size, which is fast and unique enough for most comparisons.

```js
preflight(args) 
```

**Fetch annotations from the matrix before the analysis begins.**

This function should accept an object `args` containing information about a single count matrix in the reader's desired format.
The expected structure of `args` is identical to that described for `abbreviate()`.

It should return an object containing the `genes` and `annotations` properties:

- `genes` should be an object where each key is a modality, typically `"RNA"` or `"ADT"`.
  Each modality-specific value should itself be an object where each key is a gene annotation field name and each value is an array of per-gene information.
  The latter usually contains strings with Ensembl identifiers or symbols.
- `annotations` should be an object describing the per-cell annotation fields.
  Each key should be the name of a field, and each value should be an object summarizing the contents of the field.
  Specifically, each inner object should contain:

  - a `type` property, indicating whether the annotation is `"categorical"` or "`continuous"`.
  - for categorical annotations, a `values` array containing the unique values.
    This may be truncated for brevity, in which case the `truncated` property is `true`.
  - for continuous annotations, the `min` and `max` properties containing the minimum and maximum values.

  For preflight purposes, not all annotation columns need to be listed here.

Alternatively, this method may return a promise that resolves to such an object.

```js
new Reader(args)
```

**Represent all information about a single count matrix.**

This constructor should accept an object `args` containing information about a single count matrix in the reader's desired format.
The expected structure of `args` is identical to that described for `abbreviate()`.

The structure of the class is left as an implementation detail for each reader.

```js
Reader.prototype.format()
```

**Specify the format of the count matrix.**

This should return a string specifying the format of the count matrix corresponding to this reader.
The value of the string is defined by the reader's developer.

```js
Reader.protoype.load()
```

**Load the count matrix into memory.**

This should return an object containing:

- `matrix`, a `MultiMatrix` object containing matrices for at least one modality.
  Each modality-specific matrix should be a `ScranMatrix` containing any number of rows.
  All modalities should contain data for the same number of columns.
- The object may also contain `gene`, an object where each key is the name of a modality in `matrix`.
  Each modality-specific value is another object where each keys is a gene annotation field name and each value is an array of per-gene information.
  Each array corresponds to the matrix of the corresponding modality in `matrix`;
  thus, if the rows of `matrix` were reorganized in any way (e.g., subsetting or permutation), the same reorganization should be applied to each array in `genes[<modality>]`.
- The object may also contain `annotations`, an object where each key is a cell annotation field name and each value is an array of per-cell information.
  Each array should correspond to a column in `matrix`.

Alternatively, this method may return a promise that resolves to such an object.

```js
Reader.prototype.serialize(embeddedSaver)
```

**Register the count matrix in the state file.**

This should return an array of objects where each object represents a data file used to create the count matrix.
Each object should contain `type`, a string specifying the type of file; and `name`, the name of the file.
Alternatively, this method may return a promise that resolves to this array.

If `embeddedSaver` is a function, it should be called on each file, in the same order as they occur in the returned array.
`embeddedSaver` will accept two arguments - a file and its size - and will return a Promise that resolves to an object containing `offset` and `size`.
(Here, a "file" is defined as an `ArrayBuffer` containing the contents of the file in the browser, or a path to a file in Node.js.)
The reported `offset` and `size` should then be included in the properties of the corresponding object in the returned array.

If `embeddedSaver = null`, each object in the returned array should instead contain `id`.
This should be a unique string, typically containing a link to some external database system that holds the corresponding file.
The interpretation of `id` is left to the reader, as long as it can be used to meaningfully recover the same file in `unserialize()` below.

```js
unserialize(values, embeddedLoader)
```

**Load the count matrix from the state file.**

`values` should be an array of objects where object represents a data file used to create the count matrix.
Each object should contain `type` and `name`, as described in the `Reader.serialize()` method.

If `embeddedLoader` is a function, each object in `values` will additionally contain `offset` and `size` integers.
`embeddedLoader` will accept two arguments - the offset and size - and will return the corresponding file.
(Again, a "file" is defined as an `ArrayBuffer` containing the contents of the file in the browser, or a path to a file in Node.js.)
`embeddedLoader` should be called on each entry of `values` in order - note, this should be done even if not all files are used.

If `embeddedLoader = null`, each object in `values` will additionally contain an `id` string.
This usually contains a link to some external database system as defined by `Reader.serialize()`.

This method should return an instance of the `Reader` class containing all information related to the count matrix represented by `values`.
Presumably this is done using the file contents returned by `embeddedLoader` or the linked files from the `id` values.
Alternatively, this method may return a promise that resolves to a `Reader` instance.

## Registering new readers

Applications can register custom readers by assigning to the `bakana.availableReaders` object.
The property name should be the name of the reader format (i.e., the value of `Reader.prototype.format()`) and the value should be the ES6 module object implementing the reader.
