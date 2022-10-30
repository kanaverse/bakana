# Adding custom readers

## Overview

Developers can add new data readers to supplement the defaults in `src/readers`.
This is most useful for pulling data from external databases for entry into the **bakana** workflow.
A new reader should be implemented as an ES6 class that satisfies the `Dataset` interface requirements below.

## `Dataset` interface

### Constructor

This is called by users before entry into the **bakana** workflow, and so is not subject to any particular requirements.
For consistency with the existing `Dataset` readers, developers of new file-based readers should consider accepting [`SimpleFile`](https://ltla.github.io/bakana/SimpleFile.html) objects.
This provides a convenient abstraction for files in Node.js and browser contexts.

### `format()` (static)

This is a static method that returns a string containing the name of the dataset format.
It is used to differentiate between different `Dataset` classes.
The formats implemented in _bakana_ itself are listed below:

- `10X`: for the 10X HDF5 format in `TenxHdf5Dataset`.
- `MatrixMarket`: for the 10X Matrix Market format in `TenxMatrixMarketDataset`.
- `H5AD`: for the H5AD format `H5adDataset`.
- `SummarizedExperiment`: for SummarizedExperiments saved as RDS files in `SummarizedExperimentDataset`.

### `abbreviate()`

This method should _quickly_ return an object containing some unique-ish summary of the files in this dataset.
The object will be used for comparisons between states in the `InputsState.compute()` function.
It should be easily stringified while still being unique enough to distinguish between most other instances of the same class. 
We generally recommend generating a summary based on the file name and size, which is fast and unique enough for most comparisons.
(Note that this only needs to be reasonably unique compared to other instances of the same `Dataset` class.)

### `annotations()`

This method should return an object containing a summary of the annotations in this dataset.
The object should contain:

- `features`: an object where each key is the name of a modality (e.g., RNA, ADT) and each value is a [`DataFrame`](https://ltla.github.io/bioconductor.js/DataFrame.html).
  Each `DataFrame` should contain one row per feature in that modality, along with any number of columns containing the per-feature annotation fields.
- `cells`: an object containing:
  - `number`: the number of cells in the dataset.
  - `summary`: an object where each key is the name of a per-cell annotation field and each value is an object containing a summary of that field.
    Each summary object should contain:
    - a `type` property, indicating whether the annotation is `"categorical"` or "`continuous"`.
    - for categorical annotations, a `values` array containing the unique values.
      This may be truncated for brevity, in which case the `truncated` property is `true`.
    - for continuous annotations, the `min` and `max` properties containing the minimum and maximum values.

Alternatively, this method may return a promise that resolves to such an object.

### `load()`

This method should return an object containing:

- `matrix`, a [`MultiMatrix`](https://jkanche.github.io/scran.js/MultiMatrix.html) object containing submatrices for at least one modality.
  Each modality-specific submatrix should be a [`ScranMatrix`](https://jkanche.github.io/scran.js/ScranMatrix.html) containing any number of rows.
  All modalities should contain data for the same number of columns.
- `row_ids`, an object specifying the row identities for each submatrix in `matrix`.
  Keys should correspond to the modality names in `matrix` and each value should be an integer array containing the identities of the submatrix rows.
  Feature identities should be encoded as non-negative integers that are unique within each submatrix.
  Developers should endeavor to keep identities consistent across different versions of the reader (i.e., the same feature gets the same integer identity),
  as **bakana** will automatically use this information to transform old results to match the row order in the current version of the application.
- `features`, an object where each key is the name of a modality in `matrix`.
  Each modality-specific value is a `DataFrame` with one row per feature in the corresponding entry of `matrix`.
  Columns should be per-feature annotation fields, as described for `annotations()`.
  If the rows of `matrix.get(<modality>)` were reorganized in any way, e.g., subsetting or permutation, the same reorganization should be applied to the rows of `features[<modality>]`.
- `cells`, a `DataFrame` containing one row per cell.   
  Each column should contain an array of per-cell information, corresponding to the same order of columns in each entry of `matrix`.

Alternatively, this method may return a promise that resolves to such an object.

### `serialize()`

This method should return an array of objects, where each object represents a data file used to create the count matrix.
Each object should contain:

- `type`, a string specifying the type of file.
  This is used to disambiguate multiple files for a single dataset.
- `file`, a `SimpleFile` object containing the contents of (or a reference to) the underlying data file.
  For database-based datasets, developers may wish to mock up a file containing the database identifier, e.g., as a JSON string. 

### `unserialize(files)` (static)

This is a static method that creates a new instance of the `Dataset`.

`files` will be an array of objects where object represents a data file used to create the dataset.
Each object will contain `type` and `file` as described for the `serialize()` method.

This method should return an instance of the `Dataset` class.
Presumably this is done using the file contents in `files`.
Alternatively, this method may return a promise that resolves to a `Reader` instance.

## Registering new readers

Applications can register custom readers by assigning to the `bakana.availableReaders` object.
The property name should be the name of the dataset format (i.e., the value of `format()`) and the value should be the class definition.
For example, one could re-register the 10X HDF5 reader:

```js
bakana.availableReaders[bakana.TenxHdf5Dataset.format()] = bakana.TenxHdf5Dataset;
```
