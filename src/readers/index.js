export * from "./10x.js";
export * from "./h5ad.js";
export * from "./mtx.js"
export * from "./se.js";
export * from "./ArtifactDB-abstract.js";
export * from "./ArtifactDB-zipped.js";
export * from "./utils/extract.js";
export * from "./abstract/file.js";

import { TenxHdf5Dataset } from "./10x.js";
import { H5adDataset } from "./h5ad.js";
import { TenxMatrixMarketDataset } from "./mtx.js"
import { SummarizedExperimentDataset } from "./se.js";
import { ZippedArtifactdbDataset } from "./ArtifactDB-zipped.js";

/**
 * Any class that satisfies the [Dataset contract](https://github.com/LTLA/bakana/blob/master/docs/related/custom_readers.md).
 * Each class contains methods to load data from some arbitrary data source into {@linkplain ScranMatrix} objects (for the counts)
 * and {@linkplain DataFrame} objects (for the feature or cell annotations).
 * The default set of known dataset reader classes is listed in the {@linkcode availableReaders} object
 * and includes {@linkplain TenxHdf5Dataset}, {@linkplain TenxMatrixMarketDataset}, {@linkplain H5adDataset} and {@linkplain SummarizedExperimentDataset} instances.
 *
 * @typedef Dataset
 */

/**
 * A representation of a matrix of expression values, where the values are hosted on the Wasm heap for easier compute via [**scran.js**](https://github.com/kanaverse/scran.js).
 * See [here](https://kanaverse.github.io/scran.js/ScranMatrix.html) for more details.
 *
 * @external ScranMatrix
 */ 

/**
 * A representation of multiple {@linkplain external:ScranMatrix ScranMatrix} objects, where each object contains data for the same cells but across a different feature space, e.g., for different data modalities.
 * See [here](https://kanaverse.github.io/scran.js/MultiMatrix.html) for more details.
 *
 * @external MultiMatrix
 */ 

/**
 * A DataFrame from the [**bioconductor**](https://github.com/LTLA/bioconductor.js) package, where each column is represented by some arbitrary vector-like object.
 * See [here](https://ltla.github.io/bioconductor.js/DataFrame.html) for more details.
 *
 * @external DataFrame
 */ 

/**
 * All known dataset readers.
 * Each entry contains a {@linkplain Dataset} class with the key defined as the {@linkcode Dataset#format format} return value.
 */
export const availableReaders = {
    "10X": TenxHdf5Dataset,
    "MatrixMarket": TenxMatrixMarketDataset,
    "H5AD": H5adDataset,
    "SummarizedExperiment": SummarizedExperimentDataset,
    "ArtifactDB-zipped": ZippedArtifactdbDataset
};
