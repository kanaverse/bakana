export * from "./10x.js";
export * from "./h5ad.js";
export * from "./mtx.js"
export * from "./se.js";
export * from "./base.js";
export * from "./utils/extract.js";
export * from "./abstract/file.js";

import { TenxHdf5Dataset } from "./10x.js";
import { H5adDataset } from "./h5ad.js";
import { TenxMatrixMarketDataset } from "./mtx.js"
import { SummarizedExperimentDataset } from "./se.js";

/**
 * All known readers.
 * Each entry contains a {@linkplain Dataset} class with the key defined as the {@linkcode Dataset#format format} return value.
 */
export const availableReaders = {
    "10X": TenxHdf5Dataset,
    "MatrixMarket": TenxMatrixMarketDataset,
    "H5AD": H5adDataset,
    "SummarizedExperiment": SummarizedExperimentDataset
};
