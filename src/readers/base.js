/**
 * An abstract dataset.
 * @hideconstructor
 */
export class Dataset {
    /**
     * @desc Name of the class.
     * @type {string}
     */
    static className = "Dataset";

    /**
     * @return {string} Name of the dataset format.
     */
    static format() {
        throw new Error("no known 'format' method for the base class");
    }

    /**
     * @return {Object} Object containing some unique-ish summary of the files in this dataset,
     * to be used as a key for state comparisons.
     * (Note that this only needs to be reasonably unique compared to other instances of the same dataset class.)
     */
    async abbreviate() {
        throw new Error("no available 'abbreviate' method");
    }

    /**
     * @return {Object} Object containing the following: 
     *
     * - `features`, an object where keys are the names of modalities and values contain the feature annotation for that modality.
     *   Each value is itself an object where the keys are the names of the annotation fields and the values are arrays of per-feature annotations.
     *   Each array should be of length equal to the number of features in its modalities.
     * - `cells`, an object whe keys are the name of the per-cell annotation fields and the values contain per-cell summaries.
     *   Each summary should have the same format as {@linkcode summarizeArray}.
     */
    async annotations() {
        throw new Error("no available 'annotations' method");
    }

    /**
     * @return {Object} Object containing the following: 
     *
     * - `matrix`, a [scran.MultiMatrix](https://www.jkanche.com/scran.js/MultiMatrix) object containing one [ScranMatrix](https://www.jkanche.com/scran.js/ScranMatrix) object per modality.
     *   Each ScranMatrix should have number of rows equal to the number of features in that modality. 
     *   The columns of each ScranMatrix should correspond to the same ordering of cells for each modality.
     * - `row_ids`, an object where keys are the names of modalities and values are Array/TypedArray of integers of length equal to the number of rows in the corresponding entry of `matrix`.
     *   Each integer is unique within its modality and represents some feature identity for corresponding row.
     *   The interpretation of the identity for each feature is left to the implementation but should not change across versions as it is used for serialization/unserialization.
     * - `features`, an object where keys are the names of modalities and values contain the feature annotation for that modality.
     *   Each value is itself an object where the keys are the names of the annotation fields and the values are arrays of per-feature annotations.
     *   Elements of each array should correspond to rows of the relevant `matrix` entry.
     * - `cells`, an object whe keys are the name of the per-cell annotation fields and the values contain arrays of per-cell annotations.
     *   Each array should be of length equal to the number of columns in each entry of `matrix`.
     *
     * Note that some reordering of the `features` arrays may be necessary if the ScranMatrix permutes its rows to [improve memory efficiency](https://github.com/jkanche/scran.js/blob/master/docs/related/row_permutations.md).
     * Of course, if no reordering is performed, the output of {@linkcode Dataset#features Dataset.features} can be returned as-is.
     *
     * @async 
     */
    async load() {
        throw new Error("no available 'load' method");
    }

    /**
     * Free up resources associated with this dataset.
     * @async
     */
    async free() {}

    /**
     * @return {Array} Array of objects where each object represents a file used in this dataset.
     * Each object should contain:
     *
     * - `type`: a string specifying the file type.
     * - `file`: a {@linkplain SimpleFile} object representing the file contents.
     *
     * @async
     */
    async serialize() {
        throw new Error("no available 'serialize' method");
    }

    /**
     * @param {Array} Array of objects where each object represents a file.
     * This should be equivalent to the output of {@linkcode Dataset#serialize Dataset.serialize} from the same class.
     * @return An instance of a `Dataset` of the same subclass.
     * @async
     */
    static async unserialize(files) {

    }
}
