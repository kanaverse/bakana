base=https://github.com/jkanche/random-test-files/releases/download

# Download PBMC datasets.
curl -L ${base}/10x-pbmc-v1.0.0/pbmc3k-matrix.mtx.gz > pbmc3k-matrix.mtx.gz
curl -L ${base}/10x-pbmc-v1.0.0/pbmc3k-features.tsv.gz > pbmc3k-features.tsv.gz
curl -L ${base}/10x-pbmc-v1.0.0/pbmc3k-barcodes.tsv.gz > pbmc3k-barcodes.tsv.gz
curl -L ${base}/10x-pbmc-v1.0.0/pbmc4k-tenx.h5 > pbmc4k-tenx.h5

# Download Zeisel datasets.
curl -L ${base}/zeisel-brain-v1.0.0/csc.h5ad > zeisel-brain.h5ad

