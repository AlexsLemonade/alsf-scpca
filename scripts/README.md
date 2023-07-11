# Scripts for adding and modifying `scpca-meta.json`

The files in this folder were utilized to create or adjust the `scpca-meta.json` file output by running `scpca-nf`. Specifically, these files were used to synchronize changes in `scpca-nf` with the metadata output, which allows for skipping mapping when re-processing projects through `scpca-nf`.

1. `update_scpca_checkpoints.py`: This script was used to generate the original `scpca-meta.json` files for all libraries that had been processed before `scpca-nf v0.4.0`.
In this script, the old results from mapping are copied over from the `internal` to `checkpoints` directory.
Additionally, the `scpca-meta.json` file was created and saved to the `checkpoints` directory for each run to track processing information.

2. `add-refs-scpca-meta.py`: This script specifically updates existing `scpca-meta.json` files within the `checkpoints` directory.
For any libraries processed through `scpca-nf v0.5.1` or earlier, the `scpca-meta.json` file is modified to include two additional fields - `ref_mito` and `ref_fasta_index`.
