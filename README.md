# rdeval-manuscript
Scripts used for the rdeval manuscript.

## Compute raw reads summary statistics for the Vertebrate Genomes Project
The following commands show how to access SRA accession numbers and metadata from NCBI, download the raw reads and generate summary statistics with rdeval. The process is heavily parallelized, with approximately 4 cores used by each rdeval process, and 8 experiments processed in parallel (32 cores total). SRA samples can optionally be randomly sampled.

First, we can use [NCBI datasets]((https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/) and [jq](https://jqlang.github.io/jq/) to parse the VGP umbrella BioProject and collect all SRA accessions associated with its genomes (here is jq's [manual](https://jqlang.github.io/jq/manual/)):

```
datasets summary genome accession PRJNA489243 > vgp-metadata.json
cat vgp-metadata.json | jq -r '.reports[] | .accession + "," + .assembly_info.assembly_name + "," + (.assembly_info.biosample.sample_ids[] | select(.db=="SRA").value)' > raw_data_metadata.ls
cat raw_data_metadata.ls
```

Then we can download the data (requires `fasterq-dump` in `$PATH`) and compute summary statistics by generating .rd files (requires `rdeval` in `$PATH`):

```
bash rdeval_parallel.sh raw_data_metadata.ls 0.5 // sampling fraction, 1 for full data set
```