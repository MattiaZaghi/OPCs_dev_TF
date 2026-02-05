# OPCs_dev_TF

Snakemake pipeline to process Cut&Tag / Cut&Run-like data (adapted from TbL_mintbodies).

Quick start

1. Edit `snakemake/config_CutTag.yaml` and set `RUN_ID`, `SAMPLES_JSON`, and genome/index paths.
2. Provide a samples JSON file (see `snakemake/samples_example.json`).
3. Create a conda environment or use `mamba`/`conda` to create envs from `envs/*.yml`.
4. Run:

```bash
cd Private/OPCs_dev_TF/snakemake
snakemake -s snakefile.smk --configfile config_CutTag.yaml --cores 8
```

Notes
- The provided conda env YAMLs are minimal placeholders. Adjust versions and channels to match your cluster.
- Paths in the config are relative or local; set absolute paths where needed.

Generating `samples.json`

Use `scripts/sample2json.py` to create a `samples.json` from a directory of FASTQ files named like `sample_sampleType_assay_R1.fastq.gz`.

Example:

```bash
python3 scripts/sample2json.py -i ../data -o snakemake/samples.json --relative-to snakemake
```


License: MIT
