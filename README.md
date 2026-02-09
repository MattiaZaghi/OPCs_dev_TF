# OPCs_dev_TF

Snakemake pipeline to process Cut&Tag / Cut&Run-like data (adapted from TbL_mintbodies).

Quick start

1. Edit `snakemake/config_CutTag.yaml` and set `RUN_ID`, `SAMPLES_JSON`, and genome/index paths.
2. Provide a samples JSON file (see `snakemake/samples_nanoCT.json`).
3. Run:

```bash
cd OPCs_dev_TF/snakemake
 snakemake --snakefile /home/mattia/OPCs_dev_TF/snakemake/snakefile.smk  --cores 100 --jobs 100 --profile htcondor -p --use-conda --configfile /home/mattia/OPCs_dev_TF/snakemake/config_CutTag.yaml --rerun-incomplete
```

Notes
- The provided conda env YAMLs are minimal placeholders. Adjust versions and channels to match your cluster.
- Paths in the config are relative or local; set absolute paths where needed.

Generating `samples.json`

The repository includes `scripts/sample2json.py`, a small helper that scans a directory of FASTQ files and produces a JSON mapping used by the Snakemake workflow.

How the script actually works
- It recursively searches `--fastq_dir` for files ending with `.gz` and collects their full paths.
- It reads a tab-delimited metadata file provided with `--meta` and expects at least six columns per row (no strict header parsing):
	1) sample name
	2) fastq file identifier (a substring expected to appear in the FASTQ filename)
	3) factor (e.g. a histone mark or TF)
	4) assay (e.g. ChIP, Input)
	5) sample type
	6) reference (e.g. genome)
- For each metadata row the script looks for FASTQ paths where the `fastq_name` string appears in the FASTQ path. All matching `.gz` paths are added to the output under `FILES[sample][factor][assay]`.
- The script does not require strict filename patterns (it uses substring matching), so the `fastq_name` in the metadata must match a portion of the actual FASTQ filename.
- Output is written as JSON to `samples_ChIP.json` in the current working directory (this filename is fixed in the script).

Example metadata (tab-separated, header optional):

```
sample	fastq_name	factor	assay	sample_type	reference
SampleA	SampleA_R1	H3K27ac	ChIP	IP	hg38
```

Example usage:

```bash
python3 scripts/sample2json.py --fastq_dir /full/path/to/fastqs --meta /full/path/to/meta.txt
# -> writes ./samples_ChIP.json
```

Concrete example from this workspace
-----------------------------------
Here is a real example metadata file found in the repository (`OPCs_dev_TF/snakemake/meta.txt`) — you can paste a similar table for your own data:

```
sample_name	fastq_name	factor	assay	type	reference
OPC	NM_S	SOX10	nanoCT	narrow	human
OPC	NM_SL	LUZP2	nanoCT	narrow	human
OPC	NM_ST	TSC22	nanoCT	narrow	human
```

Notes on this example:
- The second column (`fastq_name`) is a substring the script will search for inside actual FASTQ filenames. For example, a FASTQ file named `NM_S_R1_001.fastq.gz` will match `e15_nanoCTR_1`.
- The script will add any `.gz` file whose path contains that substring (both R1 and R2 if present).
- Adapt the `sample_name` and `fastq_name` columns for your dataset; keep the tab-separated format.


Notes & troubleshooting
- The script matches `fastq_name` by substring; if you see "missing" messages, check that the metadata `fastq_name` matches the actual FASTQ filenames (or include a unique substring).
- The script collects any `.gz` files that contain the `fastq_name` substring (so it will collect both R1 and R2 if both contain that substring).
- If you want the JSON in a different location or filename, move/rename `samples_ChIP.json` after the script completes (or modify the script to accept an output path).
- The JSON structure produced is nested as: sample -> factor -> assay -> [list of fastq full paths]. This is compatible with the Snakefile's expected samples mapping after minor adjustments.

Suggested improvements (optional):
- Add a `--out` argument to make the output filename configurable.
- Add stricter validation of the metadata header and a mode to match exact filenames instead of substring matching.



