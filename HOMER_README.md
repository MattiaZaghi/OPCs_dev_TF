# HOMER Motif Analysis Guide

## Overview

This guide explains how to perform known motif enrichment analysis on your BAMPE narrowPeak files using **HOMER** (Hypergeometric Optimization of Motif EnRichment).

## What is HOMER?

HOMER is a suite of tools for ChIP-seq motif discovery and analysis. It can:
- **Find known motifs** enriched in your peaks (using databases of known transcription factor binding sites)
- **Discover de novo motifs** from your peak sequences
- **Visualize** motifs and their enrichment statistics
- **Compare motifs** across samples

## Installation

HOMER requires separate installation. Follow these steps:

### 1. Install HOMER

```bash
# Download HOMER
mkdir -p ~/homer
cd ~/homer
wget http://homer.ucsd.edu/homer/configureHomer.pl
perl configureHomer.pl -install -h hg38
```

Or use conda (if available):
```bash
conda install -c bioconda homer
```

### 2. Verify Installation

```bash
findMotifsGenome.pl --help
```

If this command works, HOMER is ready to use.

### 3. Make Script Executable

```bash
chmod +x /cfs/klemming/home/m/matzag/OPCs_dev_TF/run_homer_motifs.sh
```

## Running the Analysis

### Option 1: Local/Interactive Node

```bash
cd /cfs/klemming/home/m/matzag/OPCs_dev_TF
bash run_homer_motifs.sh
```

### Option 2: HPC Compute Node (Dardel - Recommended)

Submit as a batch job:

```bash
sbatch --time=02:00:00 --nodes=1 --cpus-per-task=8 run_homer_motifs.sh
```

Or use interactive allocation:

```bash
salloc --time=02:00:00 --nodes=1 --ntasks=1 --cpus-per-task=8
module load bioinfo-tools  # If available on your cluster
bash /cfs/klemming/home/m/matzag/OPCs_dev_TF/run_homer_motifs.sh
```

## Understanding the Output

### Output Directory Structure

```
/cfs/klemming/projects/supr/uppstore2017150/Mattia/Analysis/homer_motif_analysis/
├── index.html                          # Main results index (open in browser)
├── HOMER_analysis_summary.txt          # Text summary of analysis
├── analysis_log.txt                    # Detailed run log
│
└── [Sample_Name]_motifs/
    ├── knownResults.html               # ★ Main results (browser viewable)
    ├── knownResults.txt                # Motif statistics (tab-separated)
    ├── homerResults.html              # De novo results (if generated)
    ├── motifs/
    │   ├── motif1.logo.png            # Motif logos
    │   ├── motif1.pwm                 # Position weight matrices
    │   └── ...
    └── peaks.txt                       # BED format of input peaks
```

### Key Output Files

#### 1. **knownResults.html** (Most important!)
- Interactive HTML report showing all enriched motifs
- Displays motif logos alongside statistics
- Ranked by p-value (most significant first)
- Click motif names to see peak locations and PWM details

#### 2. **knownResults.txt**
Tab-separated file with columns:
- `Rank`: Enrichment rank
- `Motif Name`: Known motif identifier
- `Consensus`: Motif consensus sequence
- `p-value`: Statistical significance (log10 scale)
- `q-value`: False discovery rate corrected p-value
- `% of Peaks`: Percentage of input peaks containing motif
- `% of Background`: Percentage in background
- `Log Odds Score`: Log2(enrichment fold-change)

Example:
```
Rank    Motif Name              p-value     % Peaks     Log Odds
1       GABPA_RUNX3             1e-150      75.2%       8.5
2       ETS1                    1e-120      62.1%       7.2
3       NFAC                    1e-85       45.3%       6.1
```

#### 3. **motifs/ directory**
Contains:
- `.logo.png` - Visual representation of motif
- `.pwm` - Position weight matrix (can be used in other tools)
- `.motif` - HOMER format motif file

### Interpreting p-values and Enrichment

- **p-value < 0.05**: Statistically significant enrichment
- **p-value < 1e-10**: Highly significant enrichment
- **% of Peaks**: If high (>50%), the motif is common in your peaks
- **Log Odds > 3**: At least 8-fold enrichment (log2(8) ≈ 3)

## Running on HPC (Dardel)

### Create SLURM Job Script

Save as `homer_job.sh`:

```bash
#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=homer_motifs
#SBATCH --output=homer_%j.log

# Load modules if available
module load bioinfo-tools

# Or activate conda environment if HOMER is installed there
# conda activate homer_env

# Run the analysis
cd /cfs/klemming/home/m/matzag/OPCs_dev_TF
bash run_homer_motifs.sh
```

Submit:
```bash
sbatch homer_job.sh
```

Check status:
```bash
squeue -u $USER
```

View output:
```bash
tail -f homer_*.log
```

## Customizing the Analysis

### Modify Script Parameters

Edit `run_homer_motifs.sh` to adjust:

```bash
# Number of parallel processors (increase for faster analysis)
NUM_PROCS=4  # Change to 8, 16, etc.

# Reference genome (if using mm10 instead of hg38)
GENOME="hg38"  # or "mm10"

# Peak size handling
# Current: -size given (uses peak sizes from file)
# Alternative: -size 200 (fixed 200bp window)
```

### Advanced Options

For de novo motif discovery (uncomment in script):

```bash
# In findMotifsGenome.pl command, add:
# -find ~/myMotifs.motif  # Find occurrences of custom motifs
# -noknown                # Skip known motif analysis
# -nogo                   # Skip Gene Ontology analysis
```

## Comparing Motifs Across Samples

After analysis, manually compare:
1. Open multiple `knownResults.html` files in separate browser tabs
2. Compare top enriched motifs between samples
3. Note differences in motif rankings and percentages
4. Export `knownResults.txt` files to Excel for statistical comparison

## Troubleshooting

### Error: "findMotifsGenome.pl: command not found"
**Solution:** HOMER not installed or not in PATH. Install as above.

### Error: "Genome files not found"
**Solution:** Run `perl configureHomer.pl -install -h hg38` to download genome databases

### Very slow analysis
**Solution:** Increase `NUM_PROCS` in script for parallel processing

### No motifs found
**Possible causes:**
- Very few peaks in input file
- Peaks don't contain common TF binding sites
- Try `-size 200` to use fixed window instead of `-size given`

## Quality Control Checklist

- [ ] Peak files found in `/cfs/klemming/projects/supr/uppstore2017150/Mattia/macs3/BAMPE`
- [ ] HOMER installed and in PATH
- [ ] Genome database (hg38) downloaded
- [ ] Script runs without errors
- [ ] Results HTML files generate successfully
- [ ] Compare results across samples

## Alternative: De Novo Motif Discovery

If you want to discover novel motifs (not just known ones), run:

```bash
findMotifsGenome.pl /path/to/peaks.narrowPeak hg38 /path/to/output -size 200
```

This generates both `knownResults.html` (known motifs) and `homerResults.html` (de novo motifs).

Warning: De novo discovery is computationally intensive (can take hours).

## References

- **HOMER Paper**: Heinz et al. Mol Cell. 2010
- **HOMER Website**: http://homer.ucsd.edu/
- **JASPAR Database**: http://jaspar.genereg.net/ (Reference motif database)
- **TF Databases**: AnimalTFDB, ReMap, CollecTF

## Example Interpretation

### Sample Result:
```
Rank | Motif        | p-value  | % Peaks | Meaning
-----|--------------|----------|---------|------------------------------------------
1    | GABPA/RUNX3  | 1e-150   | 75%     | GABPA/RUNX3 strongly enriched (75% of peaks)
2    | ELF1         | 1e-120   | 60%     | ELF1 also significantly enriched
3    | CREB         | 1e-85    | 45%     | CREB moderately enriched
```

**Interpretation**: Your peaks are enriched for GABPA, RUNX3, and ELF1 binding sites, suggesting these transcription factors may be co-binding or that your pulldown captured regions marked by these TFs.

## Support

For questions about HOMER:
- Check http://homer.ucsd.edu/homer/ngs/findMotifs.html
- Email: homer@ucsd.edu
