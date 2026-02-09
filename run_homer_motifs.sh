#!/bin/bash

# ============================================================================
# HOMER Motif Analysis Script for BAMPE narrowPeak Files
# ============================================================================
# This script performs known motif enrichment analysis on ChIP-seq peaks
# using HOMER (Hypergeometric Optimization of Motif EnRichment)
#
# Usage: bash run_homer_motifs.sh
#
# Requirements:
# - HOMER suite installed (http://homer.ucsd.edu/)
# - Genome database files (downloaded via configureHomer.pl)
# - hg38 peak files in narrowPeak format

set -e  # Exit on error

# ============================================================================
# Configuration
# ============================================================================

# HOMER installation directory (local installation)
HOMER_HOME="/cfs/klemming/home/m/matzag/homer"
export PATH="${HOMER_HOME}/bin:$PATH"

# Input and output directories
PEAK_DIR="/cfs/klemming/projects/supr/uppstore2017150/Mattia/macs3/BAMPE"
OUT_BASE_DIR="/cfs/klemming/projects/supr/uppstore2017150/Mattia/Analysis"
MOTIF_OUT_DIR="${OUT_BASE_DIR}/homer_motif_analysis"

# Create output directory
mkdir -p "$MOTIF_OUT_DIR"

# Reference genome
GENOME="hg38"

# Number of processors for parallel analysis
NUM_PROCS=4

# ============================================================================
# Check HOMER Installation
# ============================================================================

echo "Checking HOMER installation..."
echo "HOMER_HOME: $HOMER_HOME"

if [ ! -f "${HOMER_HOME}/bin/findMotifsGenome.pl" ]; then
    echo "ERROR: HOMER not found at $HOMER_HOME"
    echo "Expected location: ${HOMER_HOME}/bin/findMotifsGenome.pl"
    echo "Please check the HOMER installation path"
    exit 1
fi

echo "✓ HOMER installation found at: $HOMER_HOME"
echo "  findMotifsGenome.pl: ${HOMER_HOME}/bin/findMotifsGenome.pl"

# ============================================================================
# Verify Input Files Exist
# ============================================================================

echo ""
echo "Looking for narrowPeak files in: $PEAK_DIR"

PEAK_FILES=($(find "$PEAK_DIR" -maxdepth 1 -name "*.narrowPeak" -o -name "*_peaks.narrowPeak"))

if [ ${#PEAK_FILES[@]} -eq 0 ]; then
    echo "ERROR: No narrowPeak files found in $PEAK_DIR"
    exit 1
fi

echo "Found ${#PEAK_FILES[@]} peak file(s):"
for pf in "${PEAK_FILES[@]}"; do
    echo "  - $(basename "$pf")"
done

# ============================================================================
# Run HOMER Motif Analysis for Each Peak File
# ============================================================================

echo ""
echo "Starting HOMER motif analysis..."
echo "=============================================="

ANALYSIS_LOG="${MOTIF_OUT_DIR}/analysis_log.txt"
> "$ANALYSIS_LOG"  # Clear log file

for peak_file in "${PEAK_FILES[@]}"; do
    # Extract sample name
    sample_name=$(basename "$peak_file" | sed 's/_peaks\.narrowPeak$//; s/\.narrowPeak$//')
    
    # Create sample-specific output directory
    sample_out_dir="${MOTIF_OUT_DIR}/${sample_name}_motifs"
    mkdir -p "$sample_out_dir"
    
    echo ""
    echo "Processing: $sample_name"
    echo "Input: $peak_file"
    echo "Output: $sample_out_dir"
    echo "---"
    
    {
        echo "Sample: $sample_name"
        echo "Peak file: $peak_file"
        echo "Start time: $(date)"
        echo ""
    } >> "$ANALYSIS_LOG"
    
    # Run findMotifsGenome.pl with known motif enrichment
    # Options:
    # -size given       : Use peak sizes from the file (not fixed window)
    # -p NUM_PROCS      : Use parallel processing
    # -h                : Find Hypergeometric enrichment (uses built-in motif database)
    
    if "${HOMER_HOME}/bin/findMotifsGenome.pl" "$peak_file" "$GENOME" "$sample_out_dir" \
        -size given \
        -p "$NUM_PROCS" \
        -h 2>&1 | tee -a "$ANALYSIS_LOG"; then
        
        echo "✓ Successfully completed: $sample_name" | tee -a "$ANALYSIS_LOG"
        
        # Check if key output files were generated
        if [ -f "${sample_out_dir}/knownResults.html" ]; then
            echo "  Output files:" | tee -a "$ANALYSIS_LOG"
            echo "    - knownResults.html: Known motif enrichment results" | tee -a "$ANALYSIS_LOG"
            echo "    - knownResults.txt: Motif statistics (tab-separated)" | tee -a "$ANALYSIS_LOG"
            echo "    - motifs/ directory: Motif logos and details" | tee -a "$ANALYSIS_LOG"
        fi
    else
        echo "✗ Failed to analyze: $sample_name" | tee -a "$ANALYSIS_LOG"
        echo "  Check ${sample_out_dir}/error log for details" | tee -a "$ANALYSIS_LOG"
        continue
    fi
    
    echo "" >> "$ANALYSIS_LOG"
done

# ============================================================================
# Create Summary Report
# ============================================================================

echo ""
echo "=============================================="
echo "Generating summary report..."

SUMMARY_FILE="${MOTIF_OUT_DIR}/HOMER_analysis_summary.txt"

cat > "$SUMMARY_FILE" << 'EOF'
================================================================================
HOMER MOTIF ANALYSIS SUMMARY
================================================================================

Analysis Date: $(date)

INPUT PARAMETERS:
  Peak directory: /cfs/klemming/projects/supr/uppstore2017150/Mattia/macs3/BAMPE
  Reference genome: hg38
  Analysis type: Known motif enrichment (Hypergeometric)
  Output directory: /cfs/klemming/projects/supr/uppstore2017150/Mattia/Analysis/homer_motif_analysis

OUTPUT FILES FOR EACH SAMPLE:
  1. knownResults.html
     - Interactive HTML report of known motif enrichment
     - Shows motif logos, p-values, and peak percentages
     - Ranked by enrichment significance

  2. knownResults.txt
     - Tab-separated file with detailed motif statistics:
       * Motif name and consensus sequence
       * p-value and q-value (adjusted p-value)
       * Percentage of peaks containing the motif
       * Log odds score (fold enrichment)

  3. motifs/ directory
     - Individual motif logo files (PNG)
     - Motif consensus sequences
     - Detailed motif information

INTERPRETATION:
  - p-value: Statistical significance of motif enrichment
  - % of Peaks: Percentage of input peaks containing the motif
  - Log odds: Log2 fold-change compared to background
  - Lower p-values indicate more significant enrichments

KEY MOTIFS:
  Check the knownResults.html files (browser viewable) for:
  - Most enriched known motifs in your peaks
  - Ranked by statistical significance
  - Comparison across different TF samples

NEXT STEPS:
  1. Open knownResults.html in a web browser to explore motif enrichment
  2. Compare motif enrichment patterns across samples
  3. Validate top motifs against literature/databases
  4. Consider de novo motif discovery (generateMotifLogos.pl)

NOTES:
  - HOMER uses position weight matrices (PWMs) from established databases
  - Known motif database includes vertebrate, yeast, and other organisms
  - Peak sizes are preserved from narrowPeak files (-size given option)
  - Motif enrichment tested against background regions of equal size

For more information: http://homer.ucsd.edu/homer/ngs/findMotifs.html
================================================================================
EOF

echo "Summary saved to: $SUMMARY_FILE"

# ============================================================================
# Create Index HTML File
# ============================================================================

INDEX_HTML="${MOTIF_OUT_DIR}/index.html"

cat > "$INDEX_HTML" << 'EOF'
<!DOCTYPE html>
<html>
<head>
    <title>HOMER Motif Analysis Results</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        h1 { color: #333; border-bottom: 2px solid #007bff; }
        h2 { color: #555; margin-top: 30px; }
        .sample { background-color: #f9f9f9; padding: 15px; margin: 10px 0; border-left: 4px solid #007bff; }
        a { color: #007bff; text-decoration: none; }
        a:hover { text-decoration: underline; }
        table { border-collapse: collapse; width: 100%; margin: 20px 0; }
        th, td { border: 1px solid #ddd; padding: 10px; text-align: left; }
        th { background-color: #f2f2f2; }
        .note { background-color: #fff3cd; padding: 10px; margin: 10px 0; border-left: 4px solid #ffc107; }
    </style>
</head>
<body>
    <h1>HOMER Known Motif Enrichment Analysis</h1>
    <p><strong>Date:</strong> <span id="date"></span></p>

    <div class="note">
        <strong>📌 Note:</strong> Each sample directory contains <code>knownResults.html</code> 
        which provides the interactive motif enrichment results. Open that file in a web browser 
        to explore the enriched motifs.
    </div>

    <h2>Sample Analysis Results</h2>
    <table>
        <tr>
            <th>Sample</th>
            <th>Results</th>
            <th>Description</th>
        </tr>
EOF

# Find all knownResults.html files and add them to the index
for result_html in $(find "$MOTIF_OUT_DIR" -name "knownResults.html" -type f); do
    sample_dir=$(dirname "$result_html")
    sample_name=$(basename "$sample_dir" | sed 's/_motifs$//')
    
    echo "<tr>" >> "$INDEX_HTML"
    echo "  <td><strong>$sample_name</strong></td>" >> "$INDEX_HTML"
    echo "  <td><a href=\"${sample_name}_motifs/knownResults.html\" target=\"_blank\">📊 View Results</a></td>" >> "$INDEX_HTML"
    echo "  <td>Known motif enrichment analysis</td>" >> "$INDEX_HTML"
    echo "</tr>" >> "$INDEX_HTML"
done

cat >> "$INDEX_HTML" << 'EOF'
    </table>

    <h2>Output Files for Each Sample</h2>
    <ul>
        <li><strong>knownResults.html</strong> - Interactive report (open in browser)</li>
        <li><strong>knownResults.txt</strong> - Tab-separated motif statistics</li>
        <li><strong>motifs/</strong> - Directory with motif logos and PWM files</li>
        <li><strong>homerResults.html</strong> - De novo motif discovery results (if generated)</li>
    </ul>

    <h2>Understanding the Results</h2>
    <dl>
        <dt><strong>Motif Name</strong></dt>
        <dd>Name and ID of the known motif from HOMER database</dd>
        
        <dt><strong>p-value</strong></dt>
        <dd>Statistical significance (lower = more significant)</dd>
        
        <dt><strong>% of Peaks</strong></dt>
        <dd>Percentage of input peaks containing the motif</dd>
        
        <dt><strong>Log Odds Score</strong></dt>
        <dd>Log2 fold-enrichment compared to background sequences</dd>
    </dl>

    <h2>Resources</h2>
    <ul>
        <li><a href="http://homer.ucsd.edu/homer/ngs/findMotifs.html" target="_blank">HOMER Documentation</a></li>
        <li><a href="http://jaspar.genereg.net/" target="_blank">JASPAR Motif Database</a></li>
        <li><a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2925483/" target="_blank">HOMER Publication</a></li>
    </ul>

    <script>
        document.getElementById('date').textContent = new Date().toLocaleString();
    </script>
</body>
</html>
EOF

echo "Index HTML created: $INDEX_HTML"

# ============================================================================
# Final Summary
# ============================================================================

echo ""
echo "=============================================="
echo "✓ HOMER motif analysis complete!"
echo "=============================================="
echo ""
echo "Output directory: $MOTIF_OUT_DIR"
echo ""
echo "Key output files:"
echo "  - index.html: Start here! (open in web browser)"
echo "  - [sample]_motifs/knownResults.html: Individual sample results"
echo "  - [sample]_motifs/knownResults.txt: Motif statistics (Excel compatible)"
echo ""
echo "Analysis log: $ANALYSIS_LOG"
echo ""
echo "Next steps:"
echo "  1. Open ${MOTIF_OUT_DIR}/index.html in a web browser"
echo "  2. Click on sample results to view enriched motifs"
echo "  3. Compare motif patterns across samples"
echo ""
