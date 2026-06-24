---
name: differential-analysis
description: >
  Differential expression analysis comparing case vs control groups
  from a gene expression matrix. Supports DESeq2, limma, edgeR,
  t-test, and Wilcoxon methods. Produces volcano plot, heatmap,
  and optionally Venn diagram and chromosome location plot.
  Requires a sample-to-group map with exactly 2 groups.
type: standard
inputs:
  - name: expression_matrix.csv
    format: csv
    semantic_type: expression_matrix
    description: Gene expression matrix (genes x samples) with gene IDs in first column
  - name: sample_group_map.csv
    format: csv
    semantic_type: sample_metadata
    description: Two-column CSV mapping sample IDs (col 1) to group labels (col 2). Exactly 2 groups required.
outputs:
  - name: Diffanalysis.csv
    format: csv
    semantic_type: differential_expression_results
    columns: [gene_id, logFC, Pvalue, Padj, stat]
  - name: DEGs.csv
    format: csv
    semantic_type: filtered_degs
    columns: [gene_id, logFC, Pvalue, Padj, group]
  - name: Volcano.pdf
    format: pdf
    semantic_type: volcano_plot
  - name: Heatmap.pdf
    format: pdf
    semantic_type: heatmap_plot
  - name: Venn.pdf
    format: pdf
    semantic_type: venn_diagram
    conditional: true
  - name: Chromosome_location.pdf
    format: pdf
    semantic_type: chromosome_location_plot
    conditional: true
entry: scripts/main.R
parameters:
  - name: subcommand
    type: choice
    required: true
    bind: config
    description: Operation to perform (run, validate-input, validate-output)
  - name: --mat
    type: file
    required: true
    bind: upstream
    description: >
      Expression matrix CSV (genes x samples) with gene identifiers in first column.
      When a directory is provided, resolve: prefer "merged_expression.csv",
      fallback to first .csv matching 'expression' or 'gene'.
  - name: --map
    type: file
    required: true
    bind: upstream
    description: >
      Sample-to-group mapping CSV (sample_id, group).
      When a directory is provided, resolve: prefer "sample_group_map.csv",
      fallback to first .csv matching 'metadata' or 'sample'.
  - name: --method
    type: choice
    required: false
    default: auto
    bind: config
    description: >
      DE analysis method (deseq2, limma, edgeR, t, wilcox).
      Default auto-detects: integer count data → deseq2, normalized/log-transformed → limma.
      Explicitly passing --method bypasses auto-detection.
  - name: --p-set
    type: choice
    required: false
    default: padj
    bind: static
    description: Which p-value to use for filtering (p or padj)
  - name: --pvalue
    type: float
    required: false
    default: 0.05
    range: [0.001, 0.25]
    bind: static
    description: P-value threshold for significance
  - name: --logfc-cutoff
    type: float
    required: false
    default: 1.0
    range: [0.0, 10.0]
    bind: static
    description: Absolute log2 fold-change cutoff
  - name: --cutoff
    type: int
    required: false
    default: 10
    bind: static
    description: Minimum DEG count after filtering (below this triggers B2_FEW_DEGS)
  - name: --norm
    type: choice
    required: false
    default: TMM
    bind: static
    description: edgeR normalization method (TMM, RLE, upperquartile, none)
  - name: --model
    type: choice
    required: false
    default: glmFit
    bind: static
    description: edgeR model fitting method (glmFit, glmQLFit)
  - name: --top
    type: int
    required: false
    default: 20
    bind: static
    description: Top N genes for heatmap display and volcano labels
  - name: --force-imbalanced
    type: bool
    required: false
    default: false
    bind: static
    description: Override 10:1 sample proportion check
  - name: --rgs
    type: file
    required: false
    bind: upstream
    description: Related gene set CSV for Venn diagram intersection
  - name: --locate
    type: file
    required: false
    bind: upstream
    description: Chromosome annotation CSV (Gene, Chromosome, Start, End) for location plot
  - name: --tax-id
    type: string
    required: false
    default: "9606"
    bind: config
    description: NCBI taxonomy ID (9606=human, 10090=mouse)
  - name: --pheno-abbr
    type: string
    required: false
    bind: config
    description: Phenotype abbreviation for Venn diagram label
  - name: --gene
    type: string
    required: false
    bind: static
    description: Comma-separated gene IDs to label on volcano plot
  - name: --color-heat
    type: string
    required: false
    default: blue,white,red
    bind: static
    description: Heatmap color palette (comma-separated)
  - name: --color-panel
    type: string
    required: false
    bind: static
    description: Comma-separated colors for Venn diagram sets
  - name: --outdir
    type: file_out
    required: false
    default: .
    bind: framework
    description: Output directory for all result files
exceptions:
  - exit_code: 1
    code: B1_PROPORTION
    pattern: "sample proportion.*exceeds"
    nature: data_insufficient
    action: skip_with_warning
  - exit_code: 1
    code: B2_FEW_DEGS
    pattern: "gene count.*below.*cutoff"
    nature: data_insufficient
    action: skip_with_warning
  - exit_code: 1
    code: B3_MISSING_INPUT
    pattern: "input file not found"
    nature: data_corrupt
    action: halt
  - exit_code: 1
    code: B4_INVALID_COLUMNS
    pattern: "missing required column"
    nature: data_corrupt
    action: halt
  - exit_code: 1
    code: B5_METHOD_MISMATCH
    pattern: "not support.*matrix"
    nature: data_mismatch
    action: halt
  - exit_code: 1
    code: B6_NO_DEGS
    pattern: "no differentially expressed"
    nature: data_insufficient
    action: skip_with_warning
  - exit_code: 1
    code: B7_RGS_INTERSECTION
    pattern: "intersection.*below.*cutoff"
    nature: data_insufficient
    action: skip_with_warning
  - exit_code: 1
    code: B8_EMPTY_MATRIX
    pattern: "empty|no rows|no columns"
    nature: data_insufficient
    action: halt
  - exit_code: 1
    code: B9_SAMPLE_MISMATCH
    pattern: "sample.*mismatch|columns.*not.*match"
    nature: data_corrupt
    action: halt
  - exit_code: 1
    code: W001_DISK_FULL
    pattern: "disk full|no space"
    nature: resource
    action: halt
  - exit_code: 1
    code: W002_PERM_DENIED
    pattern: "permission denied"
    nature: resource
    action: halt
  - exit_code: 3
    code: E801_ENV_PKG
    pattern: "Missing required packages"
    nature: env_bug
    action: halt
  - exit_code: 1
    code: E802_UNSUPPORTED_TAXID
    pattern: "Unsupported tax_id"
    nature: data_mismatch
    action: halt
hardware:
  memory_gb: 4
  cpu: 2
  gpu: false
  runtime: "~2-10 minutes depending on method and matrix size"
---

# differential-analysis

## Node Function

Performs differential expression analysis comparing exactly two groups (case vs control) from a gene expression matrix. The pipeline:

1. Validates input data and sample proportion (max 10:1 ratio, overridable)
2. Runs the selected DE method (DESeq2, limma, edgeR, t-test, or Wilcoxon)
3. Filters results by log2 fold-change and p-value cutoffs
4. Generates volcano plot and heatmap
5. Optionally generates Venn diagram (if related gene set provided) and chromosome location plot (if annotation provided)
6. Writes full results (Diffanalysis.csv), filtered DEGs (DEGs.csv), and provenance (.run_result.json)

Count-based methods (DESeq2, edgeR) expect integer count data. limma and statistical tests work with normalized expression values.

## Expected Input

- **Expression matrix** (`--mat`): CSV file with gene identifiers in the first column and sample expression values in subsequent columns. Genes x samples.
- **Sample group map** (`--map`): Two-column CSV mapping sample IDs (first column) to group labels (second column). Must contain exactly 2 unique groups.

Optional inputs:
- **Related gene set** (`--rgs`): Single-column CSV of gene identifiers for Venn diagram intersection with DEGs.
- **Chromosome annotation** (`--locate`): CSV with columns Gene, Chromosome, Start, End for chromosome location plotting.

## Output Files

| File | Description |
|------|-------------|
| `Diffanalysis.csv` | Full DE results: gene_id, logFC, Pvalue, Padj [, stat] |
| `DEGs.csv` | Filtered DEGs with group column (Up/Down/Not) |
| `Volcano.pdf` | Volcano plot with significance thresholds |
| `Heatmap.pdf` | Heatmap of top DEGs, row-scaled |
| `Venn.pdf` | Venn diagram of DEGs vs related gene set (conditional) |
| `Chromosome_location.pdf` | Chromosome location plot via RCircos (conditional) |
| `.run_result.json` | Provenance: parameters, status, exceptions, file list, timestamps |

## Exceptions

**Data issues (B-codes, exit 1):**
- B1_PROPORTION: Sample ratio exceeds 10:1. Override with `--force-imbalanced`.
- B2_FEW_DEGS: DEG count below `--cutoff` threshold.
- B3_MISSING_INPUT: Required input file not found.
- B4_INVALID_COLUMNS: Missing required columns in input files.
- B5_METHOD_MISMATCH: Selected method incompatible with data.
- B6_NO_DEGS: No DEGs found at given cutoffs.
- B7_RGS_INTERSECTION: Intersection with related gene set below cutoff.
- B8_EMPTY_MATRIX: Expression matrix has zero rows or columns.
- B9_SAMPLE_MISMATCH: Sample IDs don't match between matrix and map.

**Write issues (W-codes, exit 1):**
- W001_DISK_FULL: No space left on output device.
- W002_PERM_DENIED: Cannot write to output directory.

**Environment issues (E-codes):**
- E801_ENV_PKG (exit 3): Missing required R/Bioconductor packages.
- E802_UNSUPPORTED_TAXID (exit 1): Unsupported NCBI taxonomy ID for RCircos.

## Usage Examples

```bash
# Run DESeq2 analysis with defaults
Rscript scripts/main.R run \
  --mat counts.csv \
  --map groups.csv \
  --outdir ./results

# Run limma with custom cutoffs and Venn diagram
Rscript scripts/main.R run \
  --mat expr_norm.csv \
  --map groups.csv \
  --method limma \
  --logfc-cutoff 1.5 \
  --pvalue 0.01 \
  --rgs phenotype_genes.csv \
  --pheno-abbr "AS" \
  --color-panel "#E64B35,#4DBBD5" \
  --outdir ./results

# Validate inputs before running
Rscript scripts/main.R validate-input \
  --mat counts.csv \
  --map groups.csv

# Validate outputs after running
Rscript scripts/main.R validate-output \
  --outdir ./results
```
