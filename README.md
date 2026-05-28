# IBD-seq Unified Meta-Regression

Application of the [phenome-wide unified meta-regression model](https://github.com/rivas-lab/phenome-wide-unified-model) to rare-variant summary statistics from the IBD-seq consortium (Crohn's disease, ulcerative colitis, and combined IBD; European ancestry).

The unified model integrates single-variant association statistics with functional priors — genomic **constraint** ([wgs-constraint-llm](https://github.com/rivas-lab/wgs-constraint-llm)) and **AlphaMissense** ([google-deepmind/alphamissense](https://github.com/google-deepmind/alphamissense)) — to estimate gene-level effect-size relationships via a weighted least-squares meta-regression. Background and method: [Osthoag et al., bioRxiv 2025](https://www.biorxiv.org/content/10.1101/2025.01.23.634522v1).

## Inputs

| Resource | Source | Notes |
|---|---|---|
| IBD-seq META summary stats | IBD-seq consortium | METAL format (`MarkerName`, `Allele1`, `Allele2`, `Freq1`, `Effect`, `StdErr`); per phenotype: `cd_EUR`, `uc_EUR`, `ibd_EUR` |
| Constraint probabilities (`prob_0`) | [rivas-lab/wgs-constraint-llm](https://github.com/rivas-lab/wgs-constraint-llm) | HMM RGC ≥0.9 over 20bp, WES-restricted |
| AlphaMissense (hg38) | [google-deepmind/alphamissense](https://github.com/google-deepmind/alphamissense) | `am_pathogenicity` per SNV |
| Variant annotations (gene, consequence) | Genebass / VEP | Used to attach `gene`, `pLoF`, `missense` labels |

## Pipeline

### 1. Allele alignment against the reference

The IBD-seq META files report effects on `Allele1`, which is **not** guaranteed to be the ALT allele. Before joining with constraint / AlphaMissense / annotations (all keyed on REF/ALT):

- For variants where `Allele1 == REF`: **flip the effect** — negate `Effect` (BETA), set `Freq1 → 1 − Freq1`, and swap `Allele1`/`Allele2`.
- For variants where `Allele1 == ALT`: keep as is.
- **Do not drop flip-needed variants** — flipping preserves them with the correct effect direction.
- Compute `AF` from the corrected `Freq1` (this is the MAF used downstream; the MAF ≤ 0.05 filter is applied in step 4).

### 2. Annotate with gene and variant consequence

Join each variant to its `gene` and `annotation` (`pLoF`, `LC`, `missense`, …) via the Genebass variant-level annotation table. The output of this step is a TSV with at minimum: `chr`, `pos`, `ref`, `alt`, `gene`, `annotation`, `AF`, `BETA`, `SE`.

### 3. Merge functional priors (`metareg_prep.py`)

Run [`metareg_prep.py`](https://github.com/rivas-lab/phenome-wide-unified-model/blob/main/metareg_prep.py) from the upstream repo. This adds:

- `prob_0` — constraint probability (merge on `chr`, `pos`)
- `am_pathogenicity` — AlphaMissense score (merge on `chr`, `pos`, `ref`, `alt`)
- `pLoF_indicator` — 1 if `annotation ∈ {pLoF, LC}`
- `missense_indicator` — 1 if `annotation == missense`

### 4. Unified meta-regression (`unified_reg_MAF.05.py`)

Run [`unified_reg_MAF.05.py`](https://github.com/rivas-lab/phenome-wide-unified-model/blob/main/unified_reg_MAF.05.py). For each gene this fits a weighted least-squares model:

```
BETA  ~  log_constraint + log_pathogenicity + pLoF_indicator + missense_indicator
weights = 1 / SE^2
```

restricted to variants with `AF ≤ 0.05`. The output is one row per gene with the four coefficient estimates and p-values, plus the overall model p-value.

Repeat steps 1–4 independently for `cd_EUR`, `uc_EUR`, and `ibd_EUR`.

## Dashboard (`app.py`)

A Streamlit app for interactive exploration of the gene-level results. Place one or more `*.unified_results.tsv.gz` files in the `data/` directory and run:

```bash
streamlit run app.py
```

### Sidebar controls

| Control | Description |
|---|---|
| **Dataset** | Switch between `cd_EUR`, `uc_EUR`, `ibd_EUR` result files |
| **p-value threshold** | Significance cutoff (default 0.05); highlighted in all views |
| **Predictors to show** | Toggle log constraint, log pathogenicity, pLoF indicator, missense indicator |
| **Hide genes without p_model** | Filter genes with no overall model p-value |
| **Clear cache & reload** | Force re-read of data files from disk |

### Tab 1 — Table

Sortable gene results table with coefficient and p-value columns for each active predictor. Significant cells (p < threshold) are highlighted green; nominally associated cells (p < 0.25) in yellow. A CSV download button is provided.

### Tab 2 — Rankings

Horizontal bar chart of genes ranked by overall model p-value (−log₁₀ scale). Significant genes are shown in green; non-significant in grey. A dashed red line marks the selected significance threshold.

### Tab 3 — Heatmap

Genes × predictors heatmap, sorted by model p-value (most significant at top). Two display modes:

- **−log₁₀(p-value)** — `Reds` colorscale; color intensity reflects association strength per predictor. The color range is capped at the 99th percentile of finite values so that extremely significant genes (p = 0 → −log₁₀ = ∞) do not compress the entire scale to near-zero.
- **Coefficient value** — `RdBu_r` diverging colorscale (red = positive, blue = negative); extreme outliers clipped at the 95th percentile of absolute values for readability.

### Tab 4 — Gene Detail

Per-gene view with:

- Overall model p-value and −log₁₀(p_model) metrics
- Bar chart of predictor coefficients (full opacity = significant, faded = not significant)
- Detailed table of coefficient estimates, p-values, and significance flags

## Reproduce

```text
1. Clone the upstream pipeline:
     git clone https://github.com/rivas-lab/phenome-wide-unified-model
2. Obtain constraint and AlphaMissense files (links above).
3. Obtain IBD-seq META summary stats (cd/uc/ibd_EUR.meta.all.1.txt.gz).
4. Apply allele-alignment fix-up (step 1 above) → produces *.processed.tsv.gz.
5. Attach gene + consequence annotations (step 2).
6. python metareg_prep.py  <annotated.tsv>  <prepped.tsv>
7. python unified_reg_MAF.05.py  <prepped.tsv>  <gene_results.tsv.gz>
```

## References

- Unified meta-regression model — [rivas-lab/phenome-wide-unified-model](https://github.com/rivas-lab/phenome-wide-unified-model)
- Constraint priors — [rivas-lab/wgs-constraint-llm](https://github.com/rivas-lab/wgs-constraint-llm)
- AlphaMissense — [google-deepmind/alphamissense](https://github.com/google-deepmind/alphamissense)
- Genebass — [https://app.genebass.org/](https://app.genebass.org/)
- Preprint — Osthoag et al., *A unified meta-regression model for rare-variant association studies*, bioRxiv 2025 ([10.1101/2025.01.23.634522](https://www.biorxiv.org/content/10.1101/2025.01.23.634522v1))
