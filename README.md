# Differential Gene Expression in Lung Squamous Cell Carcinoma (TCGA RNA-seq)

A paired tumour-versus-normal differential expression analysis of lung squamous cell carcinoma (LUSC). It uses TCGA RNA-seq data (FPKM) from 51 patients, each contributing both a tumour sample and a matched adjacent normal sample, followed by pathway enrichment with GSEA.

> **Course project:** Bioinformatics (SBE3031), Systems & Biomedical Engineering, Cairo University (May 2026), team of 3. This is a copy of the team repository [FaridaaHesham/LUSC-DEGs-Analysis](https://github.com/FaridaaHesham/LUSC-DEGs-Analysis), kept with its commit history.
>
> **My contribution (Nour Bahgat):**
> - Fold-change analysis and dual-threshold DEG classification (`fold_change_analysis.py`).
> - The volcano plot (`volcano_plot.py`).
> - The up/down comparison bar chart (`bar_plot.py`).
> - The z-scored heatmap of the top 30 DEGs (`heatmap_plot.py`).

![Volcano plot](outputs/figures/volcano_plot.png)

## Pipeline

| Step | Script | Method |
|---|---|---|
| 1. Load and check | `statistical_analysis.py` | 19,648 genes × 51 paired samples; check for missing values (none found) |
| 2. Transform | 〃 | log₂(FPKM + 1) to reduce right skew |
| 3. Filter | 〃 | Keep genes with mean log₂ expression ≥ 2 → **16,545 genes tested** (3,103 removed) |
| 4. Normality check | 〃 | Shapiro–Wilk on paired differences of 10,000 random genes: 82.5 % pass, which supports a parametric test |
| 5. Hypothesis test | 〃 | Two-sided paired t-test per gene (`scipy.stats.ttest_rel`) |
| 6. Multiple testing | 〃 | Benjamini–Hochberg FDR correction → 12,916 genes with FDR < 0.05 |
| 7. Effect size | `fold_change_analysis.py` | log₂FC = mean log₂(tumour) − mean log₂(normal); a DEG needs FDR < 0.05 **and** \|log₂FC\| > 1 |
| 8. Visualise | `volcano_plot.py`, `bar_plot.py`, `heatmap_plot.py` | Volcano plot; p-value-only vs dual-threshold counts; heatmap of the top 30 DEGs by \|log₂FC\| |
| 9. Enrichment | `generate_gsea_gct_cls.py` | Exports GCT/CLS files for GSEA Desktop (MSigDB Hallmark, signal-to-noise ranking, 1,000 permutations, FDR q < 0.25) |

## Results

- **5,437 DEGs: 1,850 upregulated and 3,587 downregulated** in tumour tissue.
- The top-30 heatmap separates tumour from normal samples cleanly across all 51 patients:
  - **downregulated in tumour:** alveolar and airway epithelial markers (SFTPC, SFTPA1/2, AGER, CLDN18, CA4, SCGB1A1, CYP4B1);
  - **upregulated in tumour:** squamous differentiation and carcinoma markers (KRT5, KRT6B, KRT14, KRT16, DSG3, SERPINB5, AKR1B10, GJB6).

| Up/down DEG counts | Top-30 DEG heatmap |
|---|---|
| ![Bar plot](outputs/figures/bar_plot.png) | ![Heatmap](outputs/figures/heatmap_plot.png) |

## Running

```bash
pip install -r requirements.txt
python statistical_analysis.py     # writes outputs/all_genes_results.csv and DEGs_significant.csv
python fold_change_analysis.py     # writes DEG classes and the figures
python generate_gsea_gct_cls.py    # writes the GSEA input files
```

`lusc-rsem-fpkm-tcga-t_paired.txt` (tumour) and `lusc-rsem-fpkm-tcga_paired.txt` (matched normal) are the TCGA LUSC RSEM FPKM matrices provided for the course.

## Team
Farida Abdelrahman (statistical testing and BH correction) · Nour Bahgat (fold change, DEG classification, visualisation) · Nour Allah Hassan (GSEA)
