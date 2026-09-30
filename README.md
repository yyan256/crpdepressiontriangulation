# Dietary inflammatory potential, C-reactive protein, and depression

Analysis code for the manuscript:

> **Dietary inflammatory potential, C-reactive protein, and depression: triangulation of NHANES observational mediation and Mendelian randomization**
> *Psychiatry Research*, manuscript PSY-D-26-02386

This repository contains the code used to reproduce all analyses reported in the revised manuscript. The study combines (i) cross-sectional mediation analyses of the Dietary Inflammatory Index (DII), high-sensitivity C-reactive protein (hs-CRP), and depressive symptoms in NHANES, and (ii) two-sample Mendelian randomization (MR) analyses of CRP, MDD, and related exposures.

---

## Repository structure

```
.
├── env/                      # Environment setup
│   ├── env_setup.sh          # One-shot install (system deps + R + Python)
│   └── setup_R_packages.R    # Robust R-package installer (multi-mirror)
├── data_download/            # Reference panel & GWAS summary statistics
│   ├── download_ref_panel.sh # plink2 binary + 1000 Genomes EUR panel
│   ├── download_gwas_vcf.sh  # Download GWAS VCFs from OpenGWAS/MRCIEU
│   └── download_gwas_vcf.R   # R variant of the VCF downloader
├── nhanes_mediation/         # NHANES mediation (observational arm)
│   ├── nhanes_prep.R                # Download + build the 5-cycle analytic sample
│   ├── fix_dii_rebuild.R            # Recompute DII + fix smoker/education/race coding → v2 sample
│   ├── mediation_main.R             # MAIN: lavaan WLSMV (probit) + glmnet-ridge sensitivity
│   ├── mediation_lavaan_wlsmv.R     # lavaan WLSMV (probit) standalone reproducibility check
│   ├── mediation_delta_ci.R         # Delta-method 95% CIs (sandwich vcov; logit sensitivity)
│   └── mediation_sensitivity_continuous.R  # Continuous / total-PHQ-9 / PHQ-9≥10
└── mr_analysis/              # Two-sample MR (genetic arm)
    ├── mr_main.R                    # Forward CRP→MDD and reverse MDD→CRP
    ├── mr_il6r_drug_target.R        # IL6R cis drug-target MR
    ├── mr_two_step_mediation.R      # Two-step MR (nutrient → CRP → MDD)
    └── mr_report_table.R            # Full STROBE-MR reporting battery (R1-9)
instruments/                 # Complete instrument lists (versioned; see below)
    ├── instruments_CRP_to_MDD.csv
    ├── instruments_MDD_to_CRP.csv
    ├── instruments_VitD_to_CRP.csv
    ├── instruments_VitD_to_MDD.csv
    └── instruments_IL6R_to_MDD.csv
```

---

## Data sources

All data are publicly available and are **not** stored in this repository.

| Data | Source | Identifier |
|---|---|---|
| NHANES (2005–2018) | [NHANES](https://wwwn.cdc.gov/nchs/nhanes/) via the `nhanesA` R package | 5 cycles with hs-CRP: 2005-2006, 2007-2008, 2009-2010, 2015-2016, 2017-2018 |
| CRP GWAS | [OpenGWAS](https://gwas.mrcieu.ac.uk/) | `ieu-b-35` (Ligthart et al. 2018) |
| MDD GWAS | OpenGWAS | `ieu-b-102` (PGC MDD, Howard et al. 2019) |
| Vitamin D GWAS | OpenGWAS | `ieu-b-4812` |
| Fibre intake GWAS | OpenGWAS | `ukb-b-19085` |
| LD reference | [MRCIEU](https://mrcieu.mrc.ac.uk/) | 1000 Genomes Phase 3, EUR |

Analytic sample (observational arm): **N = 24,019** across the five hs-CRP cycles.

---

## Software requirements

- **R ≥ 4.3** with packages: `nhanesA`, `survey`, `glmnet`, `lavaan`, `MendelianRandomization`, `data.table`, `dplyr`, `tidyr` (installed by `env/setup_R_packages.R`)
- **plink2** (2.0) for LD clumping
- **tabix** for VCF access
- Linux environment (the scripts assume paths under `/root/autodl-tmp/`; adjust `OUT`/`VCFDIR`/`REF` variables at the top of each script if running elsewhere)

---

## Reproduction workflow

**Step 0 — environment**

```bash
bash env/env_setup.sh
```

**Step 1 — reference panel & GWAS data**

```bash
bash data_download/download_ref_panel.sh   # plink2 + 1000G EUR panel
bash data_download/download_gwas_vcf.sh    # CRP, MDD, VitD, fibre GWAS VCFs
```

**Step 2 — NHANES mediation (observational arm)**

```bash
Rscript nhanes_mediation/nhanes_prep.R                 # build analytic sample
Rscript nhanes_mediation/fix_dii_rebuild.R             # correct DII → v2 sample
Rscript nhanes_mediation/mediation_main.R              # main results
Rscript nhanes_mediation/mediation_lavaan_wlsmv.R      # lavaan reproducibility
Rscript nhanes_mediation/mediation_delta_ci.R          # 95% CIs
Rscript nhanes_mediation/mediation_sensitivity_continuous.R  # sensitivity
```

**Step 3 — two-sample MR (genetic arm)**

```bash
Rscript mr_analysis/mr_main.R                 # forward + reverse MR
Rscript mr_analysis/mr_il6r_drug_target.R     # IL6R drug-target MR
Rscript mr_analysis/mr_two_step_mediation.R   # two-step MR
Rscript mr_analysis/mr_report_table.R         # full STROBE-MR reporting table
```

---

## Key results reproduced by this code

**Observational mediation (DII → hs-CRP → depressive symptom dimensions)**

| Model | Indirect proportion |
|---|---|
| Somatic, without BMI | 14.1% |
| Somatic, with BMI | 6.5% |
| Cognitive, without BMI | 11.6% |
| Cognitive, with BMI | 6.2% |

**Two-sample MR**

| Analysis | Instruments (clump before → after) | Mean F | MR-Egger intercept (p) | Cochran Q (p) | Steiger direction (p) |
|---|---|---|---|---|---|
| CRP → MDD | 3,951 → 54 | 180.4 | 0.0008 (0.538) | 86.8 (0.0024) | correct (0.666) |
| MDD → CRP | 4,613 → 23 | 41.5 | −0.0078 (0.179) | 42.0 (0.0062) | correct (0.021) |
| IL6R → MDD | 117 → 2 | — (2 IVs) | — | — | — |

---

## Complete instrument lists

The full, per-analysis instrument tables (harmonized effect/other alleles, exposure and outcome beta/SE/p, and per-SNP F-statistics) are provided as versioned CSV files in `instruments/`:

| File | Analysis | Harmonized instruments |
|---|---|---|
| `instruments_CRP_to_MDD.csv` | CRP → MDD (forward) | 54 |
| `instruments_MDD_to_CRP.csv` | MDD → CRP (reverse) | 23 |
| `instruments_VitD_to_CRP.csv` | Vitamin D → CRP (two-step, step 1) | 42 |
| `instruments_VitD_to_MDD.csv` | Vitamin D → MDD (total effect) | 94 |
| `instruments_IL6R_to_MDD.csv` | IL6R cis → MDD (drug-target) | 2 |

The fibre-intake GWAS (`ukb-b-19085`) contained no genome-wide significant (p < 5 × 10⁻⁸) variants, so no instrument list is provided for the fibre arm (consistent with the null fibre step reported in the manuscript). All instruments were LD-clumped at r² < 0.001 within 10 Mb using the 1000 Genomes EUR panel.

---

## Methodological notes (added in revision)

**Mediation scale.** The primary mediation results are estimated with `lavaan` **WLSMV** (probit link for the binary outcome), so that the path-a coefficient (linear, on log-CRP) and path-b coefficient (probit) are jointly estimated on a single, well-defined scale, and the indirect effect `a×b` and its proportion `a×b/(a×b+c′)` carry delta-method confidence intervals directly from `lavaan`. The earlier "linear-coefficient × logit-coefficient" mixture is avoided. The `glmnet` ridge-logistic model is retained only as a **sensitivity** check (logit scale; `intercept = FALSE` to avoid double intercepts).

**Variable coding (fixed in this revision).** The analytic sample defines *current smoker* from SMQ020 (≥100 cigarettes in lifetime) **and** SMQ040 (now smokes every day/some days): never-smoked (SMQ020 = 2) → 0; current (SMQ040 ∈ {1,2}) → 1; former (SMQ040 = 3) → 0. Education uses the correct DMDEDUC2 mapping (<High school = 1,2; High school graduate = 3; >High school = 4,5). Race is collapsed to 4 categories (Hispanic = Mexican American + Other Hispanic). Vitamin E in the DII uses `DR1TATOC` (not `DR1TVE`).

**IL6R Wald-ratio p-values.** In `instruments_IL6R_to_MDD.csv`, the `p_outcome` column is the raw MDD GWAS SNP-outcome p-value; the Wald-ratio p-values reported in the manuscript (Section 3.9) are instead computed as `2·Φ(−|β/SE|)` from the Wald ratio. Both are reproducible from the CSV columns (`beta_outcome`, `se_outcome`, `beta_exposure`).

---

## License

Code is released under the MIT License. Data are governed by their respective original licences (NHANES and OpenGWAS).
