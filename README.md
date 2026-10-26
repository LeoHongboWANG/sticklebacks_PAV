# Stickleback PAV revision analyses

This directory contains eight scripts and the available inputs for a bounded reproduction of the gene-loss revision analyses. It is an addition to the existing `sticklebacks_PAV` repository, which already contains `02.deleterious_mutation/VCF-Polar.py`.

**This is not a complete reproduction of every analysis in the manuscript.** See `docs/MISSING_ANALYSES.md` before describing the repository in Code availability. No code has been pushed by this preparation step.

Three Python scripts and one base-R check were reconstructed during revision preparation. Four other scripts were recovered from the local project and made portable. Their origins and source hashes are recorded in `docs/source_manifest.json`. Reconstructed scripts are not presented as the historical scripts used for the published figures.

## Requirements

Python 3.10 or later and NumPy (`python -m pip install -r requirements.txt`). The three reconstructed Python scripts use only NumPy and the Python standard library. The two recovered coverage scripts use the standard library and invoke R when available.

The R scripts require R and the packages `MASS`, `nlme`, `MCMCglmm`, `ape`, and `coda`. The phylogenetic script uses `parallel::mclapply`, which requires a Unix-like system for multiple workers. It retains four chains, seeds 1001–1004, 1,300,000 iterations, 300,000 burn-in iterations and thinning by 500. Do not start the full model just to check installation. The previous run took approximately 48 hours.

## Included scripts

| Script | Purpose | Origin and verification |
|---|---|---|
| `scripts/group_comparisons.R` | Individual-level unadjusted Welch tests for Fig. 2a/b | Reconstructed; text t, df and P values reproduced, identifying stale figure labels |
| `scripts/ancestral_genes.py` | Strict pooled outgroup presence >90%; optional sample loss counts | Reconstructed; 26,181 ancestral genes reproduced; all 809 sample loss counts match saved counts when the full local matrix is supplied |
| `scripts/repeatability.py` | 39-population screen and new seeded ecotype-label permutations | Reconstructed; observed 86 reproduced; a new 2,000-permutation run completed |
| `scripts/compare_candidates.py` | Fixed Table 10 list versus repeatability set and ancestral ascertainment; audit saved permutation output | Reconstructed; overlap 20; ancestral membership 8/86 and 0/44; saved corrected P reproduced |
| `scripts/phylogenetic_model.R` | Main NePi/ecotype phylogenetic Poisson model | Recovered original; portable paths and validation mode added; inputs matched; chronos issued optimization warnings; full MCMC not rerun in this preparation |
| `scripts/priority1_followups.R` | Alternative Ne estimators, count models, rarefaction, fixed/segregating partitions and locality models | Recovered original; complete script rerun on bundled source data |
| `scripts/analyze_excov_depth_robustness.py` | Coverage summaries and population depth regressions from raw `.excov` files | Recovered original; path defaults changed; raw coverage analysis not rerun here |
| `scripts/analyze_top_recurrent_exon_support.py` | Direct coverage summaries at the retained descriptive candidates | Recovered original; path defaults changed; does not reconstruct how the original 44 genes were selected |

## Run the small reproducible analyses

Run from this directory. These commands write to `results/`, which is ignored by Git.

```bash
python scripts/ancestral_genes.py \
  --outgroup-pav data/ANC_genePAV.txt.gz \
  --pav data/Ppu.genePAV.txt.gz --samples data/lost.txt \
  --outdir results/ancestral

python scripts/repeatability.py \
  --counts data/analysis_outputs/gene_loss_by_population_and_gene.csv.gz \
  --permutations 2000 --seed 20261026 \
  --outdir results/repeatability

python scripts/compare_candidates.py \
  --ancestral results/ancestral/ancestral_presence.csv \
  --screen results/repeatability/gene_screen.csv \
  --table44 data/analysis_outputs/top44_recurrent_loss_table.csv \
  --saved-null data/analysis_outputs/repeatability_permutation_distribution_noSWELUN.csv \
  --outdir results/candidate_comparison

Rscript scripts/group_comparisons.R data results/group_comparisons
Rscript scripts/priority1_followups.R data results/priority1_followups

# Validate model inputs and preprocessing without fitting the four MCMC chains.
Rscript scripts/phylogenetic_model.R data results/phylogenetic_model --validate-only

# Optional full computation, requiring substantial time and memory.
# Rscript scripts/phylogenetic_model.R data results/phylogenetic_model
```

The ancestral command also recomputes the 809 individual loss counts from the bundled full local binary PAV matrix. That matrix has 975 sample columns, whereas the manuscript growth analysis specifies 967. It is included as a source snapshot, not silently relabeled as the final manuscript matrix. The exact growth-analysis sample list still needs to be recovered.

To rebuild the coverage results, supply the raw coverage directory explicitly:

```bash
python scripts/analyze_excov_depth_robustness.py PATH_TO_EXCOV_DIRECTORY \
  --outdir results/coverage
python scripts/analyze_top_recurrent_exon_support.py \
  --excov-dir PATH_TO_EXCOV_DIRECTORY --outdir results/candidate_coverage
```

The raw `.excov` directory is not bundled. The included small summary tables permit the recovered model scripts to run independently of that directory.

## Data and statistical definitions

- `ANC_genePAV.txt.gz`: tab-separated GeneID column followed by 78 binary outgroup sample columns. One means present; zero means absent. Ascertainment uses pooled individuals, not a majority vote across five populations.
- `gene_loss_by_population_and_gene.csv.gz`: long-format saved counts with `GeneID`, `Population`, `EcotypeSimple`, `PopulationSize`, `PresentCount`, `AbsentCount`, and `AbsentRate`. A population qualifies when `2 * AbsentCount >= PopulationSize`; the rounded rate is not used for thresholding.
- The repeatability screen excludes SWE-LUN, leaving 27 freshwater and 12 marine populations. A gene qualifies in at least 14 freshwater populations and at most one marine population. This is a set-level recurrence test, not 86 separate selection tests.
- The local count table contains 38,601 genes. The manuscript reports a final 36,174-gene matrix. The latter's complete GeneID list is not recovered. The `--keep-genes` option allows an explicit final list to be applied once recovered. No list is silently inferred from the manuscript totals.
- The original 44-gene list is retained as a fixed descriptive list. It is not recreated by selecting the top 44 rows. The local file called `top44_proxy` must not replace it.
- The saved historical permutation distribution has 2,000 values, mean 1.6115, and one value at least as large as 86. The corrected P is `(1 + 1) / (2000 + 1) = 0.00099950025`.
- The historical random seed is unknown. The new reconstruction uses seed 20261026 and produces mean 1.5365 and the same rounded P of 0.001. Its output is explicitly named `new_permutation_distribution.csv`; it does not replace the original figure's saved distribution.
- The phylogenetic model scales `log10(NePi + 1)` and uses Marine as reference. The population regressions generally use unscaled log10 Ne and Freshwater as the factor reference. Coefficients from these models must not be compared as if their scales and reference levels were identical.

See `docs/VERIFICATION.md` for the verification record, including tree-preprocessing warnings.
