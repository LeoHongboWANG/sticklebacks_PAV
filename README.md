# Stickleback gene-loss analyses

Selected scripts accompanying *Demography and selection shape gene loss in freshwater sticklebacks — A pangenomic perspective*.

| Script | Analysis |
|---|---|
| `scripts/ancestral_genes.py` | Infer ancestral presence using pooled outgroup presence >90%; optionally count individual gene losses. |
| `scripts/repeatability.py` | Identify recurrent freshwater-biased absences and permute ecotype labels. |
| `scripts/phylogenetic_model.R` | Fit the phylogenetic Poisson model of gene loss, NePi and ecotype. |

The Python scripts implement the documented analysis rules and were reconstructed for this repository. The R script adapts the original model script for use with local input paths. The historical permutation seed is unknown. This repository provides selected scripts; input datasets and the full workflow are outside its scope.

## Dependencies and usage

Python: `pip install -r requirements.txt`. R: `MCMCglmm`, `ape`, `coda` (`parallel` is included with R).

```sh
python scripts/ancestral_genes.py --outgroup-pav ANC_genePAV.txt --outdir results/ancestral
python scripts/repeatability.py --counts gene_loss_by_population_and_gene.csv --outdir results/repeatability --permutations 2000 --seed 20261026
Rscript scripts/phylogenetic_model.R DATA_ROOT results/phylogenetic
```

- Ancestral input is a tab-separated binary matrix with GeneID in the first column and outgroup samples in the remaining columns. To count losses by individual, supply `--pav FULL_PAV.tsv --samples SAMPLE_TABLE.tsv`, with `SampleID` and `Population` in the sample table.
- Repeatability input: long CSV with `GeneID`, `Population`, `EcotypeSimple` (`Marine`/`Freshwater`), `PopulationSize`, `AbsentCount`. By default, SWE-LUN is excluded. A gene qualifies when at least 50% of individuals lack it in at least half of the freshwater populations and no more than one marine population. Use `--keep-genes IDs.txt` to restrict the gene set.
- Model input: `DATA_ROOT/p_dis_2.nwk` and `DATA_ROOT/analysis_outputs/excov_depth_robustness/sample_losttxt_with_coverage.csv`; the CSV requires `SampleID`, `EcotypeSimple`, `NePi`, `LostTxt`. Four chains run in parallel. Add `--validate-only` to check input preparation without fitting the model.

The repository also contains `02.deleterious_mutation/VCF-Polar.py`.
