# Stickleback gene-loss analyses

Selected scripts accompanying *Demography and selection shape gene loss in freshwater sticklebacks — A pangenomic perspective*.

| Script | Analysis |
|---|---|
| `scripts/ancestral_genes.py` | Infer ancestral presence using pooled outgroup presence >90%; optionally count individual gene losses. |
| `scripts/repeatability.py` | Screen recurrent freshwater losses and permute ecotype labels. |
| `scripts/phylogenetic_model.R` | Fit the phylogenetic Poisson model of gene loss, NePi and ecotype. |

The two Python scripts were reconstructed from the documented analysis rules; the R script is a portable version of the original model script. The historical permutation seed is unknown. Input datasets and the full analysis workflow are not included here.

## Dependencies and usage

Python: `pip install -r requirements.txt`. R: `MCMCglmm`, `ape`, `coda` (`parallel` is included with R).

```sh
python scripts/ancestral_genes.py --outgroup-pav ANC_genePAV.txt --outdir results/ancestral
python scripts/repeatability.py --counts gene_loss_by_population_and_gene.csv --outdir results/repeatability --permutations 2000 --seed 20261026
Rscript scripts/phylogenetic_model.R DATA_ROOT results/phylogenetic
```

- Ancestral input: tab-separated binary matrix, first column GeneID, remaining columns outgroup samples. Optional individual counts require `--pav FULL_PAV.tsv --samples SAMPLE_TABLE.tsv`, with `SampleID` and `Population` in the sample table.
- Repeatability input: long CSV with `GeneID`, `Population`, `EcotypeSimple` (`Marine`/`Freshwater`), `PopulationSize`, `AbsentCount`. Defaults exclude SWE-LUN, classify a population as absent at ≥50%, and require absence in ≥half of freshwater populations and ≤one marine population. Use `--keep-genes IDs.txt` to restrict the gene set.
- Model input: `DATA_ROOT/p_dis_2.nwk` and `DATA_ROOT/analysis_outputs/excov_depth_robustness/sample_losttxt_with_coverage.csv`; the CSV requires `SampleID`, `EcotypeSimple`, `NePi`, `LostTxt`. Four chains run in parallel. Add `--validate-only` to check input preparation without fitting the model.

Retain the existing `02.deleterious_mutation/VCF-Polar.py` in the repository.
