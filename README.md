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

## VCF-Polar: ancestral-allele annotation

[VCF-Polar 2](02.deleterious_mutation/README.md) provides a streaming command-line
tool for VCF/BCF ancestral-allele annotation, derived allele frequency and
per-sample derived dosage. It accepts existing INFO/AA, outgroup samples or a
reference-coordinate ancestral FASTA. The default preserves REF, ALT, GT and
phasing, and retains unresolved sites with decision reasons.

Install from the repository root and try the bundled example:

```sh
python -m pip install .
vcf-polar 02.deleterious_mutation/examples/tiny.vcf \
  --anc-list 02.deleterious_mutation/examples/outgroups.txt \
  -o tiny.annotated.vcf.gz --index --report tiny.summary.json
```

The original `02.deleterious_mutation/VCF-Polar.py` entry point remains available.
Version 2 changes the defaults and the outgroup consensus rule; see the
[migration notes](02.deleterious_mutation/README.md#migration-from-v11) and the
linked immutable v1.1 source when reproducing the original analysis.
