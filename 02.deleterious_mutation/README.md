# VCF-Polar 2

VCF-Polar annotates ancestral alleles, derived allele frequencies and per-sample
derived dosage in VCF/BCF files. The default mode keeps the original REF, ALT,
GT, phasing, AD and genotype likelihoods. It streams records in input order.

Version 2 supports **biallelic SNPs**. Indels, multiallelic sites and symbolic
variants are retained with an explicit unresolved status, unless filtering is
requested. Consensus support is a fraction of observed outgroup allele copies,
**not a calibrated probability of ancestry**.

## Install

Python 3.10 or newer and pysam are required. From the repository root:

```sh
python -m pip install .
vcf-polar --version
vcf-polar --help
```

For a new checkout:

```sh
git clone https://github.com/LeoHongboWANG/sticklebacks_PAV.git
cd sticklebacks_PAV
python -m pip install .
```

The original entry-point path also works after installing dependencies:

```sh
python 02.deleterious_mutation/VCF-Polar.py input.vcf -o annotated.vcf
```

This package is installed from this repository; it has not been published to
PyPI or Bioconda. The original v1.1 script acknowledged Ari's pipeline and
[Biostars post 173338](https://www.biostars.org/p/173338/). This repository does
not currently declare a software license; no new license is asserted by this
update.

## Try the bundled example

Run these commands from the repository root:

```sh
vcf-polar 02.deleterious_mutation/examples/tiny.vcf \
  --anc-list 02.deleterious_mutation/examples/outgroups.txt \
  -o tiny.annotated.vcf.gz --index \
  --report tiny.summary.json --sites-report tiny.decisions.tsv
```

There are two focal samples and two outgroups. The expected outgroup decisions
are:

| Position | Decision | Explanation |
|---|---|---|
| 1 | PASS, AA=A | Both outgroups are reference homozygotes |
| 2 | PASS, AA=G | Both outgroups are alternative homozygotes |
| 3 | LOW_SUPPORT | Both outgroups are heterozygous; maximum support is 0.5 |
| 4 | TOO_FEW_OUTGROUPS | One outgroup fails GQ and the other is missing |
| 5 | UNSUPPORTED_VARIANT | An indel, outside the v2 SNP inference scope |
| 6 | PASS, AA=G | Both outgroups support G |

REF, ALT and original GT values stay unchanged. The JSON report counts every
decision; the TSV contains every input site, including sites dropped by an
explicit filtering option. With the default annotate mode, all six sites are
written.

## Three evidence sources

### Existing INFO annotation

```sh
vcf-polar population.vcf.gz --source info --aa-tag AA \
  -o population.annotated.vcf.gz --report summary.json
```

The selected INFO tag must be declared in the input header. A value such as
`AA=A|additional-information` uses the first component. The original input
annotation is recorded in `VP_INPUT_AA`. Missing, ambiguous or unsupported bases
remain unresolved. Lowercase calls are rejected by default; use
`--allow-lowercase` only when appropriate for the source's confidence encoding.

### Outgroup individuals in the same VCF

```sh
vcf-polar population.with_outgroups.vcf.gz \
  --anc-list outgroups.txt --min-outgroups 2 \
  --min-dp 3 --min-gq 20 --consensus-threshold 0.8 \
  --max-outgroup-missing 0.5 \
  -o population.annotated.vcf.gz --index --report summary.json
```

Sample names must match exactly. Empty lists, unknown names, duplicates and
overlapping ingroup/outgroup lists are errors.

Every complete, quality-passing outgroup genotype contributes its allele
copies. A heterozygote contributes one copy to each allele. Haploid and
polyploid complete genotypes contribute their observed ploidy. Partial genotypes
are excluded from outgroup inference. The missing fraction is the fraction of
listed outgroup individuals failing any genotype, DP or GQ requirement.
Missing DP/GQ fails an enabled filter; `--min-dp 0 --min-gq 0` explicitly disables
those filters.

For example, one G/G plus three A/G genotypes gives G support 5/8, not 100%.
The default minimum is one passing outgroup; increase `--min-outgroups` when
replicate samples are available. This is pooled allele-copy consensus, without
species-tree weighting, substitution models or correction for reference bias,
recurrent mutation, introgression or incomplete lineage sorting. Choose
outgroups with the biological question and their divergence in mind.

An external variants-only outgroup VCF must not be interpreted as reference
homozygous wherever it has no record. Produce callable genotypes in the joint
coordinate system before using this mode.

### Ancestral FASTA in reference coordinates

```sh
samtools faidx ancestor.projected.fa
vcf-polar population.vcf.gz --ancestral-fasta ancestor.projected.fa \
  -o population.annotated.vcf.gz --report summary.json
```

The FASTA must already be aligned/projected into the **same assembly and
coordinates as the input VCF**, retaining missing reference-coordinate positions
as N. An ordinary outgroup genome is not a suitable input. Contig names must
match. Bases beyond available sequence or on absent contigs get
`FASTA_MISSING`. Lowercase and N calls are unresolved by default.

For the bundled FASTA, first run
`samtools faidx 02.deleterious_mutation/examples/ancestor.projected.fa`.

### Source selection and fallback

`--source auto` selects outgroups if `--anc-list` is present, otherwise FASTA if
supplied, otherwise INFO. If both outgroup names and FASTA are supplied, choose
`--source fasta` explicitly; the listed outgroups are then excluded from focal
frequency calculations, not used as a second inference source. Similarly,
`--source info --anc-list outgroups.txt` selects INFO evidence and excludes those
samples from focal frequencies.

There is **no implicit fallback** to old INFO/AA when outgroup or FASTA evidence
fails. Add `--fallback-info` to allow it explicitly. Such successful calls have
`VP_SOURCE=info_fallback` and record the primary failure in `VP_FALLBACK`.
A valid ancestral base outside REF/ALT is marked `AA_NOT_IN_ALLELES` and has no
derived counts or dosage in this release.

## Output annotations

These are custom tags except for the standard `AA` tag. Consumers must select
`VP_STATUS=PASS` when they need resolved polarization.

| Tag | Meaning |
|---|---|
| INFO/AA | Inferred ancestral base |
| INFO/VP_INPUT_AA | Original selected input ancestral annotation |
| INFO/VP_STATUS | PASS or a reason for leaving the site unresolved |
| INFO/VP_SOURCE | info, outgroup, fasta, or info_fallback |
| INFO/VP_N_OUT | Quality-passing outgroup individuals |
| INFO/VP_N_OUT_ALLELES | Outgroup allele copies contributing to consensus |
| INFO/VP_SUPPORT | Largest outgroup allele-copy fraction, not a posterior |
| INFO/VP_FALLBACK | Primary failure when INFO fallback was used |
| INFO/VP_DAC | Derived copies among called focal alleles |
| INFO/VP_DAN | Number of called focal allele copies |
| INFO/VP_DAF | VP_DAC / VP_DAN |
| FORMAT/VP_DOSAGE | Derived copies for a fully called sample genotype |

By default focal frequency calculations use every sample except listed
outgroups. Use `--ingroup-list ingroups.txt` to choose a subset. With INFO/FASTA
evidence, explicitly exclude outgroups if they are present in the input.
VP_DOSAGE is annotated for all samples; VP_DAC/DAN/DAF refer only to the focal
subset.

For a resolved A/G SNP with AA=G, GT=0/0 has VP_DOSAGE=2, GT=0/1 has 1 and
GT=1/1 has 0. Haploid dosage ranges from 0 to 1; polyploid dosage follows the
genotype ploidy. Partial GT contributes its observed alleles to focal
VP_DAC/DAN, but VP_DOSAGE is missing. Zero called focal copies means VP_DAN=0
and no VP_DAF. Focal counts do not apply additional DP/GQ filters; filter focal
genotypes upstream if required.

Other statuses include NO_AA, INVALID_AA, LOW_CONFIDENCE_AA,
TOO_FEW_OUTGROUPS, TOO_MISSING_OUTGROUPS, LOW_SUPPORT, FASTA_MISSING,
AA_NOT_IN_ALLELES and UNSUPPORTED_VARIANT. An old AA is saved in VP_INPUT_AA
and removed from AA if no current call was obtained, preventing accidental
reuse of stale evidence.

The tool does not infer ancestral posterior probabilities, reconstruct
orthology, or estimate fitness load. Functional annotations must be interpreted
relative to the derived allele; an ALT-effect label is not automatically an
effect label for a derived REF allele.

## Filtering, compression and streaming

```sh
vcf-polar input.vcf.gz --drop-unresolved -o resolved.vcf.gz --index
vcf-polar input.bcf -o annotated.bcf --index --threads 4
bcftools view -Ou input.vcf.gz | vcf-polar - -o annotated.vcf.gz
```

Input can be plain VCF, ordinary gzip VCF, BGZF VCF or BCF. Named .gz/.bgz
outputs use BGZF; .bcf uses BCF. Other output names and stdout use plain VCF.
Logs go to stderr. `--index` requires a named compressed output with sorted
records; it builds a TBI for VCF or CSI for BCF. `--threads` controls HTSlib I/O,
not multiprocessing of individual sites.

Input headers must declare contigs and all INFO/FORMAT tags used by the input.
For a VCF missing contig declarations, add them with
`bcftools reheader --fai reference.fa.fai` before running. Named outputs and
reports are staged next to their destination files. Processing/index failures
leave existing named results untouched; stdout is streamed and cannot be
rolled back. Input, output, reports and generated index paths must be distinct.
Replacing an output removes its previous indexes to prevent stale index reuse.

## Migration from v1.1

**Breaking defaults:** the default mode is now `annotate` and unresolved sites
are kept. Consensus now counts all quality-passing allele copies, including
heterozygotes, with DP/GQ/missingness filters. Version 2 results must not be
treated as an exact reproduction of the original paper's v1.1 analysis.

The original script remains available at the immutable
[v1.1 source commit](https://github.com/LeoHongboWANG/sticklebacks_PAV/blob/d3a690c445cd3846c758500e5cadb7b36ebf3b11/02.deleterious_mutation/VCF-Polar.py).

Explicit legacy exports remain available:

```sh
vcf-polar input.vcf.gz --mode ref -o ancestral-coded.vcf
vcf-polar input.vcf.gz --mode alt -o derived-coded.vcf
```

These are **analysis-only allele encodings**, not compatible with the original
reference FASTA. Do not use them as inputs for reference validation, VEP/SnpEff
annotation or read alignment against that FASTA. In these modes unresolved
records are dropped by default. `--keep-all` explicitly retains them unchanged
in REF/ALT/GT, so consumers must filter VP_STATUS.

For a swapped biallelic SNP, GT and PGT are recoded, Number=R vectors are
permuted and Number=G likelihood/posterior arrays are reversed. Diploid phase,
haploid GT, unphased polyploid GT and fully phased polyploid GT are supported.
Mixed polyploid phasing is rejected in legacy mode and retained by annotate
mode. Missing vectors stay missing.

For records with GT, AC/AN/AF are recomputed from all called alleles across all
samples in the file. For site-only records, AC is transformed using AN and AF
is complemented when available. MAF is preserved. Scalar DP/GQ/PS/FT/MIN_DP and
derived annotations are retained. Unknown annotations, including ALT-specific
ANN/CSQ or DS without a safe general transformation, are removed **only on
swapped records**. Names are listed in VP_DROPPED, stderr and JSON field counts.
VP_SWAPPED and VP_ORIG_REF/ALT identify the transformation.

The old `--by-chrom`, `--batch-size` and `--tmp-dir` options are accepted with a
deprecation warning. All modes now stream through a single processing path;
there are no nested process pools.

## Tests

```sh
python -m pip install ".[test]"
python -m pytest -q
```

Tests cover reference preservation, ancestral decisions, heterozygous
outgroups, missingness and quality filters, haploid/polyploid and phased GT,
legacy AD/PL/GL/GP transformations, ordinary gzip migration, working BGZF/BCF
indexes, input protection, stdin/stdout and the bundled example. GitHub Actions
runs the suite with Python 3.10 and 3.12.
