"""Reference-preserving ancestral-allele annotation for biallelic SNPs.

VCF-Polar v2 succeeds the v1.1 script in 02.deleterious_mutation.
The original workflow acknowledged Ari's pipeline and Biostars post 173338.
"""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import ExitStack
from dataclasses import dataclass, replace
import json
import logging
import math
import os
from pathlib import Path
import sys
import tempfile

import pysam

__version__ = "2.0.0"
LOG = logging.getLogger("vcf-polar")

INFO_FIELDS = {
    "AA": (1, "String", "Inferred ancestral allele; only PASS sites are resolved by VCF-Polar"),
    "VP_STATUS": (1, "String", "VCF-Polar decision; PASS means resolved"),
    "VP_SOURCE": (1, "String", "Ancestral evidence source"),
    "VP_INPUT_AA": (1, "String", "Original input ancestral annotation before VCF-Polar"),
    "VP_N_OUT": (1, "Integer", "Outgroup individuals passing genotype and quality filters"),
    "VP_N_OUT_ALLELES": (1, "Integer", "Called outgroup allele copies used for consensus"),
    "VP_SUPPORT": (1, "Float", "Largest outgroup allele-copy fraction; not a posterior probability"),
    "VP_FALLBACK": (1, "String", "Primary failure reason when explicit INFO fallback was used"),
    "VP_DAC": (1, "Integer", "Derived allele copies among called ingroup alleles"),
    "VP_DAN": (1, "Integer", "Called ingroup allele copies used for VP_DAF"),
    "VP_DAF": (1, "Float", "Derived allele frequency among called ingroup alleles"),
    "VP_SWAPPED": (0, "Flag", "REF and ALT swapped for an analysis-only legacy export"),
    "VP_ORIG_REF": (1, "String", "REF before a legacy allele swap"),
    "VP_ORIG_ALT": (1, "String", "ALT before a legacy allele swap"),
    "VP_DROPPED": (".", "String", "Annotations removed because a legacy swap cannot safely transform them"),
}
DOSAGE_FIELD = "VP_DOSAGE"


@dataclass(frozen=True)
class Call:
    aa: str | None
    status: str
    source: str
    n_samples: int = 0
    n_alleles: int = 0
    support: float | None = None
    fallback: str | None = None


def base_call(value, source: str, allow_lowercase: bool) -> Call:
    if isinstance(value, tuple):
        if len(value) != 1:
            return Call(None, "INVALID_AA", source)
        value = value[0]
    if value is None:
        return Call(None, "NO_AA", source)
    base = str(value).split("|", 1)[0]
    if base.upper() not in {"A", "C", "G", "T"}:
        return Call(None, "NO_AA" if base.upper() in {"N", "."} else "INVALID_AA", source)
    if base.islower() and not allow_lowercase:
        return Call(None, "LOW_CONFIDENCE_AA", source)
    return Call(base.upper(), "PASS", source)


def outgroup_call(record, names: list[str], args) -> Call:
    counts = [0, 0]
    n_samples = 0
    for name in names:
        sample = record.samples[name]
        gt = sample.get("GT")
        if not gt or any(a not in (0, 1) for a in gt):
            continue
        if any(
            limit > 0 and (sample.get(tag) is None or sample[tag] < limit)
            for tag, limit in (("DP", args.min_dp), ("GQ", args.min_gq))
        ):
            continue
        n_samples += 1
        for allele in gt:
            counts[allele] += 1
    n_alleles = sum(counts)
    support = max(counts) / n_alleles if n_alleles else None
    result = Call(None, "PASS", "outgroup", n_samples, n_alleles, support)
    if n_samples < args.min_outgroups:
        return replace(result, status="TOO_FEW_OUTGROUPS")
    if 1 - n_samples / len(names) > args.max_outgroup_missing + 1e-12:
        return replace(result, status="TOO_MISSING_OUTGROUPS")
    if support is None or support + 1e-12 < args.consensus_threshold:
        return replace(result, status="LOW_SUPPORT")
    allele = 0 if counts[0] > counts[1] else 1
    return replace(result, aa=record.alleles[allele].upper())


def infer_call(record, source, outgroups, fasta, args) -> Call:
    if source == "outgroup":
        call = outgroup_call(record, outgroups, args)
    elif source == "fasta":
        try:
            if record.start >= fasta.get_reference_length(record.contig):
                raise ValueError("position outside FASTA")
            value = fasta.fetch(record.contig, record.start, record.start + 1)
            call = base_call(value, "fasta", args.allow_lowercase)
        except (KeyError, ValueError):
            call = Call(None, "FASTA_MISSING", "fasta")
    else:
        call = base_call(record.info.get(args.aa_tag), "info", args.allow_lowercase)
    if call.status != "PASS" and args.fallback_info and source != "info":
        fallback = base_call(record.info.get(args.aa_tag), "info_fallback", args.allow_lowercase)
        if fallback.status == "PASS":
            call = replace(
                fallback, fallback=call.status, n_samples=call.n_samples,
                n_alleles=call.n_alleles, support=call.support,
            )
    if call.aa and call.aa not in {a.upper() for a in record.alleles}:
        call = replace(call, status="AA_NOT_IN_ALLELES")
    return call


def add_definition(header, collection, tag, number, kind, description):
    definitions = getattr(header, collection)
    if tag in definitions:
        old = definitions[tag]
        allowed_types = {kind, "Character"} if tag == "AA" else {kind}
        if old.number != number or old.type not in allowed_types:
            raise ValueError(f"Incompatible existing {collection}/{tag} header; rename or remove that tag first")
        return
    group = "INFO" if collection == "info" else "FORMAT"
    header.add_line(f'##{group}=<ID={tag},Number={number},Type={kind},Description="{description}">')


def output_header(header, mode):
    header = header.copy()
    for tag, (number, kind, description) in INFO_FIELDS.items():
        add_definition(header, "info", tag, number, kind, description)
    add_definition(
        header, "formats", DOSAGE_FIELD, 1, "Integer",
        "Derived allele copies for a fully called genotype; not ALT dosage",
    )
    if mode != "annotate":
        for tag, number, kind, description in (
            ("AC", "A", "Integer", "ALT allele count among called sample alleles"),
            ("AN", 1, "Integer", "Total called sample allele copies"),
            ("AF", "A", "Float", "ALT allele frequency among called sample alleles"),
        ):
            add_definition(header, "info", tag, number, kind, description)
    header.add_meta("vcf_polar_version", value=__version__)
    header.add_meta("vcf_polar_mode", value=mode)
    header.add_meta("vcf_polar_reference_compatible", value="true" if mode == "annotate" else "false")
    return header


def annotate_record(record, header, call, ingroups, aa_tag):
    rec = record.copy()
    rec.translate(header)
    old_aa = rec.info.get(aa_tag)
    # Do not leave a stale AA usable as an implicit fallback after inference failed.
    for tag in INFO_FIELDS:
        if tag in rec.info:
            del rec.info[tag]
    if old_aa is not None:
        rec.info["VP_INPUT_AA"] = "|".join(map(str, old_aa)) if isinstance(old_aa, tuple) else str(old_aa)
    rec.info["VP_STATUS"] = call.status
    rec.info["VP_SOURCE"] = call.source
    rec.info["VP_N_OUT"] = call.n_samples
    rec.info["VP_N_OUT_ALLELES"] = call.n_alleles
    if call.aa is not None:
        rec.info["AA"] = call.aa
    if call.support is not None:
        rec.info["VP_SUPPORT"] = call.support
    if call.fallback:
        rec.info["VP_FALLBACK"] = call.fallback

    ancestral_idx = None
    if call.status == "PASS":
        ancestral_idx = [a.upper() for a in rec.alleles].index(call.aa)
    derived_count = called_count = 0
    for name, sample in rec.samples.items():
        gt = sample.get("GT")
        complete = bool(gt) and all(a in (0, 1) for a in gt)
        sample[DOSAGE_FIELD] = (
            sum(a != ancestral_idx for a in gt) if complete and ancestral_idx is not None else None
        )
        if name in ingroups and gt and ancestral_idx is not None:
            for allele in gt:
                if allele in (0, 1):
                    called_count += 1
                    derived_count += allele != ancestral_idx
    if ancestral_idx is not None:
        rec.info["VP_DAC"] = derived_count
        rec.info["VP_DAN"] = called_count
        if called_count:
            rec.info["VP_DAF"] = derived_count / called_count
    return rec


def reverse_vector(value, expected, label):
    values = value if isinstance(value, tuple) else (value,)
    if values == (None,):
        return (None,) * expected if expected is not None else values
    if expected is not None and len(values) != expected:
        raise ValueError(f"Cannot swap malformed {label}: expected {expected} values, got {len(values)}")
    return values[::-1]


def swap_record(rec, stats):
    """Permute supported biallelic annotations; explicitly discard unsafe ones."""
    keys = list(rec.format)
    if "GT" in keys:
        gt_index = keys.index("GT")
        for raw in str(rec).rstrip("\n").split("\t")[9:]:
            gt = raw.split(":")[gt_index]
            if "/" in gt and "|" in gt:
                raise ValueError("Legacy export cannot preserve mixed polyploid phasing; use --mode annotate")
    samples = {
        name: (dict(sample.items()), sample.phased)
        for name, sample in rec.samples.items()
    }
    old_info = dict(rec.info.items())
    orig_ref, orig_alt = rec.alleles
    rec.alleles = (orig_alt, orig_ref)
    dropped = []
    safe_info = set(INFO_FIELDS) | {"AN", "DP", "MAF"}
    for tag, value in old_info.items():
        if tag in {"AC", "AF"}:
            continue
        number = rec.header.info[tag].number
        if tag in safe_info:
            continue
        if number in {"R", "G"}:
            rec.info[tag] = reverse_vector(value, 2 if number == "R" else None, f"INFO/{tag}")
        else:
            del rec.info[tag]
            dropped.append(f"INFO_{tag}")
            stats["dropped_info_fields"][tag] += 1
    safe_format = {"DP", "GQ", "PS", "FT", "MIN_DP", DOSAGE_FIELD}
    for name, (values, phased) in samples.items():
        gt = values.get("GT")
        if gt:
            rec.samples[name]["GT"] = tuple(1 - a if a in (0, 1) else a for a in gt)
            rec.samples[name].phased = phased
    for tag in keys:
        if tag == "GT" or tag in safe_format:
            continue
        number = rec.header.formats[tag].number
        if number in {"R", "G"}:
            for name, (values, _) in samples.items():
                gt = values.get("GT")
                expected = 2 if number == "R" else len(gt) + 1 if gt else None
                rec.samples[name][tag] = reverse_vector(values[tag], expected, f"FORMAT/{tag}")
        elif tag == "PGT":
            for name, (values, _) in samples.items():
                pgt = values[tag]
                if pgt is not None and any(c not in "01./|" for c in pgt):
                    raise ValueError("Invalid biallelic PGT during legacy swap")
                rec.samples[name][tag] = (
                    "".join("1" if c == "0" else "0" if c == "1" else c for c in pgt)
                    if pgt is not None else None
                )
        else:
            del rec.format[tag]
            dropped.append(f"FORMAT_{tag}")
            stats["dropped_format_fields"][tag] += 1
    if "GT" in rec.format and rec.samples:
        called = [
            a for sample in rec.samples.values()
            for a in (sample.get("GT") or ()) if a in (0, 1)
        ]
        rec.info["AC"] = (sum(called),)
        rec.info["AN"] = len(called)
        if called:
            rec.info["AF"] = (sum(called) / len(called),)
        elif "AF" in rec.info:
            del rec.info["AF"]
    else:
        ac, an, af = old_info.get("AC"), old_info.get("AN"), old_info.get("AF")
        if ac and ac[0] is not None and an is not None:
            rec.info["AC"] = (an - ac[0],)
        elif "AC" in rec.info:
            del rec.info["AC"]
        if af and af[0] is not None:
            rec.info["AF"] = (1 - af[0],)
    rec.info["VP_SWAPPED"] = True
    rec.info["VP_ORIG_REF"] = orig_ref
    rec.info["VP_ORIG_ALT"] = orig_alt
    if dropped:
        rec.info["VP_DROPPED"] = tuple(sorted(dropped))
    return rec


def sample_list(path, available):
    if not path:
        return []
    with open(path, encoding="utf-8") as handle:
        names = [s.strip() for s in handle if s.strip() and not s.lstrip().startswith("#")]
    if not names:
        raise ValueError(f"Empty sample list: {path}")
    if len(names) != len(set(names)):
        raise ValueError(f"Duplicate names in sample list: {path}")
    missing = sorted(set(names) - set(available))
    if missing:
        raise ValueError(f"Samples not found in VCF: {', '.join(missing)}")
    return names


def parser():
    p = argparse.ArgumentParser(
        description="Annotate ancestral alleles and derived dosage without changing reference alleles.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("vcf", help="Input VCF, VCF.gz or BCF; '-' reads stdin")
    p.add_argument("-o", "--output", default="-", help="Output VCF, BGZF .vcf.gz or BCF; '-' writes stdout")
    p.add_argument("--version", action="version", version=f"VCF-Polar {__version__}")
    p.add_argument("--mode", choices=["annotate", "ref", "alt"], default="annotate",
                   help="ref/alt are analysis-only legacy exports, incompatible with the original FASTA")
    p.add_argument("--source", choices=["auto", "info", "outgroup", "fasta"], default="auto")
    p.add_argument("--aa-tag", default="AA", help="Input INFO ancestral-allele tag; output uses AA")
    p.add_argument("--anc-list", help="Outgroup sample names, one per line")
    p.add_argument("--ingroup-list", help="Samples for VP_DAC/DAN/DAF; default excludes selected outgroups")
    p.add_argument("--ancestral-fasta", help="Indexed ancestral FASTA already in the input reference coordinates")
    p.add_argument("--consensus-threshold", type=float, default=0.8,
                   help="Minimum outgroup allele-copy support, greater than 0.5; not a posterior")
    p.add_argument("--min-outgroups", type=int, default=1)
    p.add_argument("--max-outgroup-missing", type=float, default=0.5,
                   help="Maximum fraction of outgroups failing GT/DP/GQ filters")
    p.add_argument("--min-dp", type=int, default=3, help="Outgroup DP filter; 0 disables, missing fails if enabled")
    p.add_argument("--min-gq", type=float, default=20, help="Outgroup GQ filter; 0 disables, missing fails if enabled")
    p.add_argument("--allow-lowercase", action="store_true", help="Allow lowercase ancestral calls")
    p.add_argument("--fallback-info", action="store_true", help="Explicitly allow INFO fallback if primary evidence fails")
    filtering = p.add_mutually_exclusive_group()
    filtering.add_argument("--drop-unresolved", dest="keep_all", action="store_false",
                           help="Drop unresolved sites; default only in legacy modes")
    filtering.add_argument("--keep-all", dest="keep_all", action="store_true",
                           help="Retain unresolved sites; default in annotate mode")
    p.set_defaults(keep_all=None)
    p.add_argument("--report", help="Write JSON decision counts and run settings")
    p.add_argument("--sites-report", help="Write per-site decision TSV, including dropped sites")
    p.add_argument("--threads", type=int, default=1, help="HTSlib compression/decompression threads")
    p.add_argument("--index", action="store_true", help="Index sorted BGZF VCF or BCF output")
    p.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], default="INFO")
    # Compatibility switches no longer select a separate processing implementation.
    p.add_argument("--by-chrom", action="store_true", help=argparse.SUPPRESS)
    p.add_argument("--batch-size", type=int, default=None, help=argparse.SUPPRESS)
    p.add_argument("--tmp-dir", default=None, help=argparse.SUPPRESS)
    return p


def validate(args):
    if args.output == "stdout":
        args.output = "-"
    for name in ("consensus_threshold", "max_outgroup_missing", "min_gq"):
        if not math.isfinite(getattr(args, name)):
            raise ValueError(f"--{name.replace('_', '-')} must be finite")
    if not 0.5 < args.consensus_threshold <= 1:
        raise ValueError("--consensus-threshold must be greater than 0.5 and at most 1")
    if not 0 <= args.max_outgroup_missing <= 1:
        raise ValueError("--max-outgroup-missing must be between 0 and 1")
    if args.min_outgroups < 1 or args.min_dp < 0 or args.min_gq < 0 or args.threads < 1:
        raise ValueError("Require min-outgroups >= 1, min-dp/min-gq >= 0 and threads >= 1")
    if args.batch_size is not None and args.batch_size < 1:
        raise ValueError("--batch-size must be positive")
    if args.by_chrom or args.batch_size is not None or args.tmp_dir:
        LOG.warning("--by-chrom/--batch-size/--tmp-dir are deprecated; input is now streamed in order")
    if args.source == "auto":
        if args.anc_list and args.ancestral_fasta:
            raise ValueError("Choose --source explicitly when both outgroups and FASTA are supplied")
        args.source = "outgroup" if args.anc_list else "fasta" if args.ancestral_fasta else "info"
    if args.source == "outgroup" and not args.anc_list:
        raise ValueError("--source outgroup requires --anc-list")
    if args.source == "fasta" and not args.ancestral_fasta:
        raise ValueError("--source fasta requires --ancestral-fasta")
    if args.ancestral_fasta and args.source != "fasta":
        raise ValueError("--ancestral-fasta is only used with --source fasta")
    if args.keep_all is None:
        args.keep_all = args.mode == "annotate"
    if args.index and (args.output == "-" or not args.output.endswith((".gz", ".bgz", ".bcf"))):
        raise ValueError("--index requires a named BGZF VCF (.gz/.bgz) or BCF output")
    outputs = [p for p in (args.output, args.report, args.sites_report) if p and p != "-"]
    # Sidecars are reserved even when not requested: stale output indexes are removed.
    sidecars = [args.output + suffix for suffix in (".tbi", ".csi", ".gzi")] if args.output != "-" else []
    if args.report == "-" or args.sites_report == "-":
        raise ValueError("Report files require a path, not stdout")
    inputs = [p for p in (args.vcf, args.anc_list, args.ingroup_list, args.ancestral_fasta) if p and p != "-"]
    if args.ancestral_fasta:
        inputs += [args.ancestral_fasta + ".fai", args.ancestral_fasta + ".gzi"]
    def same(a, b):
        return Path(a).resolve() == Path(b).resolve() or (
            os.path.exists(a) and os.path.exists(b) and os.path.samefile(a, b)
        )
    for i, output in enumerate(outputs):
        if any(same(output, path) for path in inputs + outputs[:i]):
            raise ValueError("Input, output and report paths must be distinct")
        if not Path(output).parent.is_dir():
            raise ValueError(f"Output directory does not exist: {Path(output).parent}")
    if any(same(sidecar, path) for sidecar in sidecars for path in inputs + outputs):
        raise ValueError("Output index paths must not overlap inputs or reports")
    if args.mode != "annotate":
        LOG.warning(
            "Legacy export changes REF/ALT and is incompatible with the original reference FASTA. "
            "Unsupported annotations on swapped records will be removed and reported."
        )


def run(args):
    stats = {
        "version": __version__,
        "settings": {k: v for k, v in vars(args).items() if k not in {"log_level"}},
        "input_records": 0, "output_records": 0, "resolved_records": 0,
        "swapped_records": 0, "dropped_unresolved": 0,
        "status_counts": Counter(), "source_counts": Counter(),
        "dropped_info_fields": Counter(), "dropped_format_fields": Counter(),
    }
    staged = []
    def stage(path):
        suffix = "".join(Path(path).suffixes)
        fd, name = tempfile.mkstemp(prefix=".vcf-polar-", suffix=suffix, dir=Path(path).parent)
        os.close(fd)
        staged.append((name, path))
        return name
    try:
        with ExitStack() as stack:
            # A file handle explicitly requests sequential HTSlib input. This also
            # handles ordinary gzip from v1 without a seek/tell or full-file copy.
            raw_input = stack.enter_context(open(args.vcf, "rb")) if args.vcf != "-" else "-"
            reader = stack.enter_context(pysam.VariantFile(raw_input, threads=args.threads))
            if not reader.header.contigs:
                raise ValueError("Input requires ##contig declarations; add them with bcftools reheader --fai reference.fa.fai")
            available = list(reader.header.samples)
            if args.source == "info" and args.aa_tag not in reader.header.info:
                raise ValueError(f"Input header has no INFO/{args.aa_tag}; supply outgroups or an ancestral FASTA")
            outgroups = sample_list(args.anc_list, available)
            ingroups = set(sample_list(args.ingroup_list, available)) if args.ingroup_list else set(available) - set(outgroups)
            if set(outgroups) & ingroups:
                raise ValueError("Ingroup and outgroup lists must not overlap")
            fasta = None
            if args.source == "fasta":
                fasta = stack.enter_context(pysam.FastaFile(args.ancestral_fasta))
                if reader.header.contigs and not set(reader.header.contigs) & set(fasta.references):
                    raise ValueError("Input VCF and ancestral FASTA have no matching contig names")
            stats["outgroups"] = outgroups
            stats["ingroups"] = [n for n in available if n in ingroups]
            header = output_header(reader.header, args.mode)
            target = stage(args.output) if args.output != "-" else "-"
            mode = "wb" if args.output.endswith(".bcf") else "wz" if args.output.endswith((".gz", ".bgz")) else "w"
            writer = stack.enter_context(pysam.VariantFile(target, mode, header=header, threads=args.threads))
            sites = stack.enter_context(open(stage(args.sites_report), "w", encoding="utf-8")) if args.sites_report else None
            if sites:
                sites.write("CHROM\tPOS\tREF\tALT\tAA\tSTATUS\tSOURCE\tN_OUT\tN_OUT_ALLELES\tSUPPORT\tWRITTEN\tSWAPPED\n")
            for record in reader:
                if record.contig not in header.contigs:
                    raise ValueError(f"Contig {record.contig} is not declared in the input VCF header")
                stats["input_records"] += 1
                supported = len(record.alleles or ()) == 2 and all(
                    len(a) == 1 and a.upper() in {"A", "C", "G", "T"} for a in record.alleles
                )
                call = infer_call(record, args.source, outgroups, fasta, args) if supported else Call(None, "UNSUPPORTED_VARIANT", args.source)
                stats["status_counts"][call.status] += 1
                stats["source_counts"][call.source] += 1
                resolved = call.status == "PASS"
                stats["resolved_records"] += resolved
                written = resolved or args.keep_all
                swapped = False
                if written:
                    rec = annotate_record(record, header, call, ingroups, args.aa_tag)
                    swapped = resolved and (
                        (args.mode == "ref" and call.aa != record.ref.upper()) or
                        (args.mode == "alt" and call.aa != record.alts[0].upper())
                    )
                    if swapped:
                        rec = swap_record(rec, stats)
                    writer.write(rec)
                    stats["output_records"] += 1
                    stats["swapped_records"] += swapped
                else:
                    stats["dropped_unresolved"] += 1
                if sites:
                    values = (
                        record.contig, record.pos, record.ref, ",".join(record.alts or (".",)),
                        call.aa or ".", call.status, call.source, call.n_samples, call.n_alleles,
                        "." if call.support is None else f"{call.support:.6g}", int(written), int(swapped),
                    )
                    sites.write("\t".join(map(str, values)) + "\n")
        # Build an index before replacing any existing output, so sort/index failures are atomic.
        if args.index:
            if args.output.endswith(".bcf"):
                from pysam import bcftools
                bcftools.index("-f", target)
                staged.append((target + ".csi", args.output + ".csi"))
            else:
                pysam.tabix_index(target, preset="vcf", force=True)
                staged.append((target + ".tbi", args.output + ".tbi"))
        if args.report:
            with open(stage(args.report), "w", encoding="utf-8") as handle:
                json.dump(stats, handle, indent=2, sort_keys=True, allow_nan=False)
                handle.write("\n")
        if args.output != "-":
            for suffix in (".tbi", ".csi", ".gzi"):
                old_index = args.output + suffix
                if os.path.exists(old_index):
                    os.unlink(old_index)
        for temporary, final in staged:
            os.replace(temporary, final)
    finally:
        for temporary, _ in staged:
            for leftover in (temporary, temporary + ".tbi", temporary + ".csi"):
                if os.path.exists(leftover):
                    os.unlink(leftover)
    LOG.info(
        "Processed %s sites; resolved %s; wrote %s; swapped %s; dropped unresolved %s",
        stats["input_records"], stats["resolved_records"], stats["output_records"],
        stats["swapped_records"], stats["dropped_unresolved"],
    )
    if stats["dropped_info_fields"] or stats["dropped_format_fields"]:
        LOG.warning("Removed annotations on swapped records: INFO=%s FORMAT=%s",
                    dict(stats["dropped_info_fields"]), dict(stats["dropped_format_fields"]))
    return stats


def main(argv=None):
    p = parser()
    args = p.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(levelname)s: %(message)s")
    try:
        validate(args)
        run(args)
    except (ValueError, TypeError, OSError, KeyError, pysam.utils.SamtoolsError) as exc:
        LOG.error("%s", exc, exc_info=args.log_level == "DEBUG")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
