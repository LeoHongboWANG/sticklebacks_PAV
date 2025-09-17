#!/usr/bin/env python3

"""
    Advanced VCF Recoder for Ancestral Allele Polarization
    
    Date: 2025-03-22
    Version: 1.1

    Description:
    This script processes VCF files for ancestral allele polarization. Copilot refined the code and format. It references Ari pipeline # slightly modified code from https://www.biostars.org/p/173338/
    cat > biostars_173338 << 'EOF'
    while(<>) {
        if(/^#/) {
            print $_;
            next;
        }
        chomp;
        my($chr,$site,$rs,$ref,$alt,$info) = (split(/\t/,$_))[0,1,2,3,4,7];
        my $ancestral_allele;
        $ancestral_allele = (split(/AA=/,$info))[1];
        $ancestral_allele = (split(/[\|;]/,$ancestral_allele))[0];
        $ancestral_allele = uc($ancestral_allele);
        if($ancestral_allele !~ /[ACGT]/) {
            print STDERR "Ancestral allele ($ancestral_allele) is not ACGT: $info\n";
        next;
        }
        if($alt =~ /,/ || $alt !~ /[ACGT]/) {
        print STDERR "Alternative allele ($alt) is not unique: $info\n";
        next;
        }
        if($ancestral_allele eq $ref) { #ok as it is
            print "$_\n";
            next;
        } else { #Switch alleles
            my @tmp = split(/;/,$info);
            foreach my $t ( @tmp ) {
                next if($t !~ /AF/);
                my($key,$val) = split(/=/,$t);
                my $new_AF = 1-$val;
                $_ =~ s/$key=$val/$key=$new_AF/;
            }
            my($chr,$site,$rs,$ref,$alt,$qual,$filter,$info,$format,@tmp) = split(/\t/,$_);
        my $swp = $ref;
        $ref = $alt; 
        $alt = $swp;
            my $str = join "\t",$chr,$site,$rs,$ref,$alt,$qual,$filter,$info,$format,"";
        foreach my $t ( @tmp ) {
                if($t =~ /^[01][\/][01]/) {
                    my($a,$b,$c,$d,$e) = (split(/[\/:,]/,$t));
                    if($a == 0)    {    $a = 1;    }
                    elsif($a == 1) {    $a = 0;    }
                    if($b == 0)    {    $b = 1;    }
                    elsif($b == 1) {    $b = 0;    }
                    $t = "$a/$b:$d,$c:$e";
                }
                $str .= "$t\t";
            }
            $str =~ s/\t$//;
            print "$str\n";
        }
    }
    EOF
    .
    """
import sys
import re
import argparse
import gzip
import multiprocessing
import logging
import os
import tempfile
from typing import List, Optional, Dict, Set, Tuple
from tqdm import tqdm

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

anc_sample_indices = []
cons_threshold = 0.8

def init_globals(anc_indices_arg, threshold_arg):
    global anc_sample_indices, cons_threshold
    anc_sample_indices = anc_indices_arg
    cons_threshold = threshold_arg

def validate_vcf_fields(fields: List[str]) -> bool:
    required_columns = 8
    return len(fields) >= required_columns

def recalculate_af(info: str) -> str:
    pattern = re.compile(r"(AF=)([0-9.eE+-]+)")
    def repl(match):
        af_val = float(match.group(2))
        return f"{match.group(1)}{max(0, min(1, 1 - af_val)):.6g}"
    return re.sub(pattern, repl, info)

def recode_genotype(genotype_field: str) -> str:
    if genotype_field in {"./.", ".|.", "."} or not genotype_field:
        return genotype_field
    
    parts = genotype_field.split(":", 1)
    g = parts[0]
    
    sep = "/" if "/" in g else "|" if "|" in g else None
    
    if sep:
        alleles = g.split(sep)
        swapped = [
            ("1" if a == "0" else "0") if a in {"0", "1"} else a 
            for a in alleles
        ]
        new_g = sep.join(swapped)
        return new_g + (":" + parts[1] if len(parts) > 1 else "")
    
    return genotype_field

def get_gt_indices(format_field: str) -> Optional[int]:
    if not format_field:
        return None
    
    fields = format_field.split(":")
    for i, field in enumerate(fields):
        if field == "GT":
            return i
    
    return None

def compute_consensus(fields: List[str], anc_indices: List[int]) -> Optional[str]:
    ref_allele = fields[3]
    alt_allele = fields[4]
    
    if "," in alt_allele or not re.match(r"^[ACGTacgt]+$", alt_allele):
        return None
    
    ref_count = 0
    alt_count = 0
    total = 0
    
    format_field = fields[8] if len(fields) > 8 else ""
    gt_idx = get_gt_indices(format_field)
    
    if gt_idx is None:
        logging.debug("FORMAT field does not contain GT")
        return None
    
    for idx in anc_indices:
        if idx < len(fields):
            genotype_field = fields[idx]
            if not genotype_field or genotype_field == ".":
                continue
                
            genotype_parts = genotype_field.split(":")
            if len(genotype_parts) <= gt_idx:
                continue
                
            genotype = genotype_parts[gt_idx]
            
            if "/" in genotype:
                sep = "/"
            elif "|" in genotype:
                sep = "|"
            else:
                continue
                
            alleles = genotype.split(sep)
            
            if any(a == "." for a in alleles):
                continue
                
            if all(a == "0" for a in alleles):
                ref_count += 1
                total += 1
            elif all(a == "1" for a in alleles):
                alt_count += 1
                total += 1
    
    if total == 0:
        return None
        
    if (ref_count / total) >= cons_threshold:
        return ref_allele
    elif (alt_count / total) >= cons_threshold:
        return alt_allele
    else:
        return None

def update_info_with_aa(info: str, ancestral_allele: str) -> str:
    info = re.sub(r";?AA=[^;]+", "", info)
    
    if not info:
        return f"AA={ancestral_allele}"
    else:
        return f"{info};AA={ancestral_allele}"

def recode_line(line: str, mode: str, skip_no_consensus: bool = True) -> Optional[str]:
    if line.startswith("#"):
        return line

    fields = line.strip().split("\t")
    if not validate_vcf_fields(fields):
        logging.warning(f"Invalid VCF line format: {line}")
        return None

    ref, alt, info = fields[3], fields[4], fields[7]
    
    if "," in alt or not re.match(r"^[ACGTacgt]+$", alt) or not re.match(r"^[ACGTacgt]+$", ref):
        logging.debug(f"Skipping complex variant: {line}")
        if not skip_no_consensus:
            return line
        return None

    consensus_allele = None
    global anc_sample_indices
    
    if anc_sample_indices:
        consensus_allele = compute_consensus(fields, anc_sample_indices)
        if consensus_allele:
            logging.debug(f"Consensus allele determined as {consensus_allele} from ancestral samples")
            info = update_info_with_aa(info, consensus_allele)
            fields[7] = info
    
    if not consensus_allele:
        match = re.search(r"AA=([ACGTacgt]+)", info)
        if match:
            consensus_allele = match.group(1).upper()
            logging.debug(f"Using AA={consensus_allele} from INFO field")

    if not consensus_allele:
        logging.debug(f"No ancestral allele found for position {fields[0]}:{fields[1]}")
        if skip_no_consensus:
            return None
        return line

    if consensus_allele not in {ref, alt}:
        logging.debug(f"Ancestral allele {consensus_allele} not in REF/ALT ({ref}/{alt})")
        if skip_no_consensus:
            return None
        return line

    if mode == "ref":
        swap = (consensus_allele != ref)
    else:
        swap = (consensus_allele != alt)

    if not swap:
        return "\t".join(fields)
    
    new_info = recalculate_af(info)
    new_ref, new_alt = alt, ref
    new_fields = fields[:3] + [new_ref, new_alt] + fields[5:7] + [new_info] + fields[8:]
    
    if len(new_fields) > 9:
        new_fields[9:] = [recode_genotype(gt) for gt in new_fields[9:]]
        
    return "\t".join(new_fields)

def process_chunk(chunk_id: int, lines: List[str], mode: str, skip_no_consensus: bool) -> List[Optional[str]]:
    results = []
    for line in lines:
        result = recode_line(line, mode, skip_no_consensus)
        if result is not None:
            results.append(result)
    
    return results

def extract_chromosome_data(infile_path: str, header_lines: List[str]) -> Dict[str, List[str]]:
    chroms = {}
    
    open_func = gzip.open if infile_path.endswith('.gz') else open
    mode = 'rt' if infile_path.endswith('.gz') else 'r'
    
    with open_func(infile_path, mode) as infile:
        for line in infile:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t', 2)
            if len(fields) >= 1:
                chrom = fields[0]
                if chrom not in chroms:
                    chroms[chrom] = []
                chroms[chrom].append(line.rstrip())
    
    return chroms

def process_chrom_data(chrom: str, data_lines: List[str], mode: str, anc_indices: List[int], 
                      threshold: float, threads_per_chrom: int, skip_no_consensus: bool) -> str:
    fd, temp_path = tempfile.mkstemp(suffix=f"_{chrom}.vcf")
    os.close(fd)
    
    num_threads = max(1, min(threads_per_chrom, len(data_lines) // 1000 + 1))
    
    with multiprocessing.Pool(
        num_threads, 
        initializer=init_globals, 
        initargs=(anc_indices, threshold)
    ) as pool:
        chunk_size = max(len(data_lines) // num_threads, 1000)
        chunks = [(i, data_lines[i:i + chunk_size], mode, skip_no_consensus) 
                  for i in range(0, len(data_lines), chunk_size)]
        
        logging.info(f"Processing chromosome {chrom} with {num_threads} threads")
        results = list(tqdm(
            pool.starmap(process_chunk, chunks),
            total=len(chunks),
            unit="chunk",
            desc=f"Chrom {chrom}"
        ))
    
    with open(temp_path, 'w') as outfile:
        for chunk_result in results:
            for line in chunk_result:
                outfile.write(line)
                if not line.endswith('\n'):
                    outfile.write('\n')
    
    logging.info(f"Chromosome {chrom} processing complete. Results in {temp_path}")
    return temp_path

def find_sample_indices(header_lines: List[str], sample_names: Set[str]) -> List[int]:
    for line in header_lines:
        if line.startswith('#CHROM'):
            fields = line.strip().split('\t')
            if len(fields) > 9:
                return [i for i, sample in enumerate(fields[9:], 9) if sample in sample_names]
    return []

def get_vcf_header(infile_path: str) -> List[str]:
    header_lines = []
    
    open_func = gzip.open if infile_path.endswith('.gz') else open
    mode = 'rt' if infile_path.endswith('.gz') else 'r'
    
    with open_func(infile_path, mode) as infile:
        for line in infile:
            line = line.rstrip()
            if line.startswith('#'):
                header_lines.append(line)
            else:
                break
    
    return header_lines

def ensure_aa_header(header_lines: List[str]) -> List[str]:
    aa_header = '##INFO=<ID=AA,Number=1,Type=Character,Description="Ancestral allele">'
    
    for line in header_lines:
        if line.startswith('##INFO=<ID=AA,') and 'Ancestral allele' in line:
            return header_lines
    
    chrom_idx = -1
    last_info_idx = -1
    
    for i, line in enumerate(header_lines):
        if line.startswith('#CHROM'):
            chrom_idx = i
        elif line.startswith('##INFO='):
            last_info_idx = i
    
    if last_info_idx >= 0:
        header_lines.insert(last_info_idx + 1, aa_header)
    elif chrom_idx >= 0:
        header_lines.insert(chrom_idx, aa_header)
    else:
        header_lines.append(aa_header)
    
    return header_lines

def main():
    parser = argparse.ArgumentParser(
        description="Advanced VCF Recoder for ancestral allele polarization"
    )
    parser.add_argument("vcf", help="Input VCF file (supports .gz)")
    parser.add_argument(
        "-o", "--output", 
        help="Output VCF file (default: stdout)", 
        default="stdout"
    )
    parser.add_argument(
        "--mode", 
        choices=["ref", "alt"], 
        default="ref",
        help="'ref' makes REF=ancestral allele, 'alt' makes ALT=ancestral allele"
    )
    parser.add_argument(
        "--threads", 
        type=int, 
        default=max(1, multiprocessing.cpu_count() - 1),
        help="Number of threads (default: CPU cores - 1)"
    )
    parser.add_argument(
        "--log-level", 
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        default='INFO',
        help="Set logging level"
    )
    parser.add_argument(
        "--anc-list", 
        help="File with ancestral sample names (one per line)", 
        default=None
    )
    parser.add_argument(
        "--consensus-threshold", 
        type=float, 
        default=0.8, 
        help="Consensus threshold for ancestral allele determination (default: 0.8)"
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=100000,
        help="Number of variants to process in each batch (default: 100000)"
    )
    parser.add_argument(
        "--by-chrom",
        action="store_true",
        help="Process each chromosome in parallel (recommended for large files)"
    )
    parser.add_argument(
        "--tmp-dir",
        help="Directory for temporary files (default: system temp dir)",
        default=None
    )
    parser.add_argument(
        "--keep-all",
        action="store_false",
        dest="skip_no_consensus",
        help="Keep all variants even without ancestral consensus (default: skip)"
    )
    args = parser.parse_args()

    logging.getLogger().setLevel(getattr(logging, args.log_level))

    if not args.vcf:
        logging.error("No input VCF specified")
        sys.exit(1)
    
    if args.tmp_dir:
        if not os.path.exists(args.tmp_dir):
            os.makedirs(args.tmp_dir)
        tempfile.tempdir = args.tmp_dir

    outfile = sys.stdout if args.output == "stdout" else (
        gzip.open(args.output, "wt") if args.output.endswith(".gz") else open(args.output, "w")
    )
    
    anc_samples_set = set()
    if args.anc_list:
        with open(args.anc_list, "r") as f:
            anc_samples_set = {line.strip() for line in f if line.strip()}
        logging.info(f"Loaded {len(anc_samples_set)} ancestral sample names from {args.anc_list}")
    
    header_lines = get_vcf_header(args.vcf)
    if not header_lines:
        logging.error("Failed to extract header from VCF")
        sys.exit(1)
    
    header_lines = ensure_aa_header(header_lines)
    logging.info("Ensured AA header is present in output VCF")
        
    anc_sample_indices_local = find_sample_indices(header_lines, anc_samples_set)
    if anc_samples_set and anc_sample_indices_local:
        logging.info(f"Found {len(anc_sample_indices_local)} ancestral samples in VCF")
    elif anc_samples_set:
        logging.warning("No ancestral samples from list found in VCF")
        
    for header_line in header_lines:
        outfile.write(header_line)
        if not header_line.endswith('\n'):
            outfile.write('\n')
    
    if args.by_chrom:
        logging.info(f"Extracting data by chromosome from {args.vcf}")
        chrom_data = extract_chromosome_data(args.vcf, header_lines)
        logging.info(f"Found {len(chrom_data)} chromosomes in VCF")
        
        max_chroms = min(len(chrom_data), args.threads)
        threads_per_chrom = max(1, args.threads // max_chroms)
        
        with multiprocessing.Pool(max_chroms) as pool:
            chrom_tasks = [
                (chrom, data, args.mode, anc_sample_indices_local, args.consensus_threshold, 
                 threads_per_chrom, args.skip_no_consensus)
                for chrom, data in chrom_data.items()
            ]
            
            temp_files = list(tqdm(
                pool.starmap(process_chrom_data, chrom_tasks),
                total=len(chrom_tasks),
                unit="chrom",
                desc="Processing chromosomes"
            ))
        
        logging.info(f"Combining results from {len(temp_files)} chromosomes")
        for temp_file in temp_files:
            with open(temp_file, 'r') as infile:
                for line in infile:
                    outfile.write(line)
            
            if os.path.exists(temp_file):
                os.remove(temp_file)
    
    else:
        logging.info("Processing VCF in sequential batches")
        
        def process_vcf_in_batches(infile_path: str, batch_size: int = 100000):
            data_batch = []
            
            open_func = gzip.open if infile_path.endswith('.gz') else open
            mode = 'rt' if infile_path.endswith('.gz') else 'r'
            
            with open_func(infile_path, mode) as infile:
                for line in infile:
                    line = line.rstrip()
                    if line.startswith('#'):
                        continue
                    else:
                        data_batch.append(line)
                        if len(data_batch) >= batch_size:
                            yield data_batch
                            data_batch = []
            
            if data_batch:
                yield data_batch
        
        for i, data_batch in enumerate(process_vcf_in_batches(args.vcf, args.batch_size)):
            num_threads = max(1, min(args.threads, len(data_batch) // 1000 + 1))
            logging.info(f"Processing batch {i+1} with {num_threads} threads ({len(data_batch)} variants)")
            
            chunk_size = max(len(data_batch) // num_threads, 1000)
            chunks = [(i, data_batch[i:i + chunk_size], args.mode, args.skip_no_consensus) 
                      for i in range(0, len(data_batch), chunk_size)]
            
            with multiprocessing.Pool(
                num_threads, 
                initializer=init_globals, 
                initargs=(anc_sample_indices_local, args.consensus_threshold)
            ) as pool:
                results = list(tqdm(
                    pool.starmap(process_chunk, chunks),
                    total=len(chunks),
                    unit="chunk",
                    desc=f"Batch {i+1}"
                ))
            
            for chunk_result in results:
                for line in chunk_result:
                    outfile.write(line)
                    if not line.endswith('\n'):
                        outfile.write('\n')
            
            logging.info(f"Completed batch {i+1}")
    
    logging.info("VCF processing complete")
    
    if outfile is not sys.stdout:
        outfile.close()

if __name__ == '__main__':
    main()
