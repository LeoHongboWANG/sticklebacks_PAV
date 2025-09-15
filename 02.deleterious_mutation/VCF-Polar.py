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

# what's happening
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

# not good but it works
ancestral_sample_indices = []
consensus_thresh = 0.8

def setup_worker_globals(anc_indices_arg, threshold_arg):
    """Set up global variably for worker processes."""
    global ancestral_sample_indices, consensus_thresh
    ancestral_sample_indices = anc_indices_arg
    consensus_thresh = threshold_arg

def check_vcf_line_format(fields: List[str]) -> bool:
    """Make sure we have enough columns in the VCF line."""
    min_cols = 8  # VCF spec requires at least 8 columns
    return len(fields) >= min_cols

def flip_allele_frequency(info_field: str) -> str:

    af_pattern = re.compile(r"(AF=)([0-9.eE+-]+)")
    
    def replace_af(match):
        try:
            original_af = float(match.group(2))
            # Calculate 1 - AF, but keep it within [0,1] bounds
            flipped_af = max(0, min(1, 1 - original_af))
            return f"{match.group(1)}{flipped_af:.6g}"
        except ValueError:
            logging.warning(f"Couldn't parse AF value: {match.group(2)}")
            return match.group(0)  # Return original if parsing fails
    
    return re.sub(af_pattern, replace_af, info_field)

def swap_genotype_alleles(gt_field: str) -> str:
    """
    Swap genotype alleles: 0 becomes 1, 1 becomes 0.
    Need to handle different separators and preserve other format fields.
    """
    # Handle missing data cases
    if gt_field in {"./.", ".|.", "."} or not gt_field:
        return gt_field
    
    try:
        # Split genotype from other format fields (like DP, GQ, etc.)
        format_parts = gt_field.split(":", 1)
        genotype = format_parts[0]
        
        # Figure out what separator is being used
        separator = None
        if "/" in genotype:
            separator = "/"
        elif "|" in genotype:
            separator = "|"
        
        if separator:
            alleles = genotype.split(separator)
            # Swap the alleles - only swap 0s and 1s, leave others alone
            swapped_alleles = []
            for allele in alleles:
                if allele == "0":
                    swapped_alleles.append("1")
                elif allele == "1":
                    swapped_alleles.append("0")
                else:
                    swapped_alleles.append(allele)  # Keep missing (.) or other values
            
            new_genotype = separator.join(swapped_alleles)
            
            # Put it back together with format fields
            if len(format_parts) > 1:
                return new_genotype + ":" + format_parts[1]
            else:
                return new_genotype
        
        # If no separator found, just return as-is
        return gt_field
        
    except Exception as e:
        logging.error(f"Problem swapping genotype {gt_field}: {e}")
        return gt_field  # Return original on error

def find_gt_position(format_str: str) -> Optional[int]:
    """Find where GT is in the FORMAT field."""
    if not format_str:
        return None
    
    format_fields = format_str.split(":")
    for idx, field in enumerate(format_fields):
        if field == "GT":
            return idx
    
    return None  # GT not found

def determine_ancestral_consensus(vcf_fields: List[str], ancestral_indices: List[int]) -> Optional[str]:
    """
    Figure out the ancestral allele from the specified samples.
    We count homozygous REF and ALT calls and see if either meets our threshold.
    """
    ref_allele = vcf_fields[3]
    alt_allele = vcf_fields[4]
    
    # Skip if this is a complex variant (multiple alts, non-ACGT bases)
    if "," in alt_allele or not re.match(r"^[ACGTacgt]+$", alt_allele):
        return None
    
    # Counters for homozygous calls
    ref_homo_count = 0
    alt_homo_count = 0
    total_informative = 0
    
    # Find GT position in FORMAT
    format_field = vcf_fields[8] if len(vcf_fields) > 8 else ""
    gt_position = find_gt_position(format_field)
    
    if gt_position is None:
        logging.debug("No GT field found in FORMAT")
        return None
    
    # Check each ancestral sample
    for sample_idx in ancestral_indices:
        if sample_idx < len(vcf_fields):
            sample_data = vcf_fields[sample_idx]
            if not sample_data or sample_data == ".":
                continue
                
            sample_fields = sample_data.split(":")
            if len(sample_fields) <= gt_position:
                continue
                
            gt = sample_fields[gt_position]
            
            # Determine separator
            if "/" in gt:
                sep = "/"
            elif "|" in gt:
                sep = "|"
            else:
                continue
                
            alleles = gt.split(sep)
            
            # Skip if any allele is missing
            if any(a == "." for a in alleles):
                continue
                
            # Only count homozygous calls for consensus
            if all(a == "0" for a in alleles):
                ref_homo_count += 1
                total_informative += 1
            elif all(a == "1" for a in alleles):
                alt_homo_count += 1
                total_informative += 1
            # Heterozygous calls don't contribute to consensus
    
    if total_informative == 0:
        return None
        
    # Check if either allele meets threshold
    ref_proportion = ref_homo_count / total_informative
    alt_proportion = alt_homo_count / total_informative
    
    if ref_proportion >= consensus_thresh:
        return ref_allele
    elif alt_proportion >= consensus_thresh:
        return alt_allele
    else:
        return None

def add_aa_to_info(info_field: str, ancestral_allele: str) -> str:
    """Add or update the AA field in INFO."""
    # Remove any existing AA field first
    info_field = re.sub(r";?AA=[^;]+", "", info_field)
    
    # Add the new AA field
    if not info_field:
        return f"AA={ancestral_allele}"
    else:
        return f"{info_field};AA={ancestral_allele}"

def process_vcf_line(line: str, mode: str, skip_without_consensus: bool = True) -> Optional[str]:
    """
    Process a single VCF line - this is where the magic happens.
    We determine the ancestral allele and potentially swap REF/ALT.
    """
    if line.startswith("#"):
        return line

    try:
        fields = line.strip().split("\t")
        if not check_vcf_line_format(fields):
            logging.warning(f"Malformed VCF line: {line}")
            return None

        ref_allele = fields[3]
        alt_allele = fields[4]
        info_field = fields[7]
        
        # Skip complex variants - too complicated to handle properly
        if ("," in alt_allele or 
            not re.match(r"^[ACGTacgt]+$", alt_allele) or 
            not re.match(r"^[ACGTacgt]+$", ref_allele)):
            logging.debug(f"Skipping complex variant at position {fields[0]}:{fields[1]}")
            if not skip_without_consensus:
                return line
            return None

        ancestral_allele = None
        global ancestral_sample_indices
        
        # Try to get consensus from ancestral samples first
        if ancestral_sample_indices:
            ancestral_allele = determine_ancestral_consensus(fields, ancestral_sample_indices)
            if ancestral_allele:
                logging.debug(f"Got ancestral allele {ancestral_allele} from sample consensus")
                
                # Update INFO field with AA
                info_field = add_aa_to_info(info_field, ancestral_allele)
                fields[7] = info_field
        
        # Fallback to existing AA field in INFO
        if not ancestral_allele:
            aa_match = re.search(r"AA=([ACGTacgt]+)", info_field)
            if aa_match:
                ancestral_allele = aa_match.group(1).upper()
                logging.debug(f"Using existing AA={ancestral_allele} from INFO")

        # If we still don't have an ancestral allele, decide what to do
        if not ancestral_allele:
            logging.debug(f"No ancestral allele for {fields[0]}:{fields[1]}")
            if skip_without_consensus:
                return None
            return line

        # Check that ancestral allele is actually one of our alleles
        if ancestral_allele not in {ref_allele, alt_allele}:
            logging.debug(f"Ancestral allele {ancestral_allele} not in REF/ALT ({ref_allele}/{alt_allele})")
            if skip_without_consensus:
                return None
            return line

        # Decide if we need to swap based on the mode
        need_swap = False
        if mode == "ref":
            need_swap = (ancestral_allele != ref_allele)
        else:  # mode == "alt"
            need_swap = (ancestral_allele != alt_allele)

        if not need_swap:
            # No swap needed, just return with updated INFO
            return "\t".join(fields)
        
        # Time to swap REF and ALT
        updated_info = flip_allele_frequency(info_field)
        new_ref = alt_allele
        new_alt = ref_allele
        updated_fields = fields[:3] + [new_ref, new_alt] + fields[5:7] + [updated_info] + fields[8:]
        
        # Swap genotypes for all samples
        if len(updated_fields) > 9:
            updated_fields[9:] = [swap_genotype_alleles(gt) for gt in updated_fields[9:]]
            
        return "\t".join(updated_fields)
    
    except Exception as e:
        logging.error(f"Error processing line: {line[:100]}... Error: {e}")
        return None

def process_line_chunk(chunk_num: int, lines: List[str], mode: str, skip_without_consensus: bool) -> List[Optional[str]]:
    """Process a chunk of lines - used for parallel processing."""
    results = []
    for line in lines:
        result = process_vcf_line(line, mode, skip_without_consensus)
        if result is not None:
            results.append(result)
    
    return results

def extract_chromosome_data(vcf_path: str, header_lines: List[str]) -> Dict[str, List[str]]:
    """Pull out all the data organized by chromosome."""
    chromosomes = {}
    
    # Figure out how to open the file
    if vcf_path.endswith('.gz'):
        file_opener = gzip.open
        file_mode = 'rt'
    else:
        file_opener = open
        file_mode = 'r'
    
    with file_opener(vcf_path, file_mode) as f:
        # Skip the header lines
        for line in f:
            if line.startswith('#'):
                continue
            
            # Process data lines
            line_parts = line.strip().split('\t', 2)
            if len(line_parts) >= 1:
                chrom = line_parts[0]
                if chrom not in chromosomes:
                    chromosomes[chrom] = []
                chromosomes[chrom].append(line.rstrip())
    
    return chromosomes

def process_single_chromosome(chrom: str, lines: List[str], mode: str, anc_indices: List[int], 
                             threshold: float, threads_per_chrom: int, skip_without_consensus: bool) -> str:
    """Process all variants for one chromosome and save to temp file."""
    # Create a temp file for this chromosome's results
    fd, temp_file_path = tempfile.mkstemp(suffix=f"_{chrom}.vcf")
    os.close(fd)
    
    # Figure out parallelization
    num_threads = max(1, min(threads_per_chrom, len(lines) // 1000 + 1))
    
    # Set up the worker pool
    with multiprocessing.Pool(
        num_threads, 
        initializer=setup_worker_globals, 
        initargs=(anc_indices, threshold)
    ) as pool:
        # Break lines into chunks
        chunk_size = max(len(lines) // num_threads, 1000)
        line_chunks = [(i, lines[i:i + chunk_size], mode, skip_without_consensus) 
                      for i in range(0, len(lines), chunk_size)]
        
        # Process the chunks
        logging.info(f"Processing chromosome {chrom} using {num_threads} threads")
        chunk_results = list(tqdm(
            pool.starmap(process_line_chunk, line_chunks),
            total=len(line_chunks),
            unit="chunk",
            desc=f"Chr {chrom}"
        ))
    
    # Write all results to the temp file
    with open(temp_file_path, 'w') as output:
        for chunk_result in chunk_results:
            for processed_line in chunk_result:
                output.write(processed_line)
                if not processed_line.endswith('\n'):
                    output.write('\n')
    
    logging.info(f"Chromosome {chrom} done. Saved to {temp_file_path}")
    return temp_file_path

def get_sample_column_indices(header_lines: List[str], sample_names: Set[str]) -> List[int]:
    """Find which columns contain our samples of interest."""
    for line in header_lines:
        if line.startswith('#CHROM'):
            columns = line.strip().split('\t')
            if len(columns) > 9:
                # Sample columns start at index 9
                indices = []
                for i, sample in enumerate(columns[9:], 9):
                    if sample in sample_names:
                        indices.append(i)
                return indices
    return []

def read_vcf_header(vcf_path: str) -> List[str]:
    """Extract just the header lines from a VCF file."""
    header = []
    
    if vcf_path.endswith('.gz'):
        file_opener = gzip.open
        file_mode = 'rt'
    else:
        file_opener = open
        file_mode = 'r'
    
    with file_opener(vcf_path, file_mode) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('#'):
                header.append(line)
            else:
                break  # Stop at first data line
    
    return header

def make_sure_aa_header_exists(header_lines: List[str]) -> List[str]:
    """Add AA header line if it's not already there."""
    aa_info_line = '##INFO=<ID=AA,Number=1,Type=Character,Description="Ancestral allele">'
    
    # Check if we already have an AA header
    for line in header_lines:
        if line.startswith('##INFO=<ID=AA,') and 'Ancestral allele' in line:
            return header_lines  # Already there, no need to add
    
    # Find where to insert it
    chrom_line_idx = -1
    last_info_line_idx = -1
    
    for idx, line in enumerate(header_lines):
        if line.startswith('#CHROM'):
            chrom_line_idx = idx
        elif line.startswith('##INFO='):
            last_info_line_idx = idx
    
    # Insert after last INFO line, or before #CHROM line
    if last_info_line_idx >= 0:
        header_lines.insert(last_info_line_idx + 1, aa_info_line)
    elif chrom_line_idx >= 0:
        header_lines.insert(chrom_line_idx, aa_info_line)
    else:
        # Fallback - just add at the end
        header_lines.append(aa_info_line)
    
    return header_lines

def main():
    parser = argparse.ArgumentParser(
        description="VCF Recoder - polarize variants to ancestral alleles"
    )
    parser.add_argument("vcf", help="Input VCF file (can be gzipped)")
    parser.add_argument(
        "-o", "--output", 
        help="Output VCF file (default: stdout)", 
        default="stdout"
    )
    parser.add_argument(
        "--mode", 
        choices=["ref", "alt"], 
        default="ref",
        help="Whether to make REF or ALT the ancestral allele"
    )
    parser.add_argument(
        "--threads", 
        type=int, 
        default=max(1, multiprocessing.cpu_count() - 1),
        help="Number of threads to use"
    )
    parser.add_argument(
        "--log-level", 
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        default='INFO',
        help="Logging verbosity"
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
        help="Minimum proportion for ancestral consensus"
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=100000,
        help="Variants per batch when not processing by chromosome"
    )
    parser.add_argument(
        "--by-chrom",
        action="store_true",
        help="Process chromosomes in parallel (good for big files)"
    )
    parser.add_argument(
        "--tmp-dir",
        help="Where to put temporary files",
        default=None
    )
    parser.add_argument(
        "--keep-all",
        action="store_false",
        dest="skip_no_consensus",
        help="Keep variants even without ancestral consensus"
    )
    args = parser.parse_args()

    # Set up logging
    logging.getLogger().setLevel(getattr(logging, args.log_level))

    if not args.vcf:
        logging.error("Need to specify an input VCF file")
        sys.exit(1)
    
    # Set temp directory if user specified one
    if args.tmp_dir:
        if not os.path.exists(args.tmp_dir):
            os.makedirs(args.tmp_dir)
        tempfile.tempdir = args.tmp_dir

    try:
        # Open output - either stdout or a file
        if args.output == "stdout":
            output_file = sys.stdout
        elif args.output.endswith(".gz"):
            output_file = gzip.open(args.output, "wt")
        else:
            output_file = open(args.output, "w")
        
        # Load ancestral sample names if we have them
        ancestral_samples = set()
        if args.anc_list:
            with open(args.anc_list, "r") as f:
                ancestral_samples = {line.strip() for line in f if line.strip()}
            logging.info(f"Loaded {len(ancestral_samples)} ancestral sample names")
        
        # Read the VCF header
        header_lines = read_vcf_header(args.vcf)
        if not header_lines:
            logging.error("Couldn't read VCF header")
            sys.exit(1)
        
        # Make sure we have an AA header line
        header_lines = make_sure_aa_header_exists(header_lines)
        logging.info("AA header line added to output")
            
        # Find which columns have our ancestral samples
        ancestral_sample_columns = get_sample_column_indices(header_lines, ancestral_samples)
        if ancestral_samples and ancestral_sample_columns:
            logging.info(f"Found {len(ancestral_sample_columns)} ancestral samples in the VCF")
        elif ancestral_samples:
            logging.warning("None of the ancestral samples were found in this VCF")
            
        # Write the header to output
        for header_line in header_lines:
            output_file.write(header_line)
            if not header_line.endswith('\n'):
                output_file.write('\n')
        
        if args.by_chrom:
            # Chromosome-wise processing - better for really big files
            logging.info(f"Reading data by chromosome from {args.vcf}")
            chrom_data = extract_chromosome_data(args.vcf, header_lines)
            logging.info(f"Found data for {len(chrom_data)} chromosomes")
            
            # Figure out how many threads to use per chromosome
            max_parallel_chroms = min(len(chrom_data), args.threads)
            threads_per_chrom = max(1, args.threads // max_parallel_chroms)
            
            # Process chromosomes in parallel
            with multiprocessing.Pool(max_parallel_chroms) as pool:
                chrom_processing_args = [
                    (chrom, data, args.mode, ancestral_sample_columns, args.consensus_threshold, 
                     threads_per_chrom, args.skip_no_consensus)
                    for chrom, data in chrom_data.items()
                ]
                
                temp_file_paths = list(tqdm(
                    pool.starmap(process_single_chromosome, chrom_processing_args),
                    total=len(chrom_processing_args),
                    unit="chrom",
                    desc="Processing chromosomes"
                ))
            
            # Combine all the temp files into the final output
            logging.info(f"Combining results from {len(temp_file_paths)} temp files")
            for temp_path in temp_file_paths:
                with open(temp_path, 'r') as temp_file:
                    for line in temp_file:
                        output_file.write(line)
                
                # Clean up the temp file
                try:
                    os.remove(temp_path)
                except Exception as e:
                    logging.warning(f"Couldn't remove temp file {temp_path}: {e}")
        
        else:
            # Process the VCF in batches - simpler approach
            logging.info("Processing VCF in batches")
            
            def read_vcf_in_batches(vcf_path: str, batch_size: int = 100000):
                """Read VCF file in chunks to avoid memory issues."""
                current_batch = []
                
                if vcf_path.endswith('.gz'):
                    file_opener = gzip.open
                    file_mode = 'rt'
                else:
                    file_opener = open
                    file_mode = 'r'
                
                with file_opener(vcf_path, file_mode) as f:
                    for line in f:
                        line = line.rstrip()
                        if line.startswith('#'):
                            continue  # Skip header - already written
                        else:
                            current_batch.append(line)
                            if len(current_batch) >= batch_size:
                                yield current_batch
                                current_batch = []
                
                # Don't forget the last batch
                if current_batch:
                    yield current_batch
            
            # Process each batch
            batch_counter = 0
            for batch_lines in read_vcf_in_batches(args.vcf, args.batch_size):
                batch_counter += 1
                
                # Figure out threading for this batch
                threads_for_batch = max(1, min(args.threads, len(batch_lines) // 1000 + 1))
                logging.info(f"Processing batch {batch_counter} with {threads_for_batch} threads ({len(batch_lines)} variants)")
                
                chunk_size = max(len(batch_lines) // threads_for_batch, 1000)
                batch_chunks = [(i, batch_lines[i:i + chunk_size], args.mode, args.skip_no_consensus) 
                               for i in range(0, len(batch_lines), chunk_size)]
                
                # Process chunks in parallel
                with multiprocessing.Pool(
                    threads_for_batch, 
                    initializer=setup_worker_globals, 
                    initargs=(ancestral_sample_columns, args.consensus_threshold)
                ) as pool:
                    batch_results = list(tqdm(
                        pool.starmap(process_line_chunk, batch_chunks),
                        total=len(batch_chunks),
                        unit="chunk",
                        desc=f"Batch {batch_counter}"
                    ))
                
                # Write this batch's results
                for chunk_result in batch_results:
                    for processed_line in chunk_result:
                        output_file.write(processed_line)
                        if not processed_line.endswith('\n'):
                            output_file.write('\n')
                
                logging.info(f"Batch {batch_counter} complete")
        
        logging.info("All done!")
        
    except Exception as e:
        logging.critical(f"Something went wrong: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    finally:
        # Make sure we close the output file if it's not stdout
        if 'output_file' in locals() and output_file is not sys.stdout:
            output_file.close()

if __name__ == '__main__':
    main()
if __name__ == '__main__':
    main()
