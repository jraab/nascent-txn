#!/usr/bin/env python

"""
Validate and reformat TT-seq sample sheet for Nextflow processing.

This script checks that the sample sheet has the required columns and file paths exist.
"""

import argparse
import csv
import logging
import sys
from pathlib import Path


def parse_args(args=None):
    Description = "Reformat TT-seq samplesheet file and check its contents."
    Epilog = """Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"""
    
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        Path(path).mkdir(parents=True, exist_ok=True)


def print_error(error, context="Line", context_str=""):
    error_str = f"ERROR: Please check samplesheet -> {error}"
    if context != "" and context_str != "":
        error_str = f"ERROR: Please check samplesheet -> {error}\\n{context.strip()}: '{context_str.strip()}'"
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:
    
    sample,fastq_1,fastq_2,single_end,strandedness
    SAMPLE_PE,SAMPLE_PE_RUN1_1.fastq.gz,SAMPLE_PE_RUN1_2.fastq.gz,FALSE,forward
    SAMPLE_SE,SAMPLE_SE_RUN1_1.fastq.gz,,TRUE,forward
    """
    
    sample_mapping_dict = {}
    with open(file_in, "r") as fin:
        
        ## Check header
        MIN_COLS = 2
        HEADER = ["sample", "fastq_1", "fastq_2"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        
        ## Find required columns
        if "sample" not in header:
            print_error("'sample' column not found!")
        if "fastq_1" not in header:
            print_error("'fastq_1' column not found!")
            
        sample_col = header.index("sample")
        fastq1_col = header.index("fastq_1")
        fastq2_col = header.index("fastq_2") if "fastq_2" in header else None
        single_end_col = header.index("single_end") if "single_end" in header else None
        strandedness_col = header.index("strandedness") if "strandedness" in header else None
        
        ## Check sample entries
        for line_idx, line in enumerate(fin, 2):
            if line.strip():
                lspl = [x.strip().strip('"') for x in line.strip().split(",")]
                
                # Check valid number of columns per row
                if len(lspl) < len(header):
                    print_error(
                        "Invalid number of columns (minimum = {})!".format(len(header)),
                        "Line",
                        line,
                    )
                num_cols = len([x for x in lspl if x])
                if num_cols < MIN_COLS:
                    print_error(
                        "Invalid number of populated columns (minimum = {})!".format(MIN_COLS),
                        "Line", 
                        line,
                    )
                
                ## Check sample name entries
                sample, fastq_1 = lspl[sample_col], lspl[fastq1_col]
                if sample.find(" ") != -1:
                    print_error("Sample entry contains spaces!", "Line", line)
                if not sample:
                    print_error("Sample entry has not been specified!", "Line", line)
                
                ## Check FastQ file extension
                for fastq in [fastq_1]:
                    if fastq:
                        if fastq.find(" ") != -1:
                            print_error("FastQ file contains spaces!", "Line", line)
                        if not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz"):
                            print_error(
                                "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
                                "Line",
                                line,
                            )
                
                ## Check FastQ_2 if paired-end
                fastq_2 = ""
                single_end = True
                if fastq2_col and lspl[fastq2_col]:
                    fastq_2 = lspl[fastq2_col]
                    single_end = False
                    if fastq_2.find(" ") != -1:
                        print_error("FastQ file contains spaces!", "Line", line)
                    if not fastq_2.endswith(".fastq.gz") and not fastq_2.endswith(".fq.gz"):
                        print_error(
                            "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
                            "Line",
                            line,
                        )
                
                ## Auto-detect paired-end/single-end
                if single_end_col:
                    single_end = lspl[single_end_col].upper() == "TRUE"
                else:
                    single_end = not bool(fastq_2)
                
                ## Check strandedness
                strandedness = "unstranded"
                if strandedness_col and lspl[strandedness_col]:
                    strandedness = lspl[strandedness_col]
                    if strandedness not in ["unstranded", "forward", "reverse"]:
                        print_error(
                            f"Strandedness must be 'unstranded', 'forward', or 'reverse' (got '{strandedness}')!",
                            "Line",
                            line
                        )
                
                ## Create sample mapping dictionary with metadata
                sample_info = [sample, fastq_1, fastq_2, single_end, strandedness]
                if sample not in sample_mapping_dict:
                    sample_mapping_dict[sample] = [sample_info]
                else:
                    if sample_info in sample_mapping_dict[sample]:
                        print_error("Samplesheet contains duplicate rows!", "Line", line)
                    else:
                        sample_mapping_dict[sample].append(sample_info)

    ## Write validated samplesheet with appropriate header.
    if len(sample_mapping_dict) > 0:
        out_dir = Path(file_out).parent
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(",".join(["sample", "fastq_1", "fastq_2", "single_end", "strandedness"]) + "\\n")
            for sample in sorted(sample_mapping_dict.keys()):
                ## Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
                if not all(x[3] == sample_mapping_dict[sample][0][3] for x in sample_mapping_dict[sample]):
                    print_error(f"Multiple runs of a sample must be of the same datatype!", "Sample", sample)

                for idx, sample_info in enumerate(sample_mapping_dict[sample]):
                    sample, fastq_1, fastq_2, single_end, strandedness = sample_info
                    if idx == 0:
                        fout.write(",".join([sample, fastq_1, fastq_2, str(single_end).upper(), strandedness]) + "\\n")
                    else:
                        fout.write(",".join([f"{sample}_T{idx+1}", fastq_1, fastq_2, str(single_end).upper(), strandedness]) + "\\n")
    else:
        print_error("No entries to process!", "Samplesheet", file_in)


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    main()