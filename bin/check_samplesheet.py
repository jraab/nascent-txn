#!/usr/bin/env python3

import os
import sys
import errno
import argparse


def parse_args(args=None):
    Description = "Reformat TT-seq samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context, context_str.strip()
        )
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:

    sample,fastq_1,fastq_2,single_end,strandedness
    SAMPLE_1,SAMPLE_1_R1.fastq.gz,SAMPLE_1_R2.fastq.gz,FALSE,forward
    SAMPLE_2,SAMPLE_2_R1.fastq.gz,,TRUE,forward
    """

    sample_mapping_dict = {}
    with open(file_in, "r") as fin:

        ## Check header - first check for the expected TT-seq format
        header_line = fin.readline()
        
        # Clean and parse header
        header = [x.strip().strip('"') for x in header_line.strip().split(",")]
        
        # Remove BOM from first column if present
        if header and header[0].startswith('\ufeff'):
            header[0] = header[0][1:]  # Remove the BOM character
            
        # Check for TT-seq format first: sample,fastq_1,fastq_2,single_end,strandedness
        TTSEQ_HEADER = ["sample", "fastq_1", "fastq_2", "single_end", "strandedness"]
        # Check for diffpolyA format: sample,fastq_1,fastq_2,group
        DIFFPOLYA_HEADER = ["sample", "fastq_1", "fastq_2", "group"]
        
        if len(header) >= len(TTSEQ_HEADER) and header[:len(TTSEQ_HEADER)] == TTSEQ_HEADER:
            # TT-seq format
            use_ttseq_format = True
            MIN_COLS = 3  # sample, fastq_1, single_end minimum
        elif len(header) >= len(DIFFPOLYA_HEADER) and header[:len(DIFFPOLYA_HEADER)] == DIFFPOLYA_HEADER:
            # diffpolyA format - convert to TT-seq format
            use_ttseq_format = False
            MIN_COLS = 3  # sample, fastq_1, group minimum
        else:
            print("ERROR: Please check samplesheet header -> Expected either:")
            print("  TT-seq format: {}".format(",".join(TTSEQ_HEADER)))
            print("  or diffpolyA format: {}".format(",".join(DIFFPOLYA_HEADER)))
            print("  Got: {}".format(",".join(header)))
            sys.exit(1)

        ## Check sample entries
        for line_idx, line in enumerate(fin):
            if line.strip():
                lspl = [x.strip().strip('"') for x in line.strip().split(",")]

                # Check valid number of columns per row
                min_expected = len(TTSEQ_HEADER) if use_ttseq_format else len(DIFFPOLYA_HEADER)
                if len(lspl) < min_expected:
                    print_error(
                        "Invalid number of columns (minimum = {})!".format(min_expected),
                        "Line",
                        line,
                    )

                if use_ttseq_format:
                    ## TT-seq format processing
                    sample, fastq_1, fastq_2, single_end_str, strandedness = lspl[:5]
                    
                    # Validate single_end field
                    if single_end_str.upper() not in ["TRUE", "FALSE"]:
                        print_error("single_end must be TRUE or FALSE!", "Line", line)
                    single_end = single_end_str.upper() == "TRUE"
                    
                    # Validate strandedness
                    if strandedness not in ["unstranded", "forward", "reverse"]:
                        print_error("strandedness must be 'unstranded', 'forward', or 'reverse'!", "Line", line)
                        
                else:
                    ## diffpolyA format - convert to TT-seq format
                    sample, fastq_1, fastq_2, group = lspl[:4]
                    # Auto-detect single_end and set default strandedness
                    single_end = not fastq_2.strip()
                    strandedness = "unstranded"  # Default for converted files

                ## Common validation
                sample = sample.replace(" ", "_")
                if not sample:
                    print_error("Sample entry has not been specified!", "Line", line)

                ## Check FastQ file extension
                for fastq in [fastq_1, fastq_2]:
                    if fastq:
                        if fastq.find(" ") != -1:
                            print_error("FastQ file contains spaces!", "Line", line)
                        if not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz"):
                            print_error(
                                "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
                                "Line",
                                line,
                            )

                ## Auto-detect paired-end/single-end if not explicit
                sample_info = []  ## [single_end, fastq_1, fastq_2, strandedness]
                if sample and fastq_1 and fastq_2.strip() and not single_end:  ## Paired-end
                    sample_info = [False, fastq_1, fastq_2, strandedness]
                elif sample and fastq_1 and (not fastq_2.strip() or single_end):  ## Single-end
                    sample_info = [True, fastq_1, "", strandedness]
                else:
                    print_error("Invalid combination of columns provided!", "Line", line)

                ## Create sample mapping dictionary = { sample: [ single_end, fastq_1, fastq_2, strandedness ] }
                if sample not in sample_mapping_dict:
                    sample_mapping_dict[sample] = [sample_info]
                else:
                    if sample_info in sample_mapping_dict[sample]:
                        print_error("Samplesheet contains duplicate rows!", "Line", line)
                    else:
                        sample_mapping_dict[sample].append(sample_info)

    ## Write validated samplesheet with TT-seq format
    if len(sample_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(",".join(["sample", "single_end", "fastq_1", "fastq_2", "strandedness"]) + "\n")
            for sample in sorted(sample_mapping_dict.keys()):

                ## Check that multiple runs of the same sample are of the same datatype
                if not all(x[0] == sample_mapping_dict[sample][0][0] for x in sample_mapping_dict[sample]):
                    print_error("Multiple runs of a sample must be of the same datatype!", "Sample: {}".format(sample))

                for idx, val in enumerate(sample_mapping_dict[sample]):
                    single_end, fastq_1, fastq_2, strandedness = val
                    # Use original sample name if only one entry, otherwise add suffix
                    sample_name = sample if len(sample_mapping_dict[sample]) == 1 else "{}_T{}".format(sample, idx+1)
                    single_end_str = "TRUE" if single_end else "FALSE"
                    fout.write(",".join([sample_name, single_end_str, fastq_1, fastq_2, strandedness]) + "\n")
    else:
        print_error("No entries to process!", "Samplesheet: {}".format(file_in))


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    main()