"""
>add some notices here
"""

import argparse
import concurrent.futures as cf
import multiprocessing
import os
import pathlib
from re import M
import sys

from .func import MyHelpFormatter, color
from .indexing import (
    Gffindex,
    ReadBam,
    read_override_positions,
    ReadFasta,
)
from .Outputs import WriteOutputs
from .version import __version__
from .Calls import Calls


def GetArgs(givenargs):
    def checkbam(fname):
        if os.path.isfile(fname):
            ext = "".join(pathlib.Path(fname).suffix)
            if ext != ".bam":
                parser.error(
                    f"Input file {color.YELLOW}({fname}){color.END} doesn't seem to be a BAM-file."
                )
            return fname
        print(f'"{fname}" is not a file. Exiting...')
        sys.exit(-1)

    def checkfasta(fname):
        allowedexts = [".fasta", ".fa"]
        if os.path.isfile(fname):
            ext = "".join(pathlib.Path(fname).suffix)
            if ext not in allowedexts:
                parser.error(
                    f"Reference file {color.YELLOW}({fname}){color.END} doesn't seem to be a Fasta-file."
                )
            return fname
        print(f'"{fname}" is not a file. Exiting...')
        sys.exit(1)

    def checkgff(fname):
        if os.path.isfile(fname):
            ext = "".join(pathlib.Path(fname).suffix)
            if ext != ".gff":
                parser.error(
                    f"Given file {color.YELLOW}({fname}){color.END} doesn't seem to be a GFF file."
                )
            return fname
        print(f'"{fname}" is not a file. Exiting...')
        sys.exit(1)

    def check_index_override(filenames):
        csv_fname, bam_fname = filenames
        bam_fname = checkbam(bam_fname)
        if os.path.isfile(csv_fname):
            ext = "".join(pathlib.Path(csv_fname).suffixes)
            if ".csv" not in ext:
                parser.error(
                    f"Given file {color.YELLOW}({csv_fname}){color.END} doesn't seem to be a compressed csv file."
                )
            if ".gz" not in ext:
                parser.error(
                    f"Given file {color.YELLOW}({csv_fname}){color.END} doesn't seem to be a compressed csv file."
                )
            return csv_fname, bam_fname
        print(f'"{csv_fname}" is not a file. Exiting...')
        sys.exit(1)

    parser = argparse.ArgumentParser(
        prog="TrueConsense",
        usage="%(prog)s [required options] [optional arguments]",
        description="TrueConsense: Creating biologically valid consensus sequences from reference-based alignments",
        formatter_class=MyHelpFormatter,
        add_help=False,
    )

    standard_threads = min(multiprocessing.cpu_count(), 128)

    reqs = parser.add_argument_group("Required arguments")

    reqs.add_argument(
        "--input",
        "-i",
        type=checkbam,
        metavar="File",
        help="Input file in BAM format",
        required=True,
    )

    reqs.add_argument(
        "--output",
        "-o",
        type=str,
        default=os.getcwd() + "consensus.fasta",
        metavar="File",
        help="Output consensus fasta",
        required=True,
    )

    reqs.add_argument(
        "--reference",
        "-ref",
        type=checkfasta,
        metavar="File",
        help="Reference Fasta file",
        required=True,
    )

    reqs.add_argument(
        "--features",
        "-gff",
        type=checkgff,
        metavar="File",
        help="File with genome features (GFF)",
        required=True,
    )

    reqs.add_argument(
        "--coverage-level",
        "-cov",
        type=int,
        default=30,
        metavar="100",
        help="The minimum coverage level of the consensus and variant calls",
        required=True,
    )

    reqs.add_argument(
        "--samplename",
        "-name",
        metavar="Text",
        help="Name of the sample that is being processed, will be used to create the fasta header",
        required=True,
    )

    opts = parser.add_argument_group("Optional arguments")

    opts.add_argument(
        "--variants",
        "-vcf",
        type=str,
        metavar="File",
        help="Output VCF file",
    )

    opts.add_argument(
        "--depth-of-coverage",
        "-doc",
        type=str,
        metavar="File",
        help="Output TSV file listing the coverage per position",
    )

    opts.add_argument(
        "--output-gff",
        "-ogff",
        type=str,
        metavar="File",
        help="Ouput location a corrected GFF file",
    )

    opts.add_argument(
        "--threads",
        "-t",
        default=standard_threads,
        metavar="N",
        help="Number of threads that can be used by TrueConsense",
        type=int,
    )

    opts.add_argument(
        "--noambiguity",
        "-noambig",
        action="store_true",
        help="Turn off ambiguity nucleotides in the generated consensus sequence",
    )

    opts.add_argument(
        "--index-override",
        nargs=2,
        metavar="FILE",
        help="Override the positional index of certain genome positions with 'known' information if the given alignment is not sufficient for these positions.\nFirst FILE must be a compressed csv containing the positions in the first column. The second FILE must be a bamfile to read these positions from\nPlease use with caution as this will overwrite the generated index at the given positions!\n",
    )

    opts.add_argument(
        "--version",
        "-v",
        action="version",
        version=__version__,
        help="Show the TrueConsense version and exit",
    )

    opts.add_argument(
        "--help",
        "-h",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit",
    )

    args = parser.parse_args(givenargs)
    if args.index_override:
        check_index_override(args.index_override)

    return args


def main():
    if len(sys.argv[1:]) < 1:
        print(
            "TrueConsense was called but no arguments were given, please try again.\nUse 'TrueConsense -h' to see the help document"
        )
        sys.exit(1)
    args = GetArgs(sys.argv[1:])

    with cf.ThreadPoolExecutor(max_workers=args.threads) as xc:
        IndexGff = xc.submit(Gffindex, args.features)
        ref_obj = xc.submit(ReadFasta, args.reference)
        # IndexDF = xc.submit(BuildIndex, args.input, args.reference)
        gff_obj = IndexGff.result()
        ref_obj = ref_obj.result()
        # IndexDF = IndexDF.result()

    ref_file_location = args.reference
    ref_contig = ref_obj.references[0]
    ref_seq = ref_obj.fetch(ref_contig)

    IncludeAmbig = not args.noambiguity

    call_obj = Calls(
        ref_file_location,
        ref_contig,
        ref_seq,
        IncludeAmbig=IncludeAmbig,
        significance=0.5,
    )
    call_obj.fill_positions_from_bam(bamfile=args.input)
    if args.index_override:
        override_positions_file, override_bam_file = args.index_override
        positions = list(read_override_positions(override_positions_file).index)
        call_obj.fill_positions_from_bam(override_bam_file, positions)
    call_obj.calculate_scores()

    # with cf.ThreadPoolExecutor(max_workers=args.threads) as xc:
    if args.depth_of_coverage is not None:
        call_obj.p_index[["cov"]].to_csv(args.depth_of_coverage, sep="\t", header=False)

    WriteOutputs(
        args.coverage_level,
        call_obj,
        gff_obj,  # This is a dataframe of the gff
        args.variants,
        args.samplename,
        args.output_gff,
        args.output,
    )
