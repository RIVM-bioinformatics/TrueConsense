"""
>add some notices here
"""

import argparse
import concurrent.futures as cf
import multiprocessing
import os
import pathlib
import sys

import parmap

from .Coverage import BuildCoverage
from .func import MyHelpFormatter, color
from .indexing import (
    BuildIndex,
    Gffindex,
    Override_index_positions,
    Readbam,
    read_override_index,
)
from .Outputs import WriteOutputs
from .version import __version__


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

    def check_index_override(fname):
        if os.path.isfile(fname):
            ext = "".join(pathlib.Path(fname).suffixes)
            if ".csv" not in ext:
                parser.error(
                    f"Given file {color.YELLOW}({fname}){color.END} doesn't seem to be a compressed csv file."
                )
            if ".gz" not in ext:
                parser.error(
                    f"Given file {color.YELLOW}({fname}){color.END} doesn't seem to be a compressed csv file."
                )
            return fname
        print(f'"{fname}" is not a file. Exiting...')
        sys.exit(1)

    def currentpath():
        return os.getcwd()

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
        default=currentpath(),
        metavar="DIR",
        help="Output directory where the (various) consensus fasta files will be placed",
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
        "--coverage-levels",
        "-cov",
        type=int,
        nargs="+",
        default=1,
        metavar="1 5 10",
        help="List of one or multiple 'levels' that will be used as the minimum coverage",
        required=True,
    )

    reqs.add_argument(
        "--samplename",
        "-name",
        metavar="Text",
        help="Name of the sample that is being processed, will be used to create the fasta files and fasta header(s)",
        required=True,
    )

    opts = parser.add_argument_group("Optional arguments")

    opts.add_argument(
        "--variants",
        "-vcf",
        type=str,
        metavar="DIR",
        help="Create VCF files for every given coverage level within the given directory. The created files start with the given samplename",
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
        metavar="DIR",
        help="Create a corrected GFF file for every given coverage level within the given directory. The created files start with the given samplename",
    )

    opts.add_argument(
        "--threads",
        "-t",
        default=standard_threads,
        metavar="N",
        help="Number of threads that can be used by TrueConsense.",
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
        type=check_index_override,
        metavar="File",
        help="Override the positional index of certain genome positions with 'known' information if the given alignment is not sufficient for these positions\nMust be a compressed csv.\nPlease use with caution as this will overwrite the generated index at the given positions!\n",
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

    return args


def main():
    if len(sys.argv[1:]) < 1:
        print(
            "TrueConsense was called but no arguments were given, please try again.\nUse 'TrueConsense -h' to see the help document"
        )
        sys.exit(1)
    args = GetArgs(sys.argv[1:])

    bam = Readbam(args.input)

    with cf.ThreadPoolExecutor(max_workers=args.threads) as xc:
        IndexDF = xc.submit(BuildIndex, args.input, args.reference)
        IndexGff = xc.submit(Gffindex, args.features)

        IndexDF = IndexDF.result()
        IndexGff = IndexGff.result()

    if args.index_override:
        IndexDF = Override_index_positions(
            IndexDF, read_override_index(args.index_override)
        )

    indexDict = IndexDF.to_dict("index")
    GffHeader = IndexGff.header
    GffDF = IndexGff.df
    GffDict = GffDF.to_dict("index")

    with cf.ThreadPoolExecutor(max_workers=args.threads) as xc:
        if args.depth_of_coverage is not None:
            xc.submit(BuildCoverage, indexDict, args.depth_of_coverage)

    if args.noambiguity is False:
        IncludeAmbig = True
    elif args.noambiguity is True:
        IncludeAmbig = False

    parallel(
        WriteOutputs,
        args.coverage_levels,
        indexDict,
        GffDict,
        args.input,
        IncludeAmbig,
        args.variants,
        args.samplename,
        args.reference,
        args.output_gff,
        GffHeader,
        args.output,
        args.threads,
    )


def parallel(
    function,
    covlist,
    indexDict,
    GffDict,
    inputbam,
    IncludeAmbig,
    WriteVCF,
    name,
    ref,
    gffoutdir,
    gffheader,
    outdir,
    workers,
):
    parmap.map(
        function,
        covlist,
        indexDict,
        GffDict,
        inputbam,
        IncludeAmbig,
        WriteVCF,
        name,
        ref,
        gffoutdir,
        gffheader,
        outdir,
        pm_processes=workers,
    )
