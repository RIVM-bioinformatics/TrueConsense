"""
>add some notices here
"""

import os
from os import path
import pathlib
import sys
import argparse
import multiprocessing

import concurrent.futures as cf

from .version import __version__
from .func import MyHelpFormatter, color


def GetArgs(givenargs):
    
    
    def checkbam(fname):
        if os.path.isfile(fname):
            ext = "".join(pathlib.Path(fname).suffix)
            if ext != '.bam':
                parser.error(f"Input file {color.YELLOW}({fname}){color.END} doesn't seem to be a BAM-file.")
            return fname
        else:
            print(f'"{fname}" is not a file. Exiting...')
            sys.exit(-1)
            
    def checkfasta(fname):
        allowedexts = ['.fasta', '.fa']
        if os.path.isfile(fname):
            ext = "".join(pathlib.Path(fname).suffix)
            if ext not in allowedexts:
                parser.error(f"Reference file {color.YELLOW}({fname}){color.END} doesn't seem to be a Fasta-file.")
            return fname
        else:
            print(f'"{fname}" is not a file. Exiting...')
            sys.exit(1)
            
    def checkgff(fname):
        if os.path.isfile(fname):
            ext = "".join(pathlib.Path(fname).suffix)
            if ext != '.gff':
                parser.error(f"Given file {color.YELLOW}({fname}){color.END} doesn't seem to be a GFF file.")
            return fname
        else:
            print(f'"{fname}" is not a file. Exiting...')
            sys.exit(1)
            
    def currentpath():
        return os.getcwd()
    
    parser = argparse.ArgumentParser(
        prog='TrueConsense',
        usage="%(prog)s [required options] [optional arguments]",
        description="TrueConsense: Creating biologically correct consensus sequences from reference-based alignments",
        formatter_class=MyHelpFormatter,
        add_help=False
    )
    
    standard_threads = min(multiprocessing.cpu_count(), 128)
    
    reqs = parser.add_argument_group('Required arguments')
    
    reqs.add_argument(
        '--input', '-i',
        type=lambda s: checkbam(s),
        metavar="File",
        help="Input file in BAM format",
        required=True
    )

    reqs.add_argument(
        '--ouput', '-o',
        type=str,
        default=currentpath(),
        metavar='DIR',
        help='Output directory where the (various) consensus fasta files will be placed',
        required=True
    )

    reqs.add_argument(
        '--reference', '-ref',
        type=lambda s: checkfasta(s),
        metavar="File",
        help="Reference Fasta file",
        required=True
    )

    reqs.add_argument(
        '--features', '-gff',
        type=lambda s: checkgff(s),
        metavar="File",
        help="File with genome features (GFF)",
        required=True
    )
    
    reqs.add_argument(
        '--coverage-levels', '-cov',
        type=int,
        nargs='+',
        default=1,
        metavar="1 5 10",
        help="List of one or multiple 'levels' that will be used as the minimum coverage",
        required=True
    )
    
    reqs.add_argument(
        '--samplename', '-name',
        metavar='Text',
        help='Name of the sample that is being processed, will be used to create the fasta files and fasta header(s)',
        required=True
    )

    opts = parser.add_argument_group('Optional arguments')
    
    opts.add_argument(
        '--variants', '-vcf',
        type=str,
        metavar='DIR',
        help='Create VCF files for every given coverage level within the given directory. The created files start with the given samplename'
    )
    
    opts.add_argument(
        '--depth-of-coverage', '-doc',
        type=str,
        metavar='File',
        help='Output TSV file listing the coverage per position'
    )
    
    opts.add_argument(
        '--threads', '-t',
        default=standard_threads,
        metavar='N',
        help='Number of threads that can be used by TrueConsense.',
        type=int,
        
    )
    
    opts.add_argument(
        '--help', '-h',
        action='help',
        default=argparse.SUPPRESS,
        help="Show this help message and exit"
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
    print(args)