import os
import sys
from datetime import date

from Bio import SeqIO

from .Coverage import GetCoverage
from .Events import ListInserts
from .indexing import Readbam
from .Sequences import BuildConsensus


def WriteGFF(gffheader, gffdict, output_gff, name):
    """Function takes a GFF header, a dictionary of GFF features, an output directory, and a name for
    the output file, and writes the GFF header and the GFF features to a file in the output directory

    Parameters
    ----------
    gffheader
        The header of the GFF file.
    gffdict
        a dictionary of GFF features
    outdir
        the directory where you want the output files to go
    name
        the name of the file you want to write

    """
    with open(output_gff, "w") as out:
        out.write(gffheader)

        for k, v in gffdict.items():
            for nk, nv in v.items():
                if str(nk) == str(list(v.keys())[-1]):
                    out.write(str(nv))
                else:
                    out.write(str(nv) + "\t")
            out.write("\n")


def WriteOutputs(
    mincov,
    iDict,
    uGffDict,
    inputbam,
    IncludeAmbig,
    output_vcf,
    name,
    ref,
    output_gff,
    gffheader,
    output_consensus,
):
    """
    step 1: construct the consensus sequences, both with and without inserts
    step 2: write the vcf file
    step 3: write the consensus sequence
    """
    today = date.today().strftime("%Y%m%d")

    bam = Readbam(inputbam)
    consensus, newgff = BuildConsensus(mincov, iDict, uGffDict, IncludeAmbig, bam, True)
    consensus_noinsert = BuildConsensus(
        mincov, iDict, uGffDict, IncludeAmbig, bam, False
    )[0]

    if output_gff is not None:
        WriteGFF(gffheader, newgff, output_gff, name)

    if output_vcf is not None:
        hasinserts, insertpositions = ListInserts(iDict, mincov, bam)

        q = 0
        for record in SeqIO.parse(ref, "fasta"):
            if q != 0:
                break
            q += 1
            refID = record.id
            reflist = list(record.seq)

        seqlist = list(consensus_noinsert.upper())

        with open(output_vcf, "w") as out:
            out.write(
                f"""##fileformat=VCFv4.3
##fileDate={today}
##source='TrueConsense {' '.join(sys.argv[1:])}'
##reference='{ref}'
##contig=<ID={refID}>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""
            )
            # writecontents

            delskips = []
            for i in range(len(reflist)):
                if i in delskips:
                    continue

                if reflist[i] != seqlist[i]:
                    b = i

                    if seqlist[i] == "-":

                        gapextendedreflist = []
                        while seqlist[b] == "-":
                            gapextendedreflist.append(reflist[b])
                            delskips.append(b)
                            b += 1

                        gapextension = "".join(gapextendedreflist)
                        joinedreflist = str(reflist[i - 1] + gapextension)

                        currentcov = GetCoverage(iDict, i + 1)
                        out.write(
                            f"{refID}\t{i}\t.\t{joinedreflist}\t{seqlist[i-1]}\t.\tPASS\tDP={currentcov};INDEL\n"
                        )
                    else:
                        if i == 0:
                            p = 1
                        elif i == 1:
                            p = 1
                        else:
                            p = i
                        currentcov = GetCoverage(iDict, p + 1)
                        out.write(
                            f"{refID}\t{i+1}\t.\t{reflist[i]}\t{seqlist[i]}\t.\tPASS\tDP={currentcov}\n"
                        )
                if hasinserts is True:
                    for lposition in insertpositions:
                        if i == lposition:
                            currentcov = GetCoverage(iDict, i + 1)
                            if currentcov > mincov:
                                for y in insertpositions.get(lposition):
                                    to_insert = str(
                                        insertpositions.get(lposition).get(y)
                                    )
                                    if to_insert is not None:
                                        CombinedEntry = seqlist[i] + to_insert
                                        out.write(
                                            f"{refID}\t{i}\t.\t{reflist[i]}\t{CombinedEntry}\t.\tPASS\tDP={currentcov};INDEL\n"
                                        )
                                    else:
                                        continue

    with open(output_consensus, "w") as out:
        out.write(f">{name} mincov={mincov}\n{consensus}\n")
