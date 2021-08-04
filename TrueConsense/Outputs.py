import os
import sys

from .Coverage import GetCoverage
from .Events import ListInserts
from .indexing import Readbam
from .Sequences import BuildConsensus

from datetime import date
from Bio import SeqIO


def WriteGFF(gffheader, gffdict, outdir, name, cov):
    with open(f"{outdir}/{name}_cov_ge_{cov}.gff", "w") as out:
        out.write(gffheader)

        for k, v in gffdict.items():
            for nk, nv in v.items():
                if str(nk) == str(list(v.keys())[-1]):
                    out.write(str(nv))
                else:
                    out.write(str(nv) + "\t")
            out.write("\n")
    pass


def WriteOutputs(
    cov,
    iDict,
    uGffDict,
    inputbam,
    IncludeAmbig,
    WriteVCF,
    name,
    ref,
    gffout,
    gffheader,
    outdir,
):
    """
    step 1: construct the consensus sequences, both with and without inserts
    step 2: write the vcf file
    step 3: write the consensus sequence
    """
    today = date.today().strftime("%Y%m%d")

    bam = Readbam(inputbam)
    consensus, newgff = BuildConsensus(cov, iDict, uGffDict, IncludeAmbig, bam, True)
    consensus_noinsert = BuildConsensus(cov, iDict, uGffDict, IncludeAmbig, bam, False)[
        0
    ]

    if gffout is not None:
        WriteGFF(gffheader, newgff, gffout, name, cov)

    if WriteVCF is not None:
        hasinserts, insertpositions = ListInserts(iDict, cov, bam)

        q = 0
        for record in SeqIO.parse(ref, "fasta"):
            if q != 0:
                break
            q += 1
            refID = record.id
            reflist = list(record.seq)

        seqlist = list(consensus_noinsert.upper())

        with open(f"{os.path.abspath(WriteVCF)}/{name}_cov_ge_{cov}.vcf", "w") as out:
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
                            if currentcov > cov:
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

    with open(f"{os.path.abspath(outdir)}/{name}_cov_ge_{cov}.fa", "w") as out:
        out.write(f">{name}_cov_ge_{cov}\n{consensus}\n")

    pass
