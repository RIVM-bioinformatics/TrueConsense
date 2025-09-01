import os
import sys
from datetime import date

from Bio import SeqIO

from .Coverage import GetCoverage
from .Events import ListInserts
from .indexing import Readbam
from .Sequences import BuildConsensus
from AminoExtract.gff_data import GFFColumns

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
    cols = GFFColumns.get_names()
    cols_without_attr = [col for col in cols if col != "attributes"]

    def combine_dict_into_attributes(input_dict: dict[str, str]) -> str:
        attribute_dict = {}
        for k, v in input_dict.items():
            if k == "attributes":
                for attribute in v.split(";"):
                    if attribute == "":
                        continue
                    key, value = attribute.split("=")
                    attribute_dict[key] = value
            else:
                attribute_dict[k] = v

        return ";".join(f"{k}={v}" for k, v in attribute_dict.items())

    def clean_dict(input_dict: dict[str, str]) -> dict[str, str]:
        clean_dict = {}
        attribute_dict = {}
        for k, v in input_dict.items():
            if str(k).lower() not in cols_without_attr:
                attribute_dict[str(k).lower()] = str(v)
            else:
                clean_dict[str(k).lower()] = str(v)

        attribute_str = combine_dict_into_attributes(attribute_dict)
        clean_dict["attributes"] = attribute_str
        assert list(clean_dict.keys()) == cols
        return clean_dict


    # the gffdict will have 0, 1, 2, etc as keys, for each line in the GFF file
    # the values will be dictionaries containing the GFF columns for that line, with a lot of additional columns
    # these additional columns will all be forced into the attributes column

    with open(output_gff, "w") as out:
        out.write(gffheader.raw_text)


        for line_number, gff_data in gffdict.items():
            cleaned_data = clean_dict(gff_data)
            out.write("\t".join([str(v) for v in cleaned_data.values()]) + "\n")


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
