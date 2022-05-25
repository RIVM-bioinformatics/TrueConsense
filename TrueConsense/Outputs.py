import os
import sys
from datetime import date

from .ORFs import RestoreORFS


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
    call_obj,
    gff_obj,
    WriteVCF,
    name,
    gffout,
    outdir,
):
    """
    step 1: construct the consensus sequences, both with and without inserts
    step 2: write the vcf file
    step 3: write the consensus sequence
    """
    today = date.today().strftime("%Y%m%d")

    call_obj.remove_calls_with_low_coverage(mincov)
    call_obj.sort_highest_score()
    call_obj.pick_first_in_calls()
    RestoreORFS(call_obj, gff_obj.attributes_to_columns())

    print("Making consensus")
    consensus = call_obj.consensus()
    print("Done making consensus")

    if gffout is not None:
        call_obj.update_gff_coods_with_insertions(gff_obj)
        gff_obj.to_gff3(f"{gffout}/{name}_cov_ge_{mincov}.gff")

    if WriteVCF is not None:
        with open(
            f"{os.path.abspath(WriteVCF)}/{name}_cov_ge_{mincov}.vcf", "w"
        ) as out:
            out.write(
                f"""##fileformat=VCFv4.3
##fileDate={today}
##source='TrueConsense {' '.join(sys.argv[1:])}'
##reference='{call_obj.ref_file_location}'
##contig=<ID={call_obj.ref_contig}>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""
            )
            # writecontents
            for (
                pos,
                refnuc,
                altnuc,
                cov,
            ) in call_obj.get_mutations():
                out.write(
                    f"{call_obj.ref_contig}\t{pos}\t.\t{refnuc}\t{altnuc}\t.\tPASS\tDP={cov}"
                )
                if len(refnuc) != len(altnuc):
                    out.write(";INDEL\n")
                else:
                    out.write("\n")

    with open(f"{os.path.abspath(outdir)}/{name}_cov_ge_{mincov}.fa", "w") as out:
        out.write(f">{name}_cov_ge_{mincov}\n{consensus}\n")
