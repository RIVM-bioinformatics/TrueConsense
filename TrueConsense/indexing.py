import gffpandas.gffpandas as gffpd
import pandas as pd
import pysamstats
import pysam


def Readbam(f):
    return pysam.AlignmentFile(f, "rb")


def Gffindex(file):
    return gffpd.read_gff3(file)


def BuildIndex(bamfile, ref):
    columns = ["coverage", "A", "T", "C", "G", "X", "I"]
    p_index = pd.DataFrame(columns=columns)

    for r in pysamstats.stat_pileup(
        type="variation",
        alignmentfile=Readbam(bamfile),
        stepper="nofilter",
        fafile=ref,
        pad=True,
        one_based=True,
        max_depth=1000000000,
    ):
        p_index.loc[r["pos"]] = (
            [r["reads_all"]]
            + [r["A"]]
            + [r["T"]]
            + [r["C"]]
            + [r["G"]]
            + [r["deletions"]]
            + [r["insertions"]]
        )

    return p_index
