import gffpandas.gffpandas as gffpd
import pandas as pd
import pysam


def Readbam(f):
    return pysam.AlignmentFile(f, "rb")


def Gffindex(file):
    return gffpd.read_gff3(file)


def read_override_index(f):
    return pd.read_csv(f, sep=",", compression="gzip", index_col=0)


def Override_index_positions(index, override_data):
    index.loc[override_data.index, :] = override_data[:]
    return index


def BuildIndex(bamfile, ref):
    bamfile = pysam.AlignmentFile(bamfile, "rb")
    ref_fasta = pysam.FastaFile(ref)
    ref_length = ref_fasta.lengths[0]

    pileup = bamfile.pileup(stepper="nofilter", max_depth=10000000, min_base_quality=0)

    def parse_query_sequences(l):
        coverage = a = c = t = g = x = i = 0
        for b in l:
            coverage += 1
            if b == "*":
                x += 1
            elif b[0].lower() == "a":
                a += 1
            elif b[0].lower() == "t":
                t += 1
            elif b[0].lower() == "c":
                c += 1
            elif b[0].lower() == "g":
                g += 1

            # It is important to count the insertions seperately
            if "+" in b:
                i += 1
        return coverage, a, t, c, g, x, i

    columns = ["pos", "coverage", "A", "T", "C", "G", "X", "I"]

    # 1 Is added to the position because our index starts at 1
    p_index = pd.DataFrame(
        (
            (p.pos + 1,) + parse_query_sequences(p.get_query_sequences(add_indels=True))
            for p in pileup
        ),
        columns=columns,
    )

    # Since the pileup does not return positions without any reads mapped, we have to
    # fill these with zeroes
    missing_positions = set(range(1, ref_length + 1)) - set(p_index.pos)
    missing_p_index = pd.DataFrame(
        ((i, 0, 0, 0, 0, 0, 0, 0) for i in missing_positions), columns=columns
    )
    p_index = p_index.append(missing_p_index).set_index("pos").sort_index()
    p_index.index.name = None

    return p_index
