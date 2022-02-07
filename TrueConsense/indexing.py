import gffpandas.gffpandas as gffpd
import pandas as pd
import pysam
import time


def ReadBam(f):
    return pysam.AlignmentFile(f, "rb")


def ReadFasta(f):
    return pysam.FastaFile(f)


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

    print(f"Starting building index")
    start_time = time.time()

    pileup = bamfile.pileup(stepper="nofilter", max_depth=10000000, min_base_quality=0)

    print(f"Done reading bamfile: {time.time() - start_time}")
    start_time = time.time()

    columns = ["pos", "query_sequences"]
    # 1 Is added to the position because our index starts at 1
    p_index = pd.DataFrame(
        (
            (
                p.pos + 1,
                [
                    call.split("-")[0].replace("*", "-").replace("+", "")
                    for call in p.get_query_sequences(add_indels=True)
                    if call
                ],
            )
            for p in pileup
        ),
        columns=columns,
    )
    print(f"Done parsing query sequences: {time.time() - start_time}")
    start_time = time.time()

    # Since the pileup does not return positions without any reads mapped, we have to
    # fill these with zeroes
    missing_positions = set(range(1, ref_length + 1)) - set(p_index.pos)
    missing_p_index = pd.DataFrame(
        ((i, []) for i in missing_positions), columns=columns
    )
    p_index = p_index.append(missing_p_index).set_index("pos").sort_index()
    p_index.index.name = None
    print(f"Done filling missing positions in pileup: {time.time() - start_time}")
    start_time = time.time()

    p_index = SupplyIndexWithCounts(p_index)
    print(f"Done supplying coverage levels: {time.time() - start_time}")
    start_time = time.time()

    return p_index


def SupplyIndexWithCounts(p_index):
    p_index["cov"] = p_index.query_sequences.apply(len)
    return p_index
