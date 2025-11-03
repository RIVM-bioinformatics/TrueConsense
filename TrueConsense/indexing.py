from AminoExtract import SequenceReader, GFFDataFrame
import pandas as pd
import pysam


def Readbam(f):
    """Reads in a bam file and returns a pysam.AlignmentFile object

    Parameters
    ----------
    f
        the name of the bam file

    Returns
    -------
        A pysam.AlignmentFile object

    """
    return pysam.AlignmentFile(f, "rb")


def Gffindex(file: str) -> GFFDataFrame:
    """Reads in a GFF3 file and returns a pandas dataframe

    Parameters
    ----------
    file
        the path to the gff file

    Returns
    -------
        A dataframe

    """
    reader = SequenceReader(logger=None)
    return reader.read_gff(file)


def read_override_index(f):
    """Reads a csv file, and uses the first column as the index

    Parameters
    ----------
    f
        the file to read

    Returns
    -------
        A dataframe with the index column being the first column.

    """
    return pd.read_csv(f, sep=",", compression="gzip", index_col=0)


def Override_index_positions(index, override_data):
    """Takes a dataframe and a second dataframe with the same index and columns, and replaces the values
    in the first dataframe with the values in the second dataframe

    Parameters
    ----------
    index
        the index of the dataframe you want to override
    override_data
        a dataframe with the same columns as the index, but with the values you want to override

    Returns
    -------
        Dataframe with the overridden data.

    """
    index.loc[override_data.index, :] = override_data[:]
    return index


def BuildIndex(bamfile, ref):
    """Function takes a BAM file and a reference genome, and returns a dataframe pileup contents for each position in the bamfile.

    Parameters
    ----------
    bamfile
        The path to the bam file
    ref
        The reference genome

    Returns
    -------
        A dataframe with the following columns:
        pos: position in the reference genome
        coverage: number of reads covering the position
        A: number of reads with an A at the position
        T: number of reads with a T at the position
        C: number of reads with a C at the position
        G: number of reads with a G at the position

    """
    bamfile = pysam.AlignmentFile(bamfile, "rb")
    ref_fasta = pysam.FastaFile(ref)
    ref_length = ref_fasta.lengths[0]

    pileup = bamfile.pileup(stepper="nofilter", max_depth=10000000, min_base_quality=0)

    def parse_query_sequences(l):
        """Takes a list of strings, and returns a tuple of integers

        Parameters
        ----------
        l
            the list of bases in the query sequence

        Returns
        -------
            the coverage, a, t, c, g, x, and i values.

        """
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
    p_index = pd.concat([p_index, missing_p_index]).set_index("pos").sort_index()
    p_index.index.name = None

    return p_index
