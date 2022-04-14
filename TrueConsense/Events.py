import re
from collections import Counter


def ListInserts(iDict, mincov, bam):
    """This function takes a dictionary pileup contents, a minimum coverage value, and a bam file.
    returns a boolean and a dictionary of insert positions.

    Parameters
    ----------
    iDict
        This is the dictionary that contains pileup information including coverage and insertions for each position.
    mincov
        minimum coverage to take into account
    bam
        the bam file

    Returns
    -------
        A boolean value and a dictionary of positions with the insert size and nucleotide sequence

    """
    positions = {}

    for k in iDict.keys():
        cov = iDict[k].get("coverage")
        ins = iDict[k].get("I")

        if cov < mincov:
            continue
        if cov == 0 or ins == 0:
            continue
        if cov == 0 and ins == 0:
            continue
        perc = (ins / cov) * 100
        if perc > 55:
            InsNuc, insertsize = ExtractInserts(bam, k)
            if InsNuc is None or insertsize is None:
                continue
            positions[k] = {}
            positions[k][insertsize] = InsNuc
    if not positions:
        return False, None
    return True, positions


def ExtractInserts(bam, position):
    """It takes a bam file and a position and returns the most common insert sequence and the size of the
    insert

    Parameters
    ----------
    bam
        the bam file
    position
        the position in the reference genome to extract

    Returns
    -------
        the most common insert sequence and the size of the insert.

    """
    rname = bam.references[0]
    start = position - 1
    end = position
    for pileupcolumn in bam.pileup(rname, start, end, truncate=True):
        items = pileupcolumn.get_query_sequences(add_indels=True)
        found = []
        if bool(items) is False:
            return None, None
        for i in items:
            found.append(i.upper())
        sorteddist = dict(Counter(found).most_common())
        a = next(iter(sorteddist))
        match = re.search("(\d)([a-zA-Z]+)", a)

        if match:
            bases = match.group(2)
            insertsize = match.group(1)
            return bases, insertsize
        return None, None
    return None, None


def MinorityDel(index, p):
    """If the percentage of deletions on a position is greater than 15%, then return True

    Parameters
    ----------
    index
        The dictionary with pileup contents per position
    p
        position in the genome

    Returns
    -------
        True or False

    """
    cov = index[p].get("coverage")
    dels = index[p].get("X")
    perc = (dels / cov) * 100

    if perc >= 15:
        return True
    return False
