from collections import Counter
import re


def ListInserts(iDict, mincov, bam):
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
    else:
        return True, positions


def ExtractInserts(bam, position):
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
        else:
            return None, None
    return None, None


def MinorityDel(index, p):
    cov = index[p].get("coverage")
    dels = index[p].get("X")
    perc = (dels / cov) * 100

    if perc >= 15:
        return True
    return False
