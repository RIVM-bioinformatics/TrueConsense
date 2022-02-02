import copy
import time
from itertools import chain

from .Ambig import IsAmbiguous
from .Coverage import GetCoverage
from .Events import ListInserts, MinorityDel
from .ORFs import CorrectGFF, SolveTripletLength, in_orf

import pandas as pd


def WalkForward(index, p, fixedpositions="expand"):

    if fixedpositions != "expand":
        lastposition = list(enumerate(index))[-1][1]
        p = p + 1
        targetposition = p + fixedpositions

        if targetposition >= lastposition:
            targetposition = lastposition

        nucleotrack = {}
        while p != targetposition:
            N, C = GetNucleotide(index, p, 1)
            nucleotrack[p] = N
            p += 1
        return nucleotrack
    else:
        track = []
        p = p + 1
        while True:
            N, C = GetNucleotide(index, p, 1)
            if N == "X":
                track.append(p)
                p += 1
            else:
                break
        return track


def complement_index(index, gffdict, skips):
    for i, p in enumerate(index):
        if p in skips:
            continue
        a = _orf_codonposition(gffdict, p)
        index[p]["ORF"] = a
    return index


def _orf_codonposition(gffdict, p):
    a = []
    for k in gffdict.keys():
        start = gffdict[k].get("start")
        end = gffdict[k].get("end") + 1

        if p in range(start, end):
            a.append(str(gffdict[k].get("attributes").split(";")[1].split("=")[-1]))
            a.append((p - start) % 3)
    if a:
        return tuple(a)
    return None, None


def _orf_overlapnumber(index, p):
    if _orf_hasoverlap(index, p) is True:
        return int(len(index[p].get("ORF")) / 2)
    return None


def _orf_hasoverlap(index, p):
    if len(index[p].get("ORF")) > 2:
        return True
    return False


def GetNucleotide(iDict, position, count):
    sorteddist = sorted(
        ((value, key) for key, value in GetDistribution(iDict, position).items())
    )
    return sorteddist[-count][1], sorteddist[-count][0]


def GetDistribution(iDict, position):
    dist = {}
    dist["A"] = iDict[position].get("A")
    dist["T"] = iDict[position].get("T")
    dist["C"] = iDict[position].get("C")
    dist["G"] = iDict[position].get("G")
    dist["X"] = iDict[position].get("X")
    return dist


def BuildConsensus(p_index, mincov=50, IncludeAmbig=False):

    print(f"Starting Consensus Building")
    start_time = time.time()

    p_index["calls"] = p_index[["query_sequences", "cov"]].apply(
        lambda l: call_counts(
            l.name,  # l.name is actually the query postition (the index of the df)
            l[0],
            l[1],
            mincov,
        ),
        axis=1,
    )
    print(f"Done applying call counts: {time.time() - start_time}")

    p_index.calls = p_index.calls.map(sort_highest_score)

    cons = p_index.calls.map(pick_first).map(lambda call: call["seq"] if call else "N")

    cons = "".join(cons)
    print(f"Done Building Consensus: {time.time() - start_time}")
    return cons


def call_counts(pos, query_sequences, cov, mincov):
    if len(query_sequences) < mincov:
        return []

    counts = pd.Series(query_sequences).str.upper().value_counts()
    scores = counts / cov

    return list(
        dict(pos=pos, seq=seq, n=n, score=score)
        for seq, n, score in zip(counts.index, counts, scores)
    )


def sort_highest_score(calls):
    calls.sort(key=lambda x: x["score"], reverse=True)
    return calls


def pick_first(calls):
    if not calls:
        return None
    return calls[0]
