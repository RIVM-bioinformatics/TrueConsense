import copy
import time
from itertools import *
from collections import Counter
from turtle import position
import numpy as np

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


def BuildConsensus(p_index, mincov=50, IncludeAmbig=False, gff_df=None):
    print(f"Starting Consensus Building")
    start_time = time.time()

    p_index["calls"] = p_index[["query_sequences"]].apply(
        lambda l: call_counts(
            l.name,  # l.name is actually the query postition (the index of the df)
            l["query_sequences"],
            mincov,
        ),
        axis=1,
    )
    print(f"Done determining call counts: {time.time() - start_time}")
    start_time = time.time()

    p_index["calls"] = p_index.calls.map(sort_highest_score)
    p_index["picked_call"] = p_index["calls"].map(lambda c: c[0] if c else None)

    # The significance is the minimal relative score that (combinations of) mutations need to have before being considered.
    significance = 0.5
    alt_calls = [
        alt_call
        for alt_calls in p_index.calls
        for alt_call in alt_calls[1:]
        if alt_call["rel_score"] > significance
    ]
    alt_calls.sort(key=lambda call: call["rel_score"], reverse=True)

    for _, feature in gff_df.iterrows():
        feature_score = score_feature(p_index, feature)
        print(f"Fixing {feature.Name} with score {feature_score}")

        best_calls = []
        for new_calls in significant_combinations_of_mutations(
            [c for c in alt_calls if feature.start <= c["pos"] <= feature.end],
            significance=significance,
        ):
            if feature_score > 0.99:
                break

            previous_calls = insert_calls(p_index, new_calls)
            new_feature_score = score_feature(p_index, feature)
            if (
                new_feature_score > feature_score
            ):  # TODO: Add weighing of the alternative calls to the importance of the feature
                print(
                    f"Fixed to score of {new_feature_score} with {new_calls} in place of {previous_calls}"
                )
                feature_score = new_feature_score
                best_calls = new_calls
            insert_calls(p_index, previous_calls)  # restore to default
        insert_calls(p_index, best_calls)
        print(f"Final score of {feature.Name} is {feature_score}")
    print(f"Done fixing features: {time.time() - start_time}")

    start_time = time.time()

    cons = p_index.picked_call.map(lambda call: call["seq"] if call else "N")

    cons = "".join(cons)
    print(f"Done Building Consensus: {time.time() - start_time}")
    return cons


def call_counts(pos, query_sequences, mincov):
    cov = len(query_sequences)
    if cov < mincov:
        return []

    counts = pd.Series(Counter(map(str.upper, query_sequences)))

    scores = counts / cov

    # This is (kind of) the probability that you would choose the alternative call(s) over the most occuring call(s).
    # e.g. there are 5 calls for the nucleotide C, and 3 for T. The relative score for C is 1.0 and for T 3/5=0.6.
    # These relative scores multiply in case of combinations of mutations.
    relative_scores = counts / max(counts)

    return list(
        dict(
            pos=pos,
            seq=seq,
            n=n,
            score=score,
            rel_score=rel_score,
            index=index,
        )
        for seq, n, score, rel_score, index in zip(
            counts.index,
            counts,
            scores,
            relative_scores,
            range(len(counts)),
        )
    )


def sort_highest_score(calls):
    calls.sort(key=lambda x: x["score"], reverse=True)
    return calls


def score_feature(p_index, feature, ends_with_stop_codon=True, starts_with_atg=True):
    out_of_frame_counter = 0
    in_frame_counter = 0

    encountered_stop = False
    frame_offset = 0
    position_in_codon = 0
    last_two_nucs = ""
    for i, call in enumerate(p_index[feature.start - 1 : feature.end]["picked_call"]):
        if not call:  # If there is an N
            # TODO: fix if there is a framshift (insertion or deletion) after a stretch of N's. An ins or del could shift the frame in step, in stead of out of step.
            frame_offset = 0  # Assume you are back in frame if insuffiecient info
            last_two_nucs = ""
            position_in_codon = (
                i + 1
            ) % 3  # Reset the position in codon to be in phase
            continue

        seq = call["seq"]
        if seq == "-":
            frame_offset = frame_offset - 1
            continue

        seq_with_prev = last_two_nucs + seq
        last_two_nucs = seq_with_prev[-2:]
        codons = chunks(seq_with_prev[2 - position_in_codon :], 3)
        for codon in codons:
            if len(codon) == 3:
                if (
                    encountered_stop or frame_offset % 3
                ):  # If the offset is not a multiple of 3
                    out_of_frame_counter += 1
                else:
                    in_frame_counter += 1
                    if codon in ["TAA", "TGA", "TAG"]:
                        encountered_stop = True

        position_in_codon = (position_in_codon + len(seq)) % 3

        if len(seq) > 1:
            frame_offset = frame_offset + len(seq) - 1

    score = in_frame_counter / (in_frame_counter + out_of_frame_counter)
    return score


def significant_combinations_of_mutations(ms, significance=0.5):
    key = lambda x: x["rel_score"]
    ms = sorted(ms, reverse=True, key=key)

    def score(l):
        return np.prod(list(map(key, l)))

    def pred(l):
        return score(l) > significance

    its = []
    for combination_of_length in map(
        lambda n: combinations(ms, n), range(1, len(ms) + 1)
    ):
        its.append(takewhile(pred, combination_of_length))
    return sorted(chain(*its), reverse=True, key=score)


def insert_calls(p_index, calls):
    list_of_originals = []
    for call in calls:
        list_of_originals.append(p_index.at[call["pos"], "picked_call"])
        p_index.at[call["pos"], "picked_call"] = call
    return list_of_originals


def fix_feature(p_index, feature, alternative_calls):
    pass


def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i : i + n]
