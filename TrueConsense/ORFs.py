import time
from itertools import *
import numpy as np


def RestoreORFS(call_obj, gff_df):
    start_time = time.time()

    # The significance is the minimal relative score that (combinations of) mutations need to have before being considered.
    significance = call_obj.significance

    alt_calls = [
        alt_call
        for alt_calls in call_obj.p_index.calls
        if alt_calls is not np.nan
        for alt_call in alt_calls[1:]
        # TODO: Make this work for non-sorted picked calls
        if alt_call["rel_score"] > significance
        and alt_call["n"] > 1  # TODO: This is for performance
    ]
    alt_calls.sort(key=lambda call: call["rel_score"], reverse=True)

    for _, feature in gff_df.iterrows():
        if feature["type"].upper() not in [
            "GENE",
            "CDS",
        ]:  # TODO: Switch this to only CDS in the future as gff from RefSeq is used
            continue
        # TODO: fix for overlapping features (they should be scored together, as calls in one feature can affect another)
        feature_score = call_obj.score_coding_sequence(feature)

        best_calls = []
        for new_calls in significant_combinations_of_calls(
            [c for c in alt_calls if feature.start <= c["pos"] <= feature.end],
            significance=significance,
        ):
            if feature_score > 0.99:
                break

            previous_calls = call_obj.insert_calls(new_calls)
            new_feature_score = call_obj.score_coding_sequence(feature)
            if new_feature_score > feature_score:
                # TODO: Add weighing of the significance of alternative calls vs. the importance of the feature
                print(
                    f"Fixed to score of {new_feature_score} with {new_calls} in place of {previous_calls}"
                )
                feature_score = new_feature_score
                best_calls = new_calls
            call_obj.insert_calls(previous_calls)  # restore to default
        call_obj.insert_calls(best_calls)
    print(f"Done fixing features: {time.time() - start_time}")


def significant_combinations_of_calls(calls, significance=0.5):
    """Scores and sorts alternative calls

    In the list of input calls, find combinations of calls that go over a specified significance level.

    Args:
        calls (iterable of calls): these calls should have a property rel_score.
        significance (float, optional): the product of relative scores that a combination of calls should not exceed. Defaults to 0.5.

    Returns:
        list(list(call)): A list of combinations of calls that have a minimal significance (sorted by significance)
    """
    key = lambda x: x["rel_score"]
    calls = sorted(calls, reverse=True, key=key)

    def score(l):
        """The score of a list of calls is the product of their relative score."""
        return np.prod(list(map(key, l)))

    def pred(l):
        return score(l) > significance

    iterators = []
    # Loop over all the combinations of calls of all lengths
    # TODO: Alternatives with the same position should not be combined.
    for combination_of_length in map(
        lambda n: combinations(calls, n), range(1, len(calls) + 1)
    ):
        # Stop taking combinations if pred is no longer. This only works if calls is sorted on score.
        iterators.append(takewhile(pred, combination_of_length))

    return sorted(chain(*iterators), reverse=True, key=score)


def chunks(l, n):
    """Generator for splitting a list into chunks.

    Args:
        l (list(any)): The input list
        n (int): The size of the chunks

    Yields:
        list(any): sublists of size n
    """
    for i in range(0, len(l), n):
        yield l[i : i + n]


# def in_orf(loc, gffd):
#     exists = []
#     for k in gffd.keys():
#         start = gffd[k].get("start")
#         stop = gffd[k].get("end")

#         xx = loc in range(start, stop)
#         exists.append(xx)

#     if any(exists) == True:
#         return True
#     return False


# def split_to_codons(seq):
#     return [seq[start : start + 3] for start in range(0, len(seq), 3)]


# def SolveTripletLength(uds, mds):
#     mdslen = len(mds)
#     udslen = len(uds)

#     if udslen % 3 == 0:
#         ## upcoming stretch is already divisible by 3

#         if mdslen % 3 == 0:
#             ## the minority-del group is also divisible by 3

#             return True  # return that the minority-del is good to be appended
#         ## the minority-del group is not divisible by 3
#         return False  # return that the minority-del should be ignored
#     else:
#         ## the upcoming stretch is NOT divisible by 3 -> the mds should make it so
#         if (mdslen + udslen) % 3 == 0:
#             ## the combination of minority-del group and the upcoming stretch is divisible by 3.
#             return True
#         return False


# def CorrectStartPositions(gffd, shifts, p):
#     for k in gffd.keys():
#         start = gffd[k].get("start")

#         if start > p:
#             nstart = int(start) + int(shifts)
#             update = {"start": nstart}

#             gffd[k].update(update)
#     return gffd


# def CorrectGFF(oldgffdict, newgffdict, cons, p, inserts, mincov, cov):

#     stopcodons = ["TAG", "TAA", "TGA"]
#     # rvstopcodon = ["CAT"]

#     if inserts is not None and p in inserts:
#         if cov > mincov:
#             newgffdict = CorrectStartPositions(
#                 newgffdict, list(inserts[p].keys())[0], p
#             )

#     for k in newgffdict.keys():
#         start = newgffdict[k].get("start")
#         end = newgffdict[k].get("end")
#         orient = newgffdict[k].get("strand")

#         shift = 0

#         if p in range(start, end):
#             if orient == "+":
#                 rseq = "".join(cons)[start - 1 :]
#                 shift = rseq.count("-")
#                 seq = rseq.replace("-", "")

#                 if cons[-1] == "-":
#                     achieved = False
#                     override_end = oldgffdict[k].get("end")
#                     up = {"end": override_end}
#                     newgffdict[k].update(up)

#                 else:

#                     codons = split_to_codons(seq)

#                     achieved = False

#                     it = 0
#                     for c in codons:
#                         it += 1
#                         if any(s in c for s in stopcodons) is True:
#                             achieved = True
#                             break
#                     orfsize = it * 3

#                     if shift % 3 == 0:
#                         newend = start + orfsize + shift - 1
#                     else:
#                         newend = start + orfsize + shift - 1

#                     if achieved is False:
#                         newend = newend + 1

#                     if p == newend:
#                         up = {"end": newend}
#                         newgffdict[k].update(up)

#     return newgffdict
