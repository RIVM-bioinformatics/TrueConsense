import copy

from .Ambig import IsAmbiguous
from .Coverage import GetCoverage
from .Events import ListInserts, MinorityDel
from .ORFs import CorrectGFF, SolveTripletLength, in_orf


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


def BuildConsensus(mincov, iDict, GFFdict, IncludeAmbig, bam, includeINS):
    cons = []

    p_index = complement_index(iDict, GFFdict, [])

    hasinserts, insertpositions = ListInserts(p_index, mincov, bam)

    dskips = []

    newGffdict = copy.deepcopy(GFFdict)

    for a, b in enumerate(p_index):

        cov = GetCoverage(p_index, b)
        within_orf = in_orf(b, newGffdict)

        if b in dskips:
            cons.append("-")
            newGffdict = CorrectGFF(
                GFFdict, newGffdict, cons, b, insertpositions, mincov, cov
            )
            continue

        if cov < mincov:
            cons.append("N")
            # Simply add a 'N' to the consensus at this position if the coverage is below the threshold
            newGffdict = CorrectGFF(
                GFFdict, newGffdict, cons, b, insertpositions, mincov, cov
            )
            continue
        else:
            PrimaryN, PrimaryC = GetNucleotide(p_index, b, 1)
            # get the primary nucleotide at this position

            HasAmbiguity, AmbigChar = IsAmbiguous(
                GetNucleotide(p_index, b, 1),
                GetNucleotide(p_index, b, 2),
                GetNucleotide(p_index, b, 3),
                GetNucleotide(p_index, b, 4),
                cov,
            )

            if PrimaryN != "X":

                if MinorityDel(p_index, b) is True:

                    if not WalkForward(p_index, b):

                        if MinorityDel(p_index, b + 1) is True:

                            if WalkForward(p_index, b + 1):

                                uds = WalkForward(p_index, b + 1)
                                mds = [b, b + 1]
                                if SolveTripletLength(uds, mds) is True:
                                    cons.append("-")
                                    for x in mds:
                                        dskips.append(x)
                                    for x in uds:
                                        dskips.append(x)
                                else:
                                    if IncludeAmbig is True and HasAmbiguity is True:
                                        cons.append(AmbigChar)
                                    else:
                                        if PrimaryC < mincov:
                                            cons.append(PrimaryN.lower())
                                        else:
                                            cons.append(PrimaryN.upper())
                            else:
                                if IncludeAmbig is True and HasAmbiguity is True:
                                    cons.append(AmbigChar)
                                else:
                                    if PrimaryC < mincov:
                                        cons.append(PrimaryN.lower())
                                    else:
                                        cons.append(PrimaryN.upper())
                        else:
                            if IncludeAmbig is True and HasAmbiguity is True:
                                cons.append(AmbigChar)
                            else:
                                if PrimaryC < mincov:
                                    cons.append(PrimaryN.lower())
                                else:
                                    cons.append(PrimaryN.upper())
                    else:
                        uds = WalkForward(p_index, b)
                        mds = [b]
                        if SolveTripletLength(uds, mds) is True:
                            cons.append("-")
                            dskips.append(b)
                            for x in uds:
                                dskips.append(x)
                        else:
                            if IncludeAmbig is True and HasAmbiguity is True:
                                cons.append(AmbigChar)
                            else:
                                if PrimaryC < mincov:
                                    cons.append(PrimaryN.lower())
                                else:
                                    cons.append(PrimaryN.upper())
                else:
                    if IncludeAmbig is True and HasAmbiguity is True:
                        cons.append(AmbigChar)
                    else:
                        if PrimaryC < mincov:
                            cons.append(PrimaryN.lower())
                        else:
                            cons.append(PrimaryN.upper())

            if PrimaryN == "X":
                if within_orf is True:
                    if b in dskips:  # this might be redundant, check later
                        cons.append("-")
                    else:
                        if not WalkForward(p_index, b):
                            SecondaryN, SecondaryC = GetNucleotide(p_index, b, 2)
                            if IncludeAmbig is True and HasAmbiguity is True:
                                cons.append(AmbigChar)
                            else:
                                if SecondaryC < mincov:
                                    cons.append(SecondaryN.lower())
                                else:
                                    cons.append(SecondaryN.upper())
                        else:
                            wfds = WalkForward(p_index, b)
                            if len(wfds) >= 2:
                                cons.append("-")
                                dskips.append(b)
                                for x in wfds:
                                    dskips.append(x)
                            else:
                                SecondaryN, SecondaryC = GetNucleotide(p_index, b, 2)
                                if IncludeAmbig is True and HasAmbiguity is True:
                                    cons.append(AmbigChar)
                                else:
                                    if SecondaryC < mincov:
                                        cons.append(SecondaryN.lower())
                                    else:
                                        cons.append(SecondaryN.upper())
                else:
                    cons.append("-")

        if includeINS is True:
            if cov > mincov:
                if hasinserts is True:
                    for x in insertpositions:
                        if b == x:
                            for i in insertpositions.get(x):
                                cons.append(str(insertpositions.get(x).get(i)))

        newGffdict = CorrectGFF(
            GFFdict, newGffdict, cons, b, insertpositions, mincov, cov
        )

    return "".join(cons), newGffdict
