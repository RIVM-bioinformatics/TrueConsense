from .Inserts import ListInserts, ExtractInserts
from .Coverage import GetCoverage
from .Ambig import IsAmbiguous


def Inside_ORF(loc, gff):
    exists = []
    for k in gff.keys():
        start = gff[k].get("start")
        stop = gff[k].get("end")

        in_orf = loc in range(start, stop)
        exists.append(in_orf)

    if any(exists) == True:
        return True
    else:
        return False


def BuildConsensus(mincov, iDict, UpdatedGFF, bam, withambig, includeIns):

    consensus = []

    hasinserts, insertpositions = ListInserts(iDict, mincov)

    def previousposition(loc, step):
        if loc == 1:
            return 1
        if loc == 2:
            return 1
        if loc > 2:
            return loc - step

    def nextposition(loc, step, last):
        if loc == last:
            return loc
        elif loc == (last - 1):
            return loc + 1
        elif loc == (last - 2):
            return loc + step
        else:
            return loc + step

    lastposition = list(range(len(iDict)))[-1] + 1
    for i in range(len(iDict)):
        currentposition = i + 1

        prev1pos = previousposition(currentposition, 1)
        prev2pos = previousposition(currentposition, 2)
        next1pos = nextposition(currentposition, 1, lastposition)
        next2pos = nextposition(currentposition, 2, lastposition)

        within_orf = Inside_ORF(currentposition, UpdatedGFF)

        cov = GetCoverage(iDict, currentposition)

        (
            cur_FirstNuc,
            cur_SecondNuc,
            cur_ThirdNuc,
            cur_FourthNuc,
            cur_FifthNuc,
            cur_FirstCount,
            cur_SecondCount,
            cur_ThirdCount,
            cur_FourthCount,
            cur_FifthCount,
        ) = GetProminentNucleotides(iDict, currentposition)

        (
            prv1_FirstNuc,
            prv1_SecondNuc,
            prv1_ThirdNuc,
            prv1_FourthNuc,
            prv1_FifthNuc,
            prv1_FirstCount,
            prv1_SecondCount,
            prv1_ThirdCount,
            prv1_FourthCount,
            prv1_FifthCount,
        ) = GetProminentNucleotides(iDict, prev1pos)

        (
            prv2_FirstNuc,
            prv2_SecondNuc,
            prv2_ThirdNuc,
            prv2_FourthNuc,
            prv2_FifthNuc,
            prv2_FirstCount,
            prv2_SecondCount,
            prv2_ThirdCount,
            prv2_FourthCount,
            prv2_FifthCount,
        ) = GetProminentNucleotides(iDict, prev2pos)

        (
            nxt1_FirstNuc,
            nxt1_SecondNuc,
            nxt1_ThirdNuc,
            nxt1_FourthNuc,
            nxt1_FifthNuc,
            nxt1_FirstCount,
            nxt1_SecondCount,
            nxt1_ThirdCount,
            nxt1_FourthCount,
            nxt1_FifthCount,
        ) = GetProminentNucleotides(iDict, next1pos)

        (
            nxt2_FirstNuc,
            nxt2_SecondNuc,
            nxt2_ThirdNuc,
            nxt2_FourthNuc,
            nxt2_FifthNuc,
            nxt2_FirstCount,
            nxt2_SecondCount,
            nxt2_ThirdCount,
            nxt2_FourthCount,
            nxt2_FifthCount,
        ) = GetProminentNucleotides(iDict, next2pos)

        HasAmbiguity, AmbigCharacter = IsAmbiguous(
            cur_FirstNuc,
            cur_SecondNuc,
            cur_ThirdNuc,
            cur_FourthNuc,
            cur_FirstCount,
            cur_SecondCount,
            cur_ThirdCount,
            cur_FourthCount,
            cov,
        )

        if cov < mincov:
            consensus.append("N")
        else:
            if cur_FirstNuc.upper() != "X":
                if withambig is True and HasAmbiguity is True:
                    consensus.append(AmbigCharacter)
                else:
                    if cur_FirstCount < mincov:
                        consensus.append(cur_FirstNuc.lower())
                    else:
                        consensus.append(cur_FirstNuc.upper())
            elif cur_FirstNuc.upper() == "X":
                if within_orf is False:
                    consensus.append("-")
                elif within_orf is True:

                    if (
                        IsRealDel(
                            cur_FirstNuc,
                            prv1_FirstNuc,
                            prv2_FirstNuc,
                            nxt1_FirstNuc,
                            nxt2_FirstNuc,
                        )
                        is True
                    ):
                        consensus.append("-")
                    else:
                        secondary_hasambig, sec_ambigchar = IsAmbiguous(
                            cur_SecondNuc,
                            cur_ThirdNuc,
                            cur_FourthNuc,
                            cur_FifthNuc,
                            cur_SecondCount,
                            cur_ThirdCount,
                            cur_FourthCount,
                            cur_FifthCount,
                            cov - cur_FirstCount,
                        )
                        if withambig is True and secondary_hasambig is True:
                            consensus.append(sec_ambigchar)
                        else:
                            if cur_SecondCount != 0:
                                if cur_SecondCount >= mincov:
                                    consensus.append(cur_SecondNuc)
                                else:
                                    consensus.append(cur_SecondNuc.lower())
                            elif cur_SecondNuc == 0:
                                consensus.append("N")
        if includeIns is True:
            if cov > mincov:
                if hasinserts is True:
                    for x in insertpositions:
                        if currentposition == x:
                            try:
                                InsertNuc, InsertSize = ExtractInserts(
                                    bam, currentposition
                                )
                                if InsertNuc is not None:
                                    consensus.append(InsertNuc)
                                else:
                                    continue
                            except:
                                continue

    return "".join(consensus)


def IsRealDel(c, p1, p2, n1, n2):
    if c != "X":
        return False
    if n1.upper() == "X" and p1.upper() == "X":
        return True
    elif n1.upper() == "X" and n2.upper() == "X":
        return True
    elif p1.upper() == "X" and p2.upper() == "X":
        return True
    else:
        return False


def split_to_codons(seq, num):
    return [seq[start : start + num] for start in range(0, len(seq), 3)]


def UpdateGFF(mincov, iDict, bam, GffDict):
    draftseq, insertions = DraftConsensus(mincov, iDict, bam)
    stopcodons = ["TAG", "TAA", "TGA"]

    shift = 0
    for k in GffDict.keys():
        if (GffDict[k].get("strand")) is "+":  ## 'forward' orientation of ORF
            start = GffDict[k].get("start")
            start = (
                start - 1 + shift
            )  ## adjust for 0-based + adjust for nucleotide shifts due to insertions
            try:
                nextstart = GffDict[k + 1].get("start")
                nextstart = nextstart - 1
            except:
                nextstart = None

            codons = split_to_codons(draftseq[start:], 3)

            i = 0
            for c in codons:
                i += 1
                if any(stops in c for stops in stopcodons) is True:
                    break

            orfsize = i * 3
            end = start + orfsize

            for p in insertions.keys():
                if nextstart is not None:
                    if p in range(start, nextstart):
                        localshift = int(insertions[p])
                        # GffDict[k+1]["start"] = int(nextstart) + int(localshift) + 1
                        shift += int(localshift)
                        # print(GffDict[k+1]["start"])

            GffDict[k]["start"] = start + 1
            GffDict[k]["end"] = end

        if (GffDict[k].get("strand")) is "-":  ## 'reverse' orientation of ORF
            continue

    return GffDict


def DraftConsensus(mincov, iDict, bam):

    draft = []
    insertions = {}

    hasinserts, insertpositions = ListInserts(iDict, mincov)
    for i in range(len(iDict)):
        realpos = i + 1

        cov = GetCoverage(iDict, realpos)
        (
            FirstNuc,
            SecondNuc,
            ThirdNuc,
            FourthNuc,
            FifthNuc,
            FirstCount,
            SecondCount,
            ThirdCount,
            FourthCount,
            FifthCount,
        ) = GetProminentNucleotides(iDict, realpos)

        if cov < mincov:
            draft.append("N")
        else:
            if FirstNuc == "X":
                draft.append("-")
            else:
                draft.append(FirstNuc)

        if cov > mincov:
            if hasinserts is True:
                for lposition in insertpositions:
                    if realpos == lposition:
                        try:
                            insertNucs, insertsize = ExtractInserts(bam, realpos)
                            if insertNucs is not None:
                                draft.append(insertNucs)
                                insertions[realpos] = insertsize
                            else:
                                continue
                        except:
                            continue

    return "".join(draft), insertions


def GetProminentNucleotides(iDict, position):
    sorteddist = sorted(
        ((value, key) for key, value in GetDistribution(iDict, position).items())
    )
    # sorteddist[-1][1] == most prominent nucleotide
    # sorteddist[-2][1] == second most prominent nucleotide
    # sorteddist[-3][1] == third most prominent nucleotide
    ## etc etc
    # sorteddist[-1][0] == counter for most prominent nucleotide
    # sorteddist[-2][0] == counter for second most prominent nucleotide
    ## etc etc
    return (
        sorteddist[-1][1],
        sorteddist[-2][1],
        sorteddist[-3][1],
        sorteddist[-4][1],
        sorteddist[-5][1],
        sorteddist[-1][0],
        sorteddist[-2][0],
        sorteddist[-3][0],
        sorteddist[-4][0],
        sorteddist[-5][0],
    )


def GetDistribution(iDict, position):
    dist = {}
    dist["A"] = iDict[position].get("A")
    dist["T"] = iDict[position].get("T")
    dist["C"] = iDict[position].get("C")
    dist["G"] = iDict[position].get("G")
    dist["X"] = iDict[position].get("X")
    return dist
