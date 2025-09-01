import copy

from .Ambig import IsAmbiguous
from .Coverage import GetCoverage
from .Events import ListInserts, MinorityDel
from .ORFs import CorrectGFF, SolveTripletLength, in_orf


def WalkForward(index, p, fixedpositions="expand"):
    """Function takes a dictionary of pileup data and current position, and returns a
    dictionary of future positions which contain deletions

    Parameters
    ----------
    index
        the index of the genome
    p
        the position of the nucleotide you want to start from
        fixedpositions, optional
        This is the number of nucleotides to walk forward. If you want to walk forward until you hit a
        nucleotide, set this to "expand".

    Returns
    -------
        A list of future positions which contain deletions.

    """

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
    """Takes a dictionary of pileup data, a dictionary of gff features, and a list of positions to skip, and returns a
    dictionary pileup data, with the addition of a new key, "ORF", which contains the which position in a codon the nucleotide-position has

    Parameters
    ----------
    index
        a dictionary of the form {protein_id: {'start': start, 'stop': stop, 'strand': strand, 'ORF': ORF}}
    gffdict
        a dictionary of the gff file, with the keys being the gene names and the values being the start and
        stop positions of the gene.
    skips
        a list of the names of the genes that you want to skip.

    Returns
    -------
        The index is being returned.

    """
    for i, p in enumerate(index):
        if p in skips:
            continue
        a = _orf_codonposition(gffdict, p)
        index[p]["ORF"] = a
    return index


def _orf_codonposition(gffdict, p):
    """Takes a dictionary of ORFs and a position, and returns the ORF and the codon position of that
    position

    Parameters
    ----------
    gffdict
        a dictionary of gff features
    p
        position in the genome

    Returns
    -------
        a tuple of the ORF name and the codon position of the position p.

    """
    a = []
    for k in gffdict.keys():
        start = gffdict[k].get("start")
        end = gffdict[k].get("end") + 1

        if p in range(start, end):
            # attributes itself can also be "", so both no attributes and empty attributes are skipped
            attr = gffdict[k].get("attributes", "")
            if not attr:
                continue
            a.append(str(attr.split(";")[1].split("=")[-1]))
            a.append((p - start) % 3)
    if a:
        return tuple(a)
    return None, None


def GetNucleotide(iDict, position, count):
    """Takes a dictionary of sequences, a position, and a count, and returns the nucleotide at that
    position that occurs the most, and the number of times it occurs

    Parameters
    ----------
    iDict
        The dictionary of of pileup data
    position
        the position in the sequence you want to get the nucleotide for
    count
        the number of nucleotides you want to return

    Returns
    -------
        The nucleotide and the frequency of that nucleotide at a given position.

    """
    sorteddist = sorted(
        ((value, key) for key, value in GetDistribution(iDict, position).items())
    )
    return sorteddist[-count][1], sorteddist[-count][0]


def GetDistribution(iDict, position):
    """Takes a dictionary of dictionaries and a position and returns a dictionary of the distribution of
    nucleotides at that position

    Parameters
    ----------
    iDict
        the dictionary pileup data
    position
        the position in the sequence that you want to get the distribution for

    Returns
    -------
        A dictionary of the distribution of nucleotides at a given position.

    """
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
