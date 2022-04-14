def in_orf(loc, gffd):
    """If the location is in any of the ORFs, return True. Otherwise, return False

    Parameters
    ----------
    loc
        the current position
    gffd
        a dictionary of dictionaries, where the keys are the gene IDs, and the values are dictionaries containing the gene attributes

    Returns
    -------
        A list of True/False values.

    """
    exists = []
    for k in gffd.keys():
        start = gffd[k].get("start")
        stop = gffd[k].get("end")

        xx = loc in range(start, stop)
        exists.append(xx)

    if any(exists) == True:
        return True
    return False


def split_to_codons(seq):
    """It takes a string of DNA and returns a list of codons (nucleotide triplets)

    Parameters
    ----------
    seq
        the sequence to split

    Returns
    -------
        A list of codons

    """
    return [seq[start : start + 3] for start in range(0, len(seq), 3)]


def SolveTripletLength(uds, mds):
    """Check wether the combination of the uds (upcoming stretch of deletions) and the mds (group of minority deletions) is divisible by 3. If it is, return the length of the triplet. If it is not, return None.

    Parameters
    ----------
    uds
        upcoming stretch of dels
    mds
        the minority-del group

    Returns
    -------
        A boolean value.

    """
    mdslen = len(mds)
    udslen = len(uds)

    if udslen % 3 == 0:
        ## upcoming stretch is already divisible by 3

        if mdslen % 3 == 0:
            ## the minority-del group is also divisible by 3

            return True  # return that the minority-del is good to be appended
        ## the minority-del group is not divisible by 3
        return False  # return that the minority-del should be ignored
    else:
        ## the upcoming stretch is NOT divisible by 3 -> the mds should make it so
        if (mdslen + udslen) % 3 == 0:
            ## the combination of minority-del group and the upcoming stretch is divisible by 3.
            return True
        return False


def CorrectStartPositions(gffd, shifts, p):
    """Function takes a dictionary of gff data, a list of shifts, and a position. It then iterates
    through the dictionary and updates the start position of each entry if the start position is greater
    than the position

    Parameters
    ----------
    gffd
        a dictionary of dictionaries, where each key is a gene name, and each value is a dictionary of the
        gene's attributes
    shifts
        the number of bases to shift the start positions by
    p
        the position of the first base of the insertion

    Returns
    -------
        A dictionary with the updated start positions.

    """
    for k in gffd.keys():
        start = gffd[k].get("start")

        if start > p:
            nstart = int(start) + int(shifts)
            update = {"start": nstart}

            gffd[k].update(update)
    return gffd


def CorrectGFF(oldgffdict, newgffdict, cons, p, inserts, mincov, cov):
    """This function corrects the start and end positions of the genes in the new GFF file

    Parameters
    ----------
    oldgffdict
        a dictionary of the old gff file
    newgffdict
        the new gff dictionary
    cons
        the consensus sequence
    p
        the position of the current base
    inserts
        a dictionary of insertions, where the key is the position of the insertion and the value is the
    nucleotide inserted
    mincov
        minimum coverage to consider a position as a potential insertion
    cov
        coverage of the contig

    Returns
    -------
        A dictionary of the corrected gff file.

    """

    stopcodons = ["TAG", "TAA", "TGA"]
    # rvstopcodon = ["CAT"]

    if inserts is not None and p in inserts:
        if cov > mincov:
            newgffdict = CorrectStartPositions(
                newgffdict, list(inserts[p].keys())[0], p
            )

    for k in newgffdict.keys():
        start = newgffdict[k].get("start")
        end = newgffdict[k].get("end")
        orient = newgffdict[k].get("strand")

        shift = 0

        if p in range(start, end):
            if orient == "+":
                rseq = "".join(cons)[start - 1 :]
                shift = rseq.count("-")
                seq = rseq.replace("-", "")

                if cons[-1] == "-":
                    achieved = False
                    override_end = oldgffdict[k].get("end")
                    up = {"end": override_end}
                    newgffdict[k].update(up)

                else:

                    codons = split_to_codons(seq)

                    achieved = False

                    it = 0
                    for c in codons:
                        it += 1
                        if any(s in c for s in stopcodons) is True:
                            achieved = True
                            break
                    orfsize = it * 3

                    if shift % 3 == 0:
                        newend = start + orfsize + shift - 1
                    else:
                        newend = start + orfsize + shift - 1

                    if achieved is False:
                        newend = newend + 1

                    if p == newend:
                        up = {"end": newend}
                        newgffdict[k].update(up)

    return newgffdict
