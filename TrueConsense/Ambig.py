"""
double ambiguities:
    "M": ["A", "C"],
    "R": ["A", "G"],
    "W": ["A", "T"],
    "S": ["C", "G"],
    "Y": ["C", "T"],
    "K": ["G", "T"],

triplet ambiguities:
    "V": ["A", "C", "G"],
    "H": ["A", "C", "T"],
    "D": ["A", "G", "T"],
    "B": ["C", "G", "T"],
"""


def DoubleAmbigs(n1, n2):
    """Takes primary and secondary nucleotides and returns the corresponding ambiguity nucleotide code.

    Parameters
    ----------
    n1
        the primary nucleotide
    n2
        the secondary nucleotide

    Returns
    -------
        the nucleotide ambiguity code.

    """
    nlist = [n1, n2]

    if (
        any(x == "A" for x in nlist) is True
    ):  ## possibilities are M/R/W --> A + ? (see table above)
        if any(x == "C" for x in nlist) is True:  ## only option is M --> A + C
            char = "M"
        elif any(x == "G" for x in nlist) is True:  ## only option is R --> A + G
            char = "R"
        elif any(x == "T" for x in nlist) is True:  ## only option is W --> A + T
            char = "W"
    elif (
        any(x == "C" for x in nlist) is True
    ):  ## possibilities are S/Y --> C + ? (see table above)
        if any(x == "G" for x in nlist) is True:  ## only option is S --> C + G
            char = "S"
        elif any(x == "T" for x in nlist) is True:  ## only option is Y --> C + T
            char = "Y"
    elif (
        any(x == "G" for x in nlist) is True
    ):  ## only option is K due to exclusion (see table above)
        char = "K"
    return char


def TripletAmbigs(n1, n2, n3):
    """Takes primary, secondary and tertiary nucleotides and returns the corresponding ambiguity nucleotide code.

    Parameters
    ----------
    n1
        the primary nucleotide
    n2
        the secondary nucleotide
    n3
        the tertiary nucleotide

    Returns
    -------
        the nucleotide ambiguity code.

    """
    nlist = [n1, n2, n3]
    if (
        any(x == "A" for x in nlist) is True
    ):  ## possibilities are V/H/D --> A + ? + ? (see table above)
        if (
            any(x == "C" for x in nlist) is True
        ):  ## possibilities are V/H --> A + C + ? (see table above)
            if (
                any(x == "G" for x in nlist) is True
            ):  ## only possibility left is V --> A + C + G (see table above)
                char = "V"
            elif (
                any(x == "T" for x in nlist) is True
            ):  ## only possibility left is H --> A + C + T (see table above)
                char = "H"
        elif (
            any(x == "T" for x in nlist) is True
        ):  ## only possibility left is D due to exclusion --> A + T + ? (see table above)
            char = "D"
    elif (
        any(x == "T" for x in nlist) is True
    ):  ## only possibility left is B due to exclusion --> T + ? + ? (see table above)
        char = "B"

    return char


def GetPercentages(c1, c2, c3, c4, cov):
    """Takes the number of each nucleotide and transforms it to a percentage in relation to the coverage

    Parameters
    ----------
    c1
        Integer representing the count of the first nucleotide on a certain position
    c2
        Integer representing the count of the secondary nucleotide on a certain position
    c3
        Integer representing the count of the third nucleotide on a certain position
    c4
        Integer representing the count of the fourth nucleotide on a certain position
    cov
        total coverage on a certain position

    Returns
    -------
        The percentages of each nucleotide in the sequence.

    """
    p1 = (c1 / cov) * 100
    p2 = (c2 / cov) * 100
    p3 = (c3 / cov) * 100
    p4 = (c4 / cov) * 100
    return p1, p2, p3, p4


def AmbiguityType(p1, p2, p3, p4):
    """If the distance between the first two percentages is less than 10, then check if the distance between the
    first and third percentage is less than 10, and if the distance between the second and third percentage is
    less than 10. If all of these are true, then check if the distance between the first and fourth
    percentage is less than 10, and if the distance between the second and fourth percentage is less than 10,
    and if the distance between the third and fourth percentage is less than 10. If all of these are true,
    then the ambiguity type is 4. If the first three are true, but the last three are not, then the
    ambiguity type is 3. If the first two are true, but the last three are not, then the ambiguity type
    is 2. If the first two are not true, then the ambiguity type is None

    Parameters
    ----------
    p1
        the percentage of the first nucleotide
    p2
        the percentage of the secondary nucleotide
    p3
        the percentage of the third nucleotide
    p4
        the percentage of the fourth nucleotide

    Returns
    -------
        the Ambiguity Type.

    """
    maxdistance = 10
    if (abs(p1 - p2)) <= maxdistance:
        if (abs(p1 - p3)) <= maxdistance and (abs(p2 - p3)) <= maxdistance:
            if (
                (abs(p1 - p4)) <= maxdistance
                and (abs(p2 - p4)) <= maxdistance
                and (abs(p3 - p4)) <= maxdistance
            ):
                AType = 4
            else:
                AType = 3
        else:
            AType = 2
    else:
        AType = None
    return AType


def unpack(q):
    # todo: useless function?
    return q[0], q[1]


def IsAmbiguous(one, two, three, four, cov):
    """If the coverage is 0, return False, None. If the coverage is not 0, unpack the counts and get the
    percentages. If the percentages are not ambiguous, return False, None. If the percentages are
    ambiguous, return True, the ambiguous character.

    Parameters
    ----------
    one
        the first nucleotide
    two
        the first nucleotide
    three
        the nucleotide counts for the third position in the codon
    four
        a tuple of the form (nuc, count)
    cov
        coverage

    Returns
    -------
        A tuple of two values. The first value is a boolean, which is True if the site is ambiguous, and
    False if it is not. The second value is a string, which is the ambiguous character if the site is
    ambiguous, and None if it is not.

    """
    if cov == 0:
        return False, None

    nuc1, count1 = unpack(one)
    nuc2, count2 = unpack(two)
    nuc3, count3 = unpack(three)
    nuc4, count4 = unpack(four)

    if nuc1 == "X" or nuc2 == "X":
        return False, None
    p1, p2, p3, p4 = GetPercentages(count1, count2, count3, count4, cov)
    AmbigCombination = AmbiguityType(p1, p2, p3, p4)

    if AmbigCombination is None:
        return False, None
    if AmbigCombination == 4:
        return True, "N"
    if AmbigCombination == 3:
        if any(n == "X" for n in [nuc1, nuc2, nuc3]) is True:
            return True, "N"
        char = TripletAmbigs(nuc1, nuc2, nuc3)
        return True, char
    if AmbigCombination == 2:
        char = DoubleAmbigs(nuc1, nuc2)
        return True, char
