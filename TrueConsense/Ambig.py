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
    p1 = (c1 / cov) * 100
    p2 = (c2 / cov) * 100
    p3 = (c3 / cov) * 100
    p4 = (c4 / cov) * 100
    return p1, p2, p3, p4


def AmbiguityType(p1, p2, p3, p4):
    maxdistance = 3
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


def IsAmbiguous(nuc1, nuc2, nuc3, nuc4, count1, count2, count3, count4, cov):
    if cov == 0:
        return False, None
    if nuc1 == "X" or nuc2 == "X":
        return False, None
    p1, p2, p3, p4 = GetPercentages(count1, count2, count3, count4, cov)
    AmbigCombination = AmbiguityType(p1, p2, p3, p4)

    if AmbigCombination is None:
        return False, None
    else:
        if AmbigCombination == 4:
            return True, "N"
        if AmbigCombination == 3:
            if any(n == "X" for n in [nuc1, nuc2, nuc3]) is True:
                return True, "N"
            else:
                char = TripletAmbigs(nuc1, nuc2, nuc3)
                return True, char
        if AmbigCombination == 2:
            char = DoubleAmbigs(nuc1, nuc2)
            return True, char
