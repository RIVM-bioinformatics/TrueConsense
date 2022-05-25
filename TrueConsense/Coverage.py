def GetCoverage(iDict, position):
    """Takes a dictionary and a position as input, and returns the coverage of that position

    Parameters
    ----------
    iDict
        A dictionary with pileup contents per position
    position
        the position in the genome that you want to get the coverage for

    Returns
    -------
        The coverage of the position.

    """
    return iDict[position].get("coverage")
