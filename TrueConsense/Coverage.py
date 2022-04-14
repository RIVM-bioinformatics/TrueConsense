def BuildCoverage(iDict, output):
    """This function takes a dictionary of dictionaries and an output file name as input, and writes the
    coverage of each position to the output file.

    Parameters
    ----------
    iDict
        a dictionary of with pileup contents per position
    output
        the name of the output file

    """
    with open(output, "w") as outfile:
        for i in range(len(iDict)):
            cov = iDict[i + 1].get("coverage")
            outfile.write(str(i + 1) + "\t" + str(cov) + "\n")


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
