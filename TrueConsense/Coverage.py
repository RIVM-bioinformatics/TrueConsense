

def BuildCoverage(iDict, output):
    with open(output, "w") as outfile:
        for i in range(len(iDict)):
            cov = iDict[i+1].get('coverage')
            outfile.write(str(i+1) + "\t" + str(cov) + "\n")
    pass

def GetCoverage(iDict, position):
    return iDict[position].get("coverage")