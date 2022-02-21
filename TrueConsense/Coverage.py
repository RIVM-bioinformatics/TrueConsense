def BuildCoverage(IndexDF, output):
    IndexDF["cov"] = IndexDF["query_sequences"].map(len)
    IndexDF[["cov"]].to_csv(output, sep="\t", header=False)


def GetCoverage(iDict, position):
    return iDict[position].get("coverage")
