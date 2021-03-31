from .Inserts import ListInserts, ExtractInserts
from .Coverage import GetCoverage


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

def BuildConsensus(mincov, iDict, GffDict, bam, outdir):
    
    consensus = []
    
    UpdatedGFF = UpdateGFF(mincov, iDict, bam, GffDict)
    
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
        if loc == (last - 1):
            return loc + 1
        if loc == (last - 2):
            return loc + step
    
    lastposition = list(range(len(iDict)))[-1] + 1
    for i in range(len(iDict)):
        currentposition = i + 1
        
        prev1pos = previousposition(currentposition, 1)
        prev2pos = previousposition(currentposition, 2)
        next1pos = nextposition(currentposition, 1, lastposition)
        next2pos = nextposition(currentposition, 2, lastposition)

        within_orf = Inside_ORF(currentposition, UpdatedGFF)
        
        
        
    pass
    
def split_to_codons(seq, num):
    return [seq[start : start + num] for start in range(0, len(seq), 3)]
    
def UpdateGFF(mincov, iDict, bam, GffDict):
    draftseq, insertions = DraftConsensus(mincov, iDict, bam)
    stopcodons = ["TAG", "TAA", "TGA"]
    
    shift = 0
    for k in GffDict.keys():
        if (GffDict[k].get("strand")) is "+": ## 'forward' orientation of ORF
            start = GffDict[k].get("start")
            start = start - 1 + shift ## adjust for 0-based + adjust for nucleotide shifts due to insertions
            try:
                nextstart = GffDict[k+1].get("start")
                nextstart = nextstart - 1
            except:
                nextstart = None
                
            codons = split_to_codons(draftseq[start:], 3)
            
            i = 0
            for c in codons:
                i+=1
                if any(stops in c for stops in stopcodons) is True:
                    break
            
            orfsize = i*3
            end = start+orfsize
            
            for p in insertions.keys():
                if nextstart is not None:
                    if p in range(start, nextstart):
                        localshift = int(insertions[p])
                        #GffDict[k+1]["start"] = int(nextstart) + int(localshift) + 1
                        shift += int(localshift)
                        #print(GffDict[k+1]["start"])
            
            GffDict[k]["start"] = start + 1
            GffDict[k]["end"] = end
            
        if (GffDict[k].get("strand")) is "-": ## 'reverse' orientation of ORF
            continue
    
    return GffDict
        
def DraftConsensus(mincov, iDict, bam):
    
    draft = []
    insertions = {}
    
    hasinserts, insertpositions = ListInserts(iDict, mincov)
    for i in range(len(iDict)):
        realpos = i+1
        
        cov = GetCoverage(iDict, realpos)
        FirstNuc, SecondNuc, ThirdNuc, FourthNuc, FifthNuc = GetProminentNucleotides(iDict, realpos)
        
        if cov < mincov:
            draft.append("N")
        else:
            if FirstNuc == "D":
                draft.append('-')
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
    sorteddist = sorted(((value, key) for key, value in GetDistribution(iDict, position).items()))
    # sorteddist[-1][1] == most prominent nucleotide
    # sorteddist[-2][1] == second most prominent nucleotide
    # sorteddist[-3][1] == third most prominent nucleotide
    ## etc etc
    # sorteddist[-1][0] == counter for most prominent nucleotide
    # sorteddist[-2][0] == counter for second most prominent nucleotide
    ## etc etc
    return sorteddist[-1][1], sorteddist[-2][1], sorteddist[-3][1], sorteddist[-4][1], sorteddist[-5][1], sorteddist[-1][0]
        
def GetDistribution(iDict, position):
    dist = {}
    dist['A'] = iDict[position].get("A")
    dist['T'] = iDict[position].get("T")
    dist['C'] = iDict[position].get("C")
    dist['G'] = iDict[position].get("G")
    dist["D"] = iDict[position].get("D")
    return dist