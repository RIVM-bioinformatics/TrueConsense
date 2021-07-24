from .Sequences import split_to_codons

def in_orf(loc, gffd):
    exists = []
    for k in gffd.keys():
        start = gffd[k].get("start")
        stop = gffd[k].get("end")

        xx = loc in range(start, stop)
        exists.append(xx)
    
    if any(exists) == True:
        return True
    return False

def SolveTripletLength(uds, mds):
    mdslen = len(mds)
    udslen = len(uds)

    if udslen % 3 == 0:
        ## upcoming stretch is already divisible by 3
        
        if mdslen % 3 == 0:
            ## the minority-del group is also divisible by 3

            return True # return that the minority-del is good to be appended
        else:
            ## the minority-del group is not divisible by 3
            return False # return that the minority-del should be ignored
    else:
        ## the upcoming stretch is NOT divisible by 3 -> the mds should make it so
        if (mdslen + udslen) % 3 == 0 :
            ## the combination of minority-del group and the upcoming stretch is divisible by 3.
            return True
        return False
    
def CorrectGFF(oldgffdict, newgffdict, cons, p, skips, inserts):
    
    stopcodons = ["TAG", "TAA", "TGA"]
    rvstopcodon = ["CAT"]

    for k in oldgffdict.keys():
        start = oldgffdict[k].get('start')
        end = oldgffdict[k].get('end')
        orient = oldgffdict[k].get('strand')

        shift = 0
        
        if p in range(start, end):
            if orient == '+':
                rseq = ''.join(cons)[start-1:]
                shift = rseq.count('-')
                seq = rseq.replace('-', '')

                codons = split_to_codons(seq)

                it=0
                for c in codons:
                    it += 1
                    if any(s in c for s in stopcodons) is True:
                        break
                orfsize = it * 3

                newend = start + orfsize + shift

                if p == newend:
                    up = {'end': newend}
                    newgffdict[k].update(up)

    return newgffdict