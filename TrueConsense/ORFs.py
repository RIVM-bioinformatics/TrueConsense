
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

def split_to_codons(seq):
    return [seq[start : start + 3] for start in range(0, len(seq), 3)]

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
    
def CorrectStartPositions(gffd, shifts, p):
    for k in gffd.keys():
        start = gffd[k].get('start')
        
        if start > p:
            nstart = int(start) + int(shifts)
            update = {'start': nstart}
            
            gffd[k].update(update)
    return gffd
    
def CorrectGFF(oldgffdict, newgffdict, cons, p, inserts):
    
    stopcodons = ["TAG", "TAA", "TGA"]
    #rvstopcodon = ["CAT"]
    
    if inserts is not None and p in inserts:
        newgffdict = CorrectStartPositions(newgffdict, list(inserts[p].keys())[0], p)
        

    for k in newgffdict.keys():
        start = newgffdict[k].get('start')
        end = newgffdict[k].get('end')
        orient = newgffdict[k].get('strand')

        shift = 0
        
        if p in range(start, end):
            if orient == '+':
                rseq = ''.join(cons)[start-1:]
                shift = rseq.count('-')
                seq = rseq.replace('-', '')

                if cons[-1] == '-':
                    achieved = False
                    override_end = oldgffdict[k].get('end')
                    up = {'end': override_end}
                    newgffdict[k].update(up)
                    
                else:
                    
                    codons = split_to_codons(seq)

                    achieved = False

                    it=0
                    for c in codons:
                        it += 1
                        if any(s in c for s in stopcodons) is True:
                            achieved = True
                            break
                    orfsize = it * 3

                    if shift % 3 == 0:
                        newend = start + orfsize + shift
                    else:
                        newend = start + orfsize + shift-1
                    
                    
                    if achieved is False:
                        newend = newend + 1

                    if p == newend:
                        up = {'end': newend}
                        newgffdict[k].update(up)

    return newgffdict