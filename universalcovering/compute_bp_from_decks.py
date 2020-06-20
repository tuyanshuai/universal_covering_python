import numpy as np
def compute_bp_from_decks(decks,mc,p):
    bp = np.zeros(len(decks), dtype=np.complex)
    bp[0] = p
    j = 0
    for i in range(1,len(decks)):
        bp[mc[j]+1]  =  decks[j](bp[j])
        j = mc[j]+1
    return bp
