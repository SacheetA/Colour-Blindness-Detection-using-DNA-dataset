# import python libraries
import numpy as np
from scipy.special import comb  # useful to compute probabilities
from utils import *

"""
Starting and ending locations (indices) of red and green exons in the reference sequence - Begins

1. Red Exon Locations
"""
RedExonPos = np.array([
    [149249757, 149249868], # R1
    [149256127, 149256423], # R2
    [149258412, 149258580], # R3
    [149260048, 149260213], # R4
    [149261768, 149262007], # R5
    [149264290, 149264400]  # R6
    ])
"""
2. Green Exon Locations
"""
GreenExonPos = np.array([
    [149288166, 149288277], # G1
    [149293258, 149293554], # G2
    [149295542, 149295710], # G3
    [149297178, 149297343], # G4
    [149298898, 149299137], # G5
    [149301420, 149301530]  # G6
    ])
"""
Starting and ending locations (indices) of red and green exons in the reference sequence - Ends
"""

class Query:  
    def __init__(self, LastCol):
        self.delta = 100
        self.lc = LastCol
        self.rc = self.__reduced_lastCol()      

    def __count(self, index = int, delta = int, prev_count = int):
        count = prev_count
        start = (index - 1) * delta
        end = min((index * delta), len(self.lc))     
        for element in self.lc[start:end]:
            if element == 'A':
                count[0] += 1
            elif element == 'T':
                count[1] += 1
            elif element == 'G':
                count[2] += 1
            elif element == 'C':
                count[3] += 1   
        return count

    def __reduced_lastCol(self):
        dim = (len(self.lc) // self.delta) + 1
        dict = {}
        dict['A'] = [0]
        dict['T'] = [0]
        dict['G'] = [0]
        dict['C'] = [0]
        for i in range(dim):
            latest_count = [dict['A'][i], dict['T'][i], dict['G'][i], dict['C'][i]]
            cumulative_count = self.__count(index = i+1 , delta = self.delta, prev_count = latest_count)
            dict['A'].append(cumulative_count[0])
            dict['T'].append(cumulative_count[1])
            dict['G'].append(cumulative_count[2])
            dict['C'].append(cumulative_count[3])
        # for item in dict.values():
        #     item.remove(item[0])
        return dict

    def rank(self, char = str, index = int):
        interval_lower = index // self.delta
        count = self.rc[char][interval_lower]
        for element in self.lc[(interval_lower*self.delta) : index]:
            if element == char:
                count += 1
        return count


def matchString(read):
    read = read[::-1]
    p = query.rank('A', len(LastCol))
    q = query.rank('C', len(LastCol))
    r = query.rank('G', len(LastCol))
    idx1 = 0
    idx2 = len(LastCol)
    k = 0
    mismatch = False
    output = []
    
    for ele in read:                         # ele --> element
        l = query.rank(ele, idx1)            # lower index
        u = query.rank(ele, idx2)            # upper index  

        # Discontinue if there is mismatch
        if l == u:
            mismatch = True
            break

        # Compute the indices in First Column
        if ele == 'A':
            idx1 = l
            idx2 = u
        elif ele == 'C':
            idx1 = p + l           # no. of A's + first index of C in Last Column of BWT
            idx2 = p + u           # no. of A's + last index of C in Last Column of BWT, similarly follows for G and T
        elif ele == 'G': 
            idx1 = (p+q) + l
            idx2 = (p+q) + u
        elif ele == 'T':
            idx1 = (p+q+r) + l
            idx2 = (p+q+r) + u

    if not(mismatch):
        output = Map[idx1 : idx2]

    return output, mismatch


def ReverseComplement(read):
    new_read = ''
    for r in read:
        if r == 'A':
            new_read += 'T'
        elif r == 'C':
            new_read += 'G'
        elif r == 'G':
            new_read += 'C'
        elif r == 'T':
            new_read += 'A'
    return new_read


def traverse_match(read, op, sub_len):
    i = 0
    possible_locs_idx =  []
    possible_locs = []
    for ele in op:
        if len(ele) != 0:
            possible_locs_idx.append(i)
        i += 1

    for pl in possible_locs_idx:
        locs = op[pl]
        if pl == 0:
            for loc in locs:
                k = 0
                if (len(read) + loc) > len(RefSeq):       # handles index error towards the end of RefSeq
                    continue
                for i in range(len(read)):
                    if read[i] != RefSeq[loc + i]:
                        k += 1
                if k <= 2:
                    possible_locs.append(loc)
        elif pl == 1:
            for loc in locs:
                k = 0
                if (len(read) + loc - sub_len) > len(RefSeq):       # handles index error towards the end of RefSeq
                    continue
                for i in range(len(read)):
                    if read[i] != RefSeq[loc- sub_len + i]:
                        k += 1
                if k <= 2:
                    possible_locs.append(loc - sub_len)
        elif pl == 2:
            for loc in locs:
                k = 0
                if (len(read) + loc - 2*sub_len) > len(RefSeq):       # handles index error towards the end of RefSeq
                    continue
                for i in range(len(read)):
                    if read[i] != RefSeq[loc- 2*sub_len + i]:
                        k += 1
                if k <= 2:
                    possible_locs.append(loc - 2 * sub_len)

    return [*set(possible_locs)]        #removes duplicate elements from possible_locs


def match_single_read(read):
    # Initializing
    sub_len = len(read) // 3
    read1 = read[:sub_len]
    read2 = read[sub_len: 2*sub_len]
    read3 = read[2*sub_len:]
    rev_complement = False
    op = [0, 0, 0]
    loc = []

    # Read Matching
    op[0], mismatch1 = matchString(read1)
    op[1], mismatch2 = matchString(read2)
    op[2], mismatch3 = matchString(read3)

    if (mismatch1 + mismatch2 + mismatch3) > 2:
        rev_complement = True
    else:
        loc = traverse_match(read, op, sub_len)
        if len(loc) == 0:
            rev_complement = True
        
    return loc, rev_complement


def convert_to_prob(configurations):
    conf = []
    for configuration in configurations:               
        c = []
        for ele in configuration:
            ele /= 100                          #convert from percetange to fraction
            c.append(ele /(1 + ele))
        conf.append(c)
    return conf[0], conf[1], conf[2], conf[3]


def MatchReadToLoc(read):
    """
    Input: a read (string)
    Output: list of potential locations at which the read may match the reference sequence. 
    Refer to example in the README file.
    IMPORTANT: This function must also handle the following:
        1. cases where the read does not match with the reference sequence
        2. any other special case that you may encounter
    """
    # function body - Begins
    positions, rev_comp = match_single_read(read)
    if rev_comp:
        reverse_complement = ReverseComplement(read)
        positions, _ = match_single_read(reverse_complement)       
    # function body - Ends
    return positions # list of potential locations at which the read may match the reference sequence.



def WhichExon(positions):
    """
    Input: list of potential locations at which the read may match the reference sequence.
    Output: Update(increment) to the counts of the 12 exons
    IMPORTANT: This function must also handle the following:
        1. cases where the read does not match with the reference sequence
        2. cases where there are more than one matches (match at two exons)
        3. any other special case that you may encounter
    """
    r1,r2,r3,r4,r5,r6,g1,g2,g3,g4,g5,g6 = 0,0,0,0,0,0,0,0,0,0,0,0

    # function body - Begins
    r = read
    l = len(read)
    r = np.zeros(len(RedExonPos), dtype = np.float16)
    g = np.zeros(len(GreenExonPos), dtype = np.float16)
    red = np.zeros(len(RedExonPos), dtype = np.float16)
    green = np.zeros(len(GreenExonPos), dtype = np.float16)

    for position in positions:
        sub_parts = [position, position + len(read) - 1]
        idx = []
        idx1 = []
        for sp in sub_parts:
            for i in range(len(RedExonPos)):
                if RedExonPos[i, 0] <= sp <= RedExonPos[i, 1]:
                    r[i] = 1
                    idx.append(i) 
            if any(r == 1):
                red += r
                r[idx] = 0
                break  
                  
            for i in range(len(GreenExonPos)):
                if GreenExonPos[i, 0] <= sp <= GreenExonPos[i, 1]:
                    g[i] = 1
                    idx1.append(i)
            if any(g == 1):
                green += g
                g[idx1] = 0
                break    

    #Exon1
    if red[0] == 1 and green[0] == 1:
        r1 += 0.5
        g1 += 0.5
    elif red[0] == 1:
        r1 += 1
    elif green[0] == 1:
        g1 += 1

    #Exon2
    if red[1] == 1 and green[1] == 1:
        r2 += 0.5
        g2 += 0.5
    elif red[1] == 1:
        r2 += 1
    elif green[1] == 1:
        g2 += 1

    #Exon3
    if red[2] == 1 and green[2] == 1:
        r3 += 0.5
        g3 += 0.5
    elif red[2] == 1:
        r3 += 1
    elif green[2] == 1:
        g3 += 1

    #Exon4
    if red[3] == 1 and green[3] == 1:
        r4 += 0.5
        g4 += 0.5
    elif red[3] == 1:
        r4 += 1
    elif green[3] == 1:
        g4 += 1

    #Exon5
    if red[4] == 1 and green[4] == 1:
        r5 += 0.5
        g5 += 0.5
    elif red[4] == 1:
        r5 += 1
    elif green[4] == 1:
        g5 += 1

    #Exon6
    if red[5] == 1 and green[5] == 1:
        r6 += 0.5
        g6 += 0.5
    elif red[5] == 1:
        r6 += 1
    elif green[5] == 1:
        g6 += 1    
    # function body - Ends     
    return np.array([r1,r2,r3,r4,r5,r6,g1,g2,g3,g4,g5,g6])



def ComputeProb(ExonMatchCounts):
    """
    Input: The counts for each exon
    Output: Probabilities of each of the four configurations (a list of four real numbers)
    """
    # function body - Begins
    # Four Configurations:
    C0 = [50, 50, 50, 50]
    C1 = [100, 100, 0, 0]
    C2 = [33, 33, 100, 100]
    C3 = [33, 33, 33, 100] 
    C0, C1, C2, C3 = convert_to_prob([C0, C1, C2, C3])

    # counts of exon r2, r3, r4, r5, g2, g3, g4, g5
    r2 = ExonMatchCounts[1]
    r3 = ExonMatchCounts[2]
    r4 = ExonMatchCounts[3]
    r5 = ExonMatchCounts[4]
    g2 = ExonMatchCounts[7]
    g3 = ExonMatchCounts[8]
    g4 = ExonMatchCounts[9]
    g5 = ExonMatchCounts[10]

    # Computing Probabilities
    # Probability of generating the counts given configuration C0   
    # used scipy.special.comb() to compute combinations 
    P01 = (comb(r2+g2, r2) * (C0[0]**r2) * ((1-C0[0])**g2))
    P02 = (comb(r3+g3, r3) * (C0[1]**r3) * ((1-C0[1])**g3)) 
    P03 = (comb(r4+g4, r4) * (C0[2]**r4) * ((1-C0[2])**g4))
    P04 = (comb(r5+g5, r5) * (C0[3]**r5) * ((1-C0[3])**g5))  
    P0 = P01 * P02 * P03 * P04

    # Probability of generating the counts given configuration C1
    # used scipy.special.comb() to compute combinations 
    P11 = (comb(r2+g2, r2) * (C1[0]**r2) * ((1-C1[0])**g2))
    P12 = (comb(r3+g3, r3) * (C1[1]**r3) * ((1-C1[1])**g3)) 
    P13 = (comb(r4+g4, r4) * (C1[2]**r4) * ((1-C1[2])**g4))
    P14 = (comb(r5+g5, r5) * (C1[3]**r5) * ((1-C1[3])**g5))  
    P1 = P11 * P12 * P13 * P14

    # Probability of generating the counts given configuration C2   
    # used scipy.special.comb() to compute combinations 
    P21 = (comb(r2+g2, r2) * (C2[0]**r2) * ((1-C2[0])**g2))
    P22 = (comb(r3+g3, r3) * (C2[1]**r3) * ((1-C2[1])**g3)) 
    P23 = (comb(r4+g4, r4) * (C2[2]**r4) * ((1-C2[2])**g4))
    P24 = (comb(r5+g5, r5) * (C2[3]**r5) * ((1-C2[3])**g5))  
    P2 = P21 * P22 * P23 * P24

    # Probability of generating the counts given configuration C3   
    # used scipy.special.comb() to compute combinations 
    P31 = (comb(r2+g2, r2) * (C3[0]**r2) * ((1-C3[0])**g2))
    P32 = (comb(r3+g3, r3) * (C3[1]**r3) * ((1-C3[1])**g3)) 
    P33 = (comb(r4+g4, r4) * (C3[2]**r4) * ((1-C3[2])**g4))
    P34 = (comb(r5+g5, r5) * (C3[3]**r5) * ((1-C3[3])**g5))  
    P3 = P31 * P32 * P33 * P34

    # function body - ends
    return [P0, P1, P2, P3]



def BestMatch(ListProb):
    """
    Input: Probabilities of each of the four configurations (a list of four real numbers)
    Output: Most likely configuration (an integer). Refer to lecture slides
    """
    # function body - Begins
    i = 0
    temp = - np.inf
    for ele in ListProb:
        if ele > temp:
            temp = ele
            i += 1
    MostLikelyConfiguration = i
    # function body - ends
    return MostLikelyConfiguration # it holds 0, 1, 2, or 3



if __name__ == "__main__":
    # load all the data files
    LastCol = loadLastCol("../data/chrX_last_col.txt") # loads the last column
    RefSeq = loadRefSeq("../data/chrX.fa") # loads the reference sequence
    Reads = loadReads("../data/reads") # loads the reads
    Map = loadMapToRefSeq("../data/chrX_map.txt") # loads the mapping to the reference sequence
    query = Query(LastCol)      # Utility object of Query class

    # run the functions
    ExonMatchCounts = np.zeros(12) # initialize the counts for exons
    for read in Reads: # update the counts for exons
        positions = MatchReadToLoc(read) # get the list of potential match locations
        ExonMatchCounts += WhichExon(positions) # update the counts of exons, if applicable
    
    ListProb = ComputeProb(ExonMatchCounts) # compute probabilities of each of the four configurations
    MostLikely = BestMatch(ListProb) # find the most likely configuration
    print("Configuration %d is the best match"%MostLikely)
