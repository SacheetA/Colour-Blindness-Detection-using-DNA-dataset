import numpy as np


def loadLastCol(filename):
    """
    Input: Path of the file corresponding the last column (BWT).
    Output: The last column (BWT) in string format.
    """
    # function body - Begins
    LastCol = ''
    f = open(filename, 'r')
    lines = f.readlines()
    for line in lines:
        LastCol = LastCol + line.strip()            # strip() removes \n (newlines)
    # function body - Ends
    return LastCol                                  #string data type

def loadRefSeq(filename):
    """
    Input: Path of the file containing the reference sequence.
    Output: The reference sequence in string format.
    """
    # function body - Begins
    RefSeq = ''
    f = open(filename, 'r')
    lines = f.readlines()[1:]           # skips 1st line containing header
    for line in lines:
        RefSeq = RefSeq + line.strip()
    # function body - Ends
    return RefSeq                       # string data type

def loadReads(filename):
    """
    Input: Path of the file containing all the reads.
    Output: A list containing all the reads.
    """
    # function body - Begins
    Reads = []
    f = open(filename, 'r')
    lines = f.readlines()
    for line in lines:
        line = line.replace('N', 'A')               # There are 'N's in some of the reads, replace them by 'A's
        Reads.append(line.strip())                  # strip() removes \n (newlines)
    # function body - Ends
    return Reads                                    # list of strings

def loadMapToRefSeq(filename):
    """
    Input: Path of the file containing mapping of the first column to the reference sequence.
    Output: numpy integer array containing the mapping.
    """
    # function body - Begins
    MapToRefSeq = np.loadtxt(filename,dtype = np.int32)
    # function body - Ends
    return MapToRefSeq                      # numpy integer array