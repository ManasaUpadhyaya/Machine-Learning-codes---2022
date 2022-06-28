import numpy as np

min_loop_length = 4


def pair_check(tup):
    if tup in [('A', 'U'), ('U', 'A'), ('C', 'G'), ('G', 'C')]:
        return True
    return False

def matrix(seq, x, y):
    if x >= y-min_loop_length:
        return 0
    else:
        not_pair = matrix(seq, x, y-1)
        pair = [1+matrix(seq, x, a-1)+matrix(seq, a+1,y-1) for a in range(x, y-min_loop_length) if pair_check((seq[a], seq[y]))]
        if not pair:
            pair = [0]
        paired_seq = max(pair)
        return max(not_pair, paired_seq)


def cal_OPT(sequence):
    """\
    Description:
    ------------
        Implement Nussinov algorithm(dynamic programming) and construct the OPT matrix that stores the optimal score. You may use bottom up approach or top down approach.
    Parameter:
    ------------
        sequence: RNA sequence of length n, type `list`.
    Return:
    ------------
        OPT: OPT chart, an np.array of the shape (n,n)
    """
    # Implement your algorithm here
    n = len(sequence)
    OPT = np.zeros((n,n))
    for k in range(min_loop_length, n):
        for i in range(n-k):
            j = i+k
            OPT[i][j] = matrix(sequence, i, j)
    return OPT

def get_base_pairs(i, j, OPT, base_pairs, seq):
    if j<=i:
        return None#loop has run its course and we are the end of the sequence
    elif OPT[i][j] == OPT[i][j-1]:
        get_base_pairs(i, j-1, OPT, base_pairs, seq)
    else:
        for x in [a for a in range(i, j-min_loop_length) if pair_check((seq[a], seq[j]))]:
            if x-1 < 0:
                if OPT[i][j] == OPT[x+1][j-1] + 1:
                    base_pairs.append((x,j))
                    get_base_pairs(x+1, j-1, OPT, base_pairs,seq)
                    break
            elif OPT[i][j] == OPT[i][x-1] + OPT[x+1][j-1] + 1:
                base_pairs.append((x,j))
                get_base_pairs(i, x-1, OPT, base_pairs, seq)
                get_base_pairs(x+1, j-1, OPT, base_pairs, seq)
                break
    return base_pairs

def traceback(OPT, sequence):
    """\
    Description:
    ------------
        Backtracking algorithm, to find the pairing of bases in RNA sequence according to DP chart.
        You can add new utility functions to help you.    
    Parameter:
    ------------
        OPT: OPT matrix from cal_OPT(sequence)
        sequence: RNA sequence of length n, type `list`.
        fill feel to include more parameters into the function aside from `OPT` and `sequence`...
    Return:
    ------------
        structure: 
            a list of tuples that stores the pairing of the optimal solution, 
            e.g. 
            For a sequence with connection of bases 1 and 2, 3 and 4(1, 2, 3, 4 are indices, counting from 0 to n-1 for sequence of length n),
            the returned structure should be [(1, 2), (3, 4)]
    """
    # Implement your algorithm here
    i, j, base_pairs = 0, len(sequence)-1, []
    base_pairs = get_base_pairs(i, j, OPT, base_pairs, sequence) 
    return base_pairs

sequences = []
with open("test_data", "r") as fp:
    for seq in fp:
        sequences.append(seq.strip("\n"))

for i, seq in enumerate(sequences):
    OPT = cal_OPT(seq)

    structure = traceback(OPT, seq)
    print(structure)
    with open("result_" + str(i) + ".txt", "w") as fp:
        for pair in structure:
            fp.write(str(pair) + "\n")
