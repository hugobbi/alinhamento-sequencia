import numpy as np

def calcula_score(termo1, termo2, match, mismatch):
    if termo1 == termo2:
        return match
    else:
        return mismatch

def smith_waterman(seq1, seq2, match, mismatch, gap):
    n_seq1 = len(seq1)
    n_seq2 = len(seq2)

    M = np.zeros((n_seq1+1, n_seq2+1))
    M[:,0] = np.linspace(0, gap*n_seq1, n_seq1+1)
    M[0,:] = np.linspace(0, gap*n_seq2, n_seq2+1)

    possiveis_scores = np.zeros(3)
    for i in range(1, n_seq1+1):
        for j in range(1, n_seq2+1):
            possiveis_scores[0] = M[i-1, j-1] + calcula_score(seq1[i-1], seq2[j-1], match, mismatch)
            possiveis_scores[1] = M[i-1, j] + gap
            possiveis_scores[2] = M[i, j-1] + gap
            M[i,j] = np.max(possiveis_scores)

    print(M)