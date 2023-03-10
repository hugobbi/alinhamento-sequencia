import numpy as np

def calcula_score(termo1, termo2, match, mismatch):
    if termo1 == termo2:
        return match
    else:
        return mismatch

def needleman_wunsch(seq1, seq2, match, mismatch, gap):
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

    alinhamento_1 = ""
    alinhamento_2 = ""
    i = n_seq1
    j = n_seq2
    while i > 0 and j > 0:
        score = M[i,j]
        score_diag = M[i-1, j-1]
        score_cima = M[i, j-1]
        score_esq = M[i-1, j]

        if score == score_diag + calcula_score(seq1[i-1], seq2[j-1], match, mismatch):
            alinhamento_1 += seq1[i-1]
            alinhamento_2 += seq2[j-1]
            i -= 1
            j -= 1

        elif score == score_esq + gap:
            alinhamento_1 += seq1[i-1]
            alinhamento_2 += '-'
            i -= 1
        elif score == score_cima + gap:
            alinhamento_1 += '-'
            alinhamento_2 += seq2[j-1]
            j -= 1
    
    while i > 0:
        alinhamento_1 += seq1[i-1]
        alinhamento_2 += '-'
        i -= 1

    while j > 0:
        alinhamento_1 += '-'
        alinhamento_2 += seq2[j-1]
        j -= 1
    
    alinhamento_1_final = alinhamento_1[::-1]
    alinhamento_2_final = alinhamento_2[::-1]

    print("Seq 1:")
    print(alinhamento_1_final)
    print("\nSeq 2:")
    print(alinhamento_2_final)
    print(f"\nScore = {M[n_seq1][n_seq2]}")

    count_match = 0
    for i, j in zip(alinhamento_1_final, alinhamento_2_final):
        if i == j:
            count_match += 1

    print(f"\nIdentidade = {(count_match / len(alinhamento_1_final))*100}%")

seq_carpa_comum = "ATGAGTCTCTCTGATAAGGACAAGGCTGCTGTGAAAGCCCTATGGGCTAAGATCAGCCCCAAAGCCGATGATATCGGCGCTGAAGCTCTCGGCAGAATGCTGACCGTCTACCCTCAGACCAAGACCTACTTCGCTCACTGGGATGACCTGAGCCCTGGGTCCGGTCCTGTGAAGAAGCATGGCAAGGTTATCATGGGTGCAGTGGCCGATGCCGTTTCAAAAATAGACGACCTTGTGGGAGGTCTGGCCTCCCTGAGCGAACTTCATGCTTCCAAGCTGCGTGTTGACCCGGCCAACTTCAAGATCCTCGCACACAATGTCATCGTGGTCATCGGCATGCTCTTCCCTGGAGACTTCCCCCCAGAGGTTCACATGTCAGTTGACAAGTTTTTCCAGAACTTGGCTCTGGCTCTCTCTGAGAAGTACCGCTAA"
seq_homo_sapiens = "ATGGTGCTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGGTCGGCGCGCACGCTGGCGAGTATGGTGCGGAGGCCCTGGAGAGGATGTTCCTGTCCTTCCCCACCACCAAGACCTACTTCCCGCACTTCGACCTGAGCCACGGCTCTGCCCAAGTTAAGGGCCACGGCAAGAAGGTGGCCGACGCGCTGACCAACGCCGTGGCGCACGTGGACGACATGCCCAACGCGCTGTCCGCCCTGAGCGACCTGCACGCGCACAAGCTTCGGGTGGACCCGGTCAACTTCAAGCTCCTAAGCCACTGCCTGCTGGTGACCCTGGCCGCCCACCTCCCCGCCGAGTTCACCCCTGCGGTGCACGCTTCCCTGGACAAGTTCCTGGCTTCTGTGAGCACCGTGCTGACCTCCAAATACCGTTAA"
seq_camundongo = "ATGGTGCTCTCTGGGGAAGACAAAAGCAACATCAAGGCTGCCTGGGGGAAGATTGGTGGCCATGGTGCTGAATATGGAGCTGAAGCCCTGGAAAGGATGTTTGCTAGCTTCCCCACCACCAAGACCTACTTTCCTCACTTTGATGTAAGCCACGGCTCTGCCCAGGTCAAGGGTCACGGCAAGAAGGTCGCCGATGCGCTGGCCAGTGCTGCAGGCCACCTCGATGACCTGCCCGGTGCCTTGTCTGCTCTGAGCGACCTGCATGCCCACAAGCTGCGTGTGGATCCCGTCAACTTCAAGCTCCTGAGCCACTGCCTGCTGGTGACCTTGGCTAGCCACCACCCTGCCGATTTCACCCCCGCGGTACATGCCTCTCTGGACAAATTCCTTGCCTCTGTGAGCACCGTGCTGACCTCCAAGTACCGTTAA"
seq_cabra = "ATGTCTCTGACCAGGACTGAGAGGACCATCATCCTGTCCCTGTGGAGCAAGATCTCCACACAGGCAGACGTCATTGGCACCGAGACCCTGGAGAGGCTCTTCTCCTGCTACCCGCAGGCCAAGACCTACTTCCCGCACTTCGACCTGCACTCGGGCTCCGCGCAGCTGCGCGCGCACGGCTCCAAGGTGGTGGCCGCCGTGGGCGACGCGGTCAAGAGCATCGACAACGTGACGAGCGCGCTGTCCAAGCTGAGCGAGCTGCACGCCTACGTGCTGCGCGTGGACCCGGTCAACTTCAAGTTCCTGTCCCACTGCCTGCTGGTCACGTTGGCCTCGCACTTCCCCGCCGACTTCACGGCCGACGCGCACGCCGCCTGGGACAAGTTCCTGTCCATCGTGTCCGGCGTCCTGACGGAGAAGTACCGCTGA"

GAP = -5
MATCH = 5
MISMATCH = -4
needleman_wunsch("ACGT", "AGCT", MATCH, MISMATCH, GAP)



