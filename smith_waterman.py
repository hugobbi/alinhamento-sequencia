import numpy as np

def calcula_score(termo1, termo2, match, mismatch):
    return match if termo1 == termo2 else mismatch

def smith_waterman(seq1, seq2, match, mismatch, gap):
    n_seq1 = len(seq1)
    n_seq2 = len(seq2)
    M = np.zeros((n_seq1+1, n_seq2+1))

    max_score = 0
    max_score_positions = [(0,0)]
    possiveis_scores = np.zeros(4)
    for i in range(1, n_seq1+1):
        for j in range(1, n_seq2+1):
            possiveis_scores[0] = M[i-1, j-1] + calcula_score(seq1[i-1], seq2[j-1], match, mismatch)
            possiveis_scores[1] = M[i-1, j] + gap
            possiveis_scores[2] = M[i, j-1] + gap
            possiveis_scores[3] = 0
            M[i,j] = np.max(possiveis_scores)
            if M[i,j] > max_score:
                max_score = M[i,j]
                max_score_positions = [(i,j)]
            elif M[i,j] == max_score:
                max_score_positions.append((i,j))
    print(M)
    print(f"{max_score=}")
    print(f"{max_score_positions=}")

    for max_pos in max_score_positions:    
        al1 = ""
        al2 = ""
        i,j = max_pos
        while M[i,j] != 0:
            score = M[i,j]
            score_diag = M[i-1, j-1]
            score_cima = M[i, j-1]
            score_esq = M[i-1, j]
            if score == score_diag + calcula_score(seq1[i-1], seq2[j-1], match, mismatch):
                al1 += seq1[i-1]
                al2 += seq2[j-1]
                i -= 1
                j -= 1
            elif score == score_esq + gap:
                al1 += seq1[i-1]
                al2 += '-'
                i -= 1
            elif score == score_cima + gap:
                al1 += '-'
                al2 += seq2[j-1]
                j -= 1

        al1 = al1[::-1]
        al2 = al2[::-1]
            
        print(f"Al1: \n{al1}")
        print(f"\nAl2: \n{al2}")
        print(f"\n{max_score=}")
        match_count = 0
        for i,j in zip(al1, al2):
            if i == j:
                match_count += 1
        print(f"Identidade: {match_count/len(al1)*100}%")

CHIMPANZEE = "MTENSTSAPAAKPKRAKASKKSTDHPKYSDMIVAAIQAEKNRAGSSRQSIQKYIKSHYKVGENADSQIKLSIKRLVTTGVLKQTKGVGASGSFRLAKSDEPKKSVAFKKTKKEIKKVATPKKASKPKKAASKAPTKKPKATPVKKAKKKLAATPKKAKKPKTVKAKPVKASKPKKAKPVKPKAKSSAKRAGKKK"
COW = "MTENSTSTPAAKPKRAKASKKSTDHPKYSDMIVAAIQAEKNRAGSSRQSIQKYIKSHYKVGENADSQIKLSIKRLVTTGVLKQTKGVGASGSFRLAKSDEPKRSVAFKKTKKEVKKVATPKKAAKPKKAASKAPSKKPKATPVKKAKKKPAATPKKTKKPKTVKAKPVKASKPKKTKPVKPKAKSSAKRTGKKK"

GAP = -2
MATCH = 1
MISMATCH = -1
smith_waterman(CHIMPANZEE, COW, MATCH, MISMATCH, GAP)