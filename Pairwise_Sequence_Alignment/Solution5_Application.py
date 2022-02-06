from Pairwise_Sequence_Alignment import Solution2_SubMaxtrix_GapPenalties
from Pairwise_Sequence_Alignment import Solution3_Needleman_Wunsch
from Pairwise_Sequence_Alignment import Solution4_Local_Alig_Smith_Waterman

# Find the e longest common sub-sequence between two sequences
def longest_common_subseq (sequence1, sequence2, alphabet = "ACGT"):
    submatrix = Solution2_SubMaxtrix_GapPenalties.create_submat(1,0,alphabet)
    S,T = Solution3_Needleman_Wunsch.needleman_Wunsh_Algorithsm(sequence1, sequence2, submatrix, 0)
    alingment = Solution3_Needleman_Wunsch.align_recover(T, sequence1, sequence2)

    sizeOfAling = len(alingment[0])
    largest_subsq = ""
    for i in range(sizeOfAling):
        if alingment[0][i] == alingment[1][i]:
            largest_subsq += alingment[1][i]
    return largest_subsq

# Calculation of the edit distance, the minimum number of operation
# required to transform one string to another
# input matches to have a score of 0, while gaps and mismatches are scored with −1
def edit_distance(seq1, seq2, alphabet ="ACGT"):
    submatrix = Solution2_SubMaxtrix_GapPenalties.create_submat(0,-1,alphabet)
    S,T = Solution3_Needleman_Wunsch.needleman_Wunsh_Algorithsm(seq1, seq2, submatrix, -1)
    res =-1*S[len(seq1)][len(seq2)]
    return res

# Find the longest common string in two sequences
def longest_common_string (seq1, seq2, alphabet = "ACGT"):
    m = max( len (seq1), len (seq2))
    pen = -1 * (m+1)
    sm = Solution2_SubMaxtrix_GapPenalties.create_submat(1, pen, alphabet)
    S,T,b = Solution4_Local_Alig_Smith_Waterman.smith_Water_Algorithsm(seq1, seq2, sm, pen)
    alinL= Solution4_Local_Alig_Smith_Waterman.recover_align_local(S, T, seq1, seq2)
    return alinL[0]

