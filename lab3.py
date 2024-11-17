import streamlit as st
import numpy as np

# Function to calculate the alignment score and matrix for global alignment
def needleman_wunsch(seq1, seq2, match_score=1, mismatch_penalty=-1, gap_penalty=2):
    n, m = len(seq1), len(seq2)
    matrix = np.zeros((n + 1, m + 1), dtype=int)
    backtrace = np.zeros((n + 1, m + 1), dtype=tuple)

    # Initialize the scoring matrix and backtrace pointers
    for i in range(n + 1):
        matrix[i][0] = i * gap_penalty
    for j in range(m + 1):
        matrix[0][j] = j * gap_penalty

    # Fill the matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty)
            delete = matrix[i - 1][j] + gap_penalty
            insert = matrix[i][j - 1] + gap_penalty
            matrix[i][j] = max(match, delete, insert)

            if matrix[i][j] == match:
                backtrace[i][j] = (i - 1, j - 1)
            elif matrix[i][j] == delete:
                backtrace[i][j] = (i - 1, j)
            else:
                backtrace[i][j] = (i, j - 1)

    # Traceback
    aligned_seq1, aligned_seq2 = "", ""
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and backtrace[i][j] == (i - 1, j - 1):
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            i, j = i - 1, j - 1
        elif i > 0 and backtrace[i][j] == (i - 1, j):
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            j -= 1

    return matrix, aligned_seq1, aligned_seq2

# Function to calculate the alignment score and matrix for local alignment
def smith_waterman(seq1, seq2, match_score, mismatch_penalty, gap_penalty):
    n, m = len(seq1), len(seq2)
    matrix = np.zeros((n + 1, m + 1), dtype=int)
    backtrace = np.zeros((n + 1, m + 1), dtype=tuple)

    # Fill the matrix
    max_score = 0
    max_pos = None
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty)
            delete = matrix[i - 1][j] + gap_penalty
            insert = matrix[i][j - 1] + gap_penalty
            matrix[i][j] = max(0, match, delete, insert)

            if matrix[i][j] == match:
                backtrace[i][j] = (i - 1, j - 1)
            elif matrix[i][j] == delete:
                backtrace[i][j] = (i - 1, j)
            elif matrix[i][j] == insert:
                backtrace[i][j] = (i, j - 1)

            if matrix[i][j] > max_score:
                max_score = matrix[i][j]
                max_pos = (i, j)

    # Traceback
    aligned_seq1, aligned_seq2 = "", ""
    i, j = max_pos
    while matrix[i][j] > 0:
        if backtrace[i][j] == (i - 1, j - 1):
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            i, j = i - 1, j - 1
        elif backtrace[i][j] == (i - 1, j):
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            j -= 1

    return matrix, aligned_seq1, aligned_seq2

# Streamlit application
st.title("Sequence Alignment")
st.sidebar.header("Alignment Parameters")
seq1 = st.text_input("Enter the first sequence:", "ACGTG")
seq2 = st.text_input("Enter the second sequence:", "ACTG")
match_score = st.sidebar.number_input("Match score:", value=1, step=1)
mismatch_penalty = st.sidebar.number_input("Mismatch penalty:", value=-1, step=1)
gap_penalty = st.sidebar.number_input("Gap penalty:", value=-2, step=1)

alignment_algorithm = st.selectbox("Choose the alignment algorithm:", ("Needleman-Wunsch (Global)", "Smith-Waterman (Local)"))

if st.button("Run Alignment"):
    if alignment_algorithm == "Needleman-Wunsch (Global)":
        matrix, aligned_seq1, aligned_seq2 = needleman_wunsch(seq1, seq2, match_score, mismatch_penalty, gap_penalty)
        alignment_type = "Global Alignment (Needleman-Wunsch)"
    else:
        matrix, aligned_seq1, aligned_seq2 = smith_waterman(seq1, seq2, match_score, mismatch_penalty, gap_penalty)
        alignment_type = "Local Alignment (Smith-Waterman)"

    st.subheader("Alignment Results")
    st.write(f"**{alignment_type}**")
    st.write("**Aligned Sequences:**")
    st.text(aligned_seq1)
    st.text(aligned_seq2)
    score = matrix[len(seq1)][len(seq2)] if alignment_algorithm == "Needleman-Wunsch (Global)" else np.max(matrix)
    st.write(f"**Alignment Score:** {score}")

    st.write("**Scoring Matrix:**")
    st.dataframe(matrix)
