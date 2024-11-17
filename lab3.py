import streamlit as st
import numpy as np
import pandas as pd

# Function to calculate the alignment score and matrix for global alignment (Needleman-Wunsch)
def needleman_wunsch(seq1, seq2, match_score, mismatch_penalty, gap_penalty):
    n, m = len(seq1), len(seq2)
    matrix = np.zeros((n + 1, m + 1), dtype=int)
    backtrace = np.zeros((n + 1, m + 1), dtype=tuple)

    # Initialize the scoring matrix
    for i in range(n + 1):
        matrix[i][0] = i * gap_penalty
    for j in range(m + 1):
        matrix[0][j] = j * gap_penalty

    # Fill the scoring matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty)
            delete = matrix[i - 1][j] + gap_penalty
            insert = matrix[i][j - 1] + gap_penalty
            matrix[i][j] = max(match, delete, insert)

            # Record backtrace path
            if matrix[i][j] == match:
                backtrace[i][j] = (i - 1, j - 1)
            elif matrix[i][j] == delete:
                backtrace[i][j] = (i - 1, j)
            else:
                backtrace[i][j] = (i, j - 1)

    return matrix, backtrace

# Function to calculate the alignment score and matrix for local alignment (Smith-Waterman)
def smith_waterman(seq1, seq2, match_score, mismatch_penalty, gap_penalty):
    n, m = len(seq1), len(seq2)
    matrix = np.zeros((n + 1, m + 1), dtype=int)
    backtrace = np.zeros((n + 1, m + 1), dtype=tuple)

    max_score = 0
    max_pos = None

    # Fill the scoring matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty)
            delete = matrix[i - 1][j] + gap_penalty
            insert = matrix[i][j - 1] + gap_penalty
            matrix[i][j] = max(0, match, delete, insert)

            # Record backtrace path
            if matrix[i][j] == match:
                backtrace[i][j] = (i - 1, j - 1)
            elif matrix[i][j] == delete:
                backtrace[i][j] = (i - 1, j)
            elif matrix[i][j] == insert:
                backtrace[i][j] = (i, j - 1)

            if matrix[i][j] > max_score:
                max_score = matrix[i][j]
                max_pos = (i, j)

    return matrix, backtrace, max_pos

# Function to perform traceback
def traceback(matrix, backtrace, seq1, seq2, start_pos=None):
    aligned_seq1, aligned_seq2 = "", ""
    traceback_path = []
    
    if start_pos:
        i, j = start_pos
    else:
        i, j = len(seq1), len(seq2)

    while (i > 0 or j > 0) and (start_pos is None or matrix[i][j] > 0):
        traceback_path.append((i, j))
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

    return aligned_seq1, aligned_seq2, traceback_path

# Streamlit app
st.title("Sequence Alignment Tool")
st.sidebar.header("Alignment Settings")
alignment_type = st.sidebar.selectbox("Select Alignment Type", ["Needleman-Wunsch (Global)", "Smith-Waterman (Local)"])
match_score = st.sidebar.number_input("Match Score", value=1, step=1)
mismatch_penalty = st.sidebar.number_input("Mismatch Penalty", value=-1, step=1)
gap_penalty = st.sidebar.number_input("Gap Penalty", value=-2, step=1)

seq1 = st.text_input("Enter First Sequence", "ACGTG")
seq2 = st.text_input("Enter Second Sequence", "ACTG")

if st.button("Run Alignment"):
    if alignment_type == "Needleman-Wunsch (Global)":
        matrix, backtrace = needleman_wunsch(seq1, seq2, match_score, mismatch_penalty, gap_penalty)
        aligned_seq1, aligned_seq2, traceback_path = traceback(matrix, backtrace, seq1, seq2)
    else:
        matrix, backtrace, max_pos = smith_waterman(seq1, seq2, match_score, mismatch_penalty, gap_penalty)
        aligned_seq1, aligned_seq2, traceback_path = traceback(matrix, backtrace, seq1, seq2, max_pos)

    st.subheader("Results")
    st.write("**Aligned Sequences:**")
    st.text(aligned_seq1)
    st.text(aligned_seq2)
    st.write(f"**Alignment Score:** {matrix[len(seq1)][len(seq2)] if alignment_type == 'Needleman-Wunsch (Global)' else np.max(matrix)}")

    st.write("**Scoring Matrix with Traceback Path:**")
    df_matrix = pd.DataFrame(matrix, index=["-"] + list(seq1), columns=["-"] + list(seq2))
    styled_matrix = df_matrix.style.apply(
        lambda x: ["background-color: yellow" if (i, j) in traceback_path else "" for i, j in enumerate(x.index)],
        axis=1,
    )
    st.dataframe(styled_matrix)
