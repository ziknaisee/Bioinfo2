import numpy as np
import pandas as pd
import streamlit as st

# Function to initialize the matrix
def initialize_matrix(seq1, seq2, gap_penalty=-2):
    rows = len(seq1) + 1
    cols = len(seq2) + 1
    matrix = np.zeros((rows, cols), dtype=int)

    for i in range(1, rows):
        matrix[i][0] = gap_penalty * i
    for j in range(1, cols):
        matrix[0][j] = gap_penalty * j

    return matrix

# Function to fill the matrix for Needleman-Wunsch
def fill_matrix_needleman_wunsch(seq1, seq2, matrix, match_score=1, mismatch_score=-1, gap_penalty=-2):
    rows, cols = matrix.shape

    for i in range(1, rows):
        for j in range(1, cols):
            match = matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score)
            delete = matrix[i-1][j] + gap_penalty
            insert = matrix[i][j-1] + gap_penalty
            matrix[i][j] = max(match, delete, insert)

    return matrix

# Traceback function for Needleman-Wunsch
def traceback_needleman_wunsch(seq1, seq2, matrix, match_score=1, mismatch_score=-1, gap_penalty=-2):
    alignment_a, alignment_b = "", ""
    traceback_path = []
    i, j = len(seq1), len(seq2)

    while i > 0 and j > 0:
        traceback_path.append((i, j))
        score_current = matrix[i][j]
        score_diagonal = matrix[i-1][j-1]
        score_up = matrix[i-1][j]
        score_left = matrix[i][j-1]

        if score_current == score_diagonal + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score):
            alignment_a = seq1[i-1] + alignment_a
            alignment_b = seq2[j-1] + alignment_b
            i -= 1
            j -= 1
        elif score_current == score_up + gap_penalty:
            alignment_a = seq1[i-1] + alignment_a
            alignment_b = "-" + alignment_b
            i -= 1
        else:
            alignment_a = "-" + alignment_a
            alignment_b = seq2[j-1] + alignment_b
            j -= 1

    while i > 0:
        traceback_path.append((i, 0))
        alignment_a = seq1[i-1] + alignment_a
        alignment_b = "-" + alignment_b
        i -= 1

    while j > 0:
        traceback_path.append((0, j))
        alignment_a = "-" + alignment_a
        alignment_b = seq2[j-1] + alignment_b
        j -= 1

    traceback_path.reverse()
    return alignment_a, alignment_b, traceback_path

# Highlight function for traceback path
def highlight_traceback(data, traceback_path):
    def style_cell(row, col):
        return 'background-color: lightgreen;' if (row, col) in traceback_path else ''

    # Create a DataFrame of styles
    styles = pd.DataFrame(
        [[style_cell(row, col) for col in range(data.shape[1])] for row in range(data.shape[0])],
        index=data.index,
        columns=data.columns
    )
    return styles

# Streamlit app
st.title("Pairwise Sequence Alignment with Traceback")
st.write("Compare two sequences using Needleman-Wunsch algorithm and visualize the traceback.")

# Input sequences
seq1 = st.text_input("Enter Sequence 1").upper()
seq2 = st.text_input("Enter Sequence 2").upper()

# Parameters
match_score = st.number_input("Match score", value=1, step=1)
mismatch_score = st.number_input("Mismatch score", value=-1, step=1)
gap_penalty = st.number_input("Gap penalty", value=-2, step=1)

if st.button("Align Sequences"):
    if seq1 and seq2:
        # Initialize the matrix
        matrix = initialize_matrix(seq1, seq2, gap_penalty)

        # Fill the matrix
        matrix = fill_matrix_needleman_wunsch(seq1, seq2, matrix, match_score, mismatch_score, gap_penalty)

        # Perform traceback
        alignment_a, alignment_b, traceback_path = traceback_needleman_wunsch(
            seq1, seq2, matrix, match_score, mismatch_score, gap_penalty
        )

        # Display the alignment matrix
        st.write("Alignment Matrix:")
        df_matrix = pd.DataFrame(matrix, index=["-"] + list(seq1), columns=["-"] + list(seq2))

        # Highlight the traceback
        styled_df = df_matrix.style.apply(
            lambda x: highlight_traceback(df_matrix, traceback_path),
            axis=None
        )

        # Display styled DataFrame
        st.dataframe(styled_df)

        # Display optimal alignment
        st.write("Optimal Alignment:")
        st.text(f"Sequence 1: {alignment_a}")
        st.text(f"Sequence 2: {alignment_b}")

        # Traceback path
        st.write("Traceback Path (Row, Column):")
        st.write(traceback_path)

    else:
        st.warning("Please enter both sequences.")
