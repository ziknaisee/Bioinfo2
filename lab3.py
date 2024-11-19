import streamlit as st
import numpy as np

# Function to calculate the alignment score and matrix for global alignment (Needleman-Wunsch)
def needleman_wunsch(seq1, seq2, match_score, mismatch_penalty, gap_penalty):
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
    path = []
    while i > 0 or j > 0:
        path.append((i, j))
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

    return matrix, aligned_seq1, aligned_seq2, path

# Streamlit application
st.title("Sequence Alignment Visualization")
st.sidebar.header("Alignment Parameters")

# Input parameters
seq1 = st.text_input("Enter the first sequence (Sequence 1):", "FKHMEDPLE")
seq2 = st.text_input("Enter the second sequence (Sequence 2):", "FMDTPLNE")
match_score = st.sidebar.number_input("Match score:", value=1, step=1)
mismatch_penalty = st.sidebar.number_input("Mismatch penalty:", value=-1, step=1)
gap_penalty = st.sidebar.number_input("Gap penalty:", value=-2, step=1)

if st.button("Run Needleman-Wunsch Alignment"):
    matrix, aligned_seq1, aligned_seq2, path = needleman_wunsch(seq1, seq2, match_score, mismatch_penalty, gap_penalty)

    # Display the scoring matrix
    st.subheader("Scoring Matrix")
    st.write("Sequences displayed outside the table")

    # HTML table for visualization
    html = "<table style='border-collapse: collapse; text-align: center;'>"
    html += "<tr><th></th><th></th>"  # Empty corner
    for char in seq2:
        html += f"<th style='border: 1px solid black; padding: 5px;'>{char}</th>"
    html += "</tr>"

    for i in range(len(matrix)):
        html += "<tr>"
        if i == 0:
            html += "<th></th>"  # Empty top-left corner
        else:
            html += f"<th style='border: 1px solid black; padding: 5px;'>{seq1[i - 1]}</th>"

        for j in range(len(matrix[i])):
            color = "background-color: white;"
            if (i, j) in path:
                color = "background-color: lightcoral;"  # Highlight path
            html += f"<td style='border: 1px solid black; padding: 5px; {color}'>{matrix[i][j]}</td>"
        html += "</tr>"

    html += "</table>"
    st.markdown(html, unsafe_allow_html=True)

    st.subheader("Aligned Sequences")
    st.text(f"Sequence 1: {aligned_seq1}")
    st.text(f"Sequence 2: {aligned_seq2}")

    st.write(f"**Alignment Score:** {matrix[len(seq1)][len(seq2)]}")
