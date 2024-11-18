import streamlit as st
import numpy as np

# Function to calculate Needleman-Wunsch (Global Alignment)
def needleman_wunsch(seq1, seq2, match_score, mismatch_penalty, gap_penalty):
    n, m = len(seq1), len(seq2)
    matrix = np.zeros((n + 1, m + 1), dtype=int)
    backtrace = np.zeros((n + 1, m + 1), dtype=tuple)

    # Initialize scoring matrix
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

            if matrix[i][j] == match:
                backtrace[i][j] = (i - 1, j - 1)
            elif matrix[i][j] == delete:
                backtrace[i][j] = (i - 1, j)
            else:
                backtrace[i][j] = (i, j - 1)

    # Traceback to get alignment
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

# Function to calculate Smith-Waterman (Local Alignment)
def smith_waterman(seq1, seq2, match_score, mismatch_penalty, gap_penalty):
    n, m = len(seq1), len(seq2)
    matrix = np.zeros((n + 1, m + 1), dtype=int)
    backtrace = np.zeros((n + 1, m + 1), dtype=tuple)

    # Fill the scoring matrix
    max_score = 0
    max_pos = None
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty)
            delete = matrix[i - 1][j] + gap_penalty
            insert = matrix[i][j - 1] + gap_penalty
            matrix[i][j] = max(0, match, delete, insert)

            if matrix[i][j] > max_score:
                max_score = matrix[i][j]
                max_pos = (i, j)

    # Traceback to get alignment
    aligned_seq1, aligned_seq2 = "", ""
    i, j = max_pos
    path = []
    while matrix[i][j] > 0:
        path.append((i, j))
        if matrix[i][j] == matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty):
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            i, j = i - 1, j - 1
        elif matrix[i][j] == matrix[i - 1][j] + gap_penalty:
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

seq1 = st.sidebar.text_input("Enter Sequence 1:", "FKHMPLNE")
seq2 = st.sidebar.text_input("Enter Sequence 2:", "FMDTPLNE")
match_score = st.sidebar.number_input("Match score:", value=1, step=1)
mismatch_penalty = st.sidebar.number_input("Mismatch penalty:", value=-1, step=1)
gap_penalty = st.sidebar.number_input("Gap penalty:", value=-2, step=1)

alignment_algorithm = st.selectbox("Choose Alignment Algorithm:", ["Needleman-Wunsch (Global)", "Smith-Waterman (Local)"])

if st.button("Run Alignment"):
    if alignment_algorithm == "Needleman-Wunsch (Global)":
        matrix, aligned_seq1, aligned_seq2, path = needleman_wunsch(seq1, seq2, match_score, mismatch_penalty, gap_penalty)
        alignment_type = "Global Alignment (Needleman-Wunsch)"
    else:
        matrix, aligned_seq1, aligned_seq2, path = smith_waterman(seq1, seq2, match_score, mismatch_penalty, gap_penalty)
        alignment_type = "Local Alignment (Smith-Waterman)"

    st.subheader(f"{alignment_type} Results")
    st.write("**Aligned Sequences:**")
    st.write("Sequence 1: " + aligned_seq1)
    st.write("Sequence 2: " + aligned_seq2)

    st.write("**Scoring Matrix:**")
    matrix_html = "<table style='border-collapse: collapse;'>"
    matrix_html += "<tr><th></th><th></th>" + "".join([f"<th>{c}</th>" for c in seq2]) + "</tr>"
    for i in range(len(matrix)):
        row_html = f"<tr><th>{seq1[i - 1] if i > 0 else ''}</th>" if i > 0 else "<tr><th></th>"
        for j in range(len(matrix[i])):
            color = "background-color: red;" if (i, j) in path else "background-color: white;"
            row_html += f"<td style='border: 1px solid black; {color}'>{matrix[i][j]}</td>"
        row_html += "</tr>"
        matrix_html += row_html
    matrix_html += "</table>"
    st.markdown(matrix_html, unsafe_allow_html=True)
