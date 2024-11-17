import requests
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import streamlit as st


match_score=1
miss_score=-1
gap_penalty=-2

def initialMatrix (Sequence1,Sequence2):
    
    row=len(Sequence1)+1
    column=len(Sequence2)+1
    matrix=[[0]*(column+1) for_in range(row+1)] #intialize for matrix dimensions (m+1)*(n+1)
    
    return initialMatrix