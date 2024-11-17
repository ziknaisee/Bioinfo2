import requests
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import streamlit as st


# Function to retrieve PPI data from BioGRID
def retrieve_ppi_biogrid(target_protein):
    biogrid_url = "https://webservice.thebiogrid.org/interactions"
    params = {
        "accessKey": "af87dc666d66f0268669603d727a058f",
        "format": "json",
        "searchNames": True,
        "geneList": target_protein,
        "organism": 9606,
        "searchbiogridids": True,
        "includeInteractors": True
    }
    response = requests.get(biogrid_url, params=params)
    network = response.json()
    network_df = pd.DataFrame.from_dict(network, orient='index')
    return network_df


# Function to retrieve PPI data from STRING
def retrieve_ppi_string(target_protein):
    string_url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": target_protein,
        "species": 9606
    }
    response = requests.get(string_url, params=params)
    network = response.json()
    network_df = pd.json_normalize(network)
    return network_df


# Function to generate the network graph
def generate_network(dataframe, protein_column_a, protein_column_b):
    dataframe[protein_column_a] = [gene.upper() for gene in dataframe[protein_column_a]]
    dataframe[protein_column_b] = [gene.upper() for gene in dataframe[protein_column_b]]
    network_graph = nx.from_pandas_edgelist(dataframe, protein_column_a, protein_column_b)
    return network_graph


# Function to compute centralities
def get_centralities(network_graph):
    degree_centrality = nx.degree_centrality(network_graph)
    betweenness_centrality = nx.betweenness_centrality(network_graph)
    closeness_centrality = nx.closeness_centrality(network_graph)

    # Handle disconnected components for eigenvector centrality by using the largest connected component
    largest_component = max(nx.connected_components(network_graph), key=len)
    largest_subgraph = network_graph.subgraph(largest_component)

    try:
        # Increase max_iter for eigenvector centrality to improve convergence
        eigenvector_centrality = nx.eigenvector_centrality(largest_subgraph, max_iter=500)
    except nx.PowerIterationFailedConvergence:
        st.warning("Eigenvector centrality did not converge. Using default values.")
        eigenvector_centrality = {node: 0 for node in largest_subgraph.nodes}

    pagerank_centrality = nx.pagerank(network_graph)

    # Prepare centrality DataFrame
    centralities_df = pd.DataFrame({
        "Node": list(degree_centrality.keys()),
        "Degree Centrality": list(degree_centrality.values()),
        "Betweenness Centrality": list(betweenness_centrality.values()),
        "Closeness Centrality": list(closeness_centrality.values()),
        "Eigenvector Centrality": [eigenvector_centrality.get(node, 0) for node in degree_centrality],
        "PageRank Centrality": [pagerank_centrality.get(node, 0) for node in degree_centrality]
    }).sort_values(by="Degree Centrality", ascending=False).head(5)

    return centralities_df



st.title("Protein-Protein Interaction (PPI) Network and Centrality Measures")
st.header("Zikry Danial A22EC0298")
protein_id = st.text_input("Enter the Protein ID (e.g., TP53):")
database_choice = st.selectbox("Select Database", ["BioGRID", "STRING"])


if st.button("Retrieve PPI Data"):
    if protein_id:
        if database_choice == "BioGRID":
            ppi_data = retrieve_ppi_biogrid(protein_id)
            protein_column_a = "OFFICIAL_SYMBOL_A"
            protein_column_b = "OFFICIAL_SYMBOL_B"
            
            # Display BioGRID PPI Data and Network Info in the first column
            with st.expander("PPI Data Information (BioGRID)"):
                st.write("PPI Data (BioGRID):")
                st.dataframe(ppi_data.head())  # Display PPI Data

                # Create network graph and display info
                network_graph = generate_network(ppi_data, protein_column_a, protein_column_b)
                st.write("Number of edges:", network_graph.number_of_edges())
                st.write("Number of nodes:", network_graph.number_of_nodes())
                
                # BioGRID Network visualization (no labels on nodes, small size)
                slayout = nx.spring_layout(network_graph, seed=123)
                plt.figure(figsize=(10, 10))
                
                # Draw BioGRID network with no labels, small node size
                nx.draw(network_graph, slayout, with_labels=False, node_size=100, node_color='lightblue')

                # Highlight the nodes with high centrality
                high_centrality_nodes = [node for node, centrality in sorted(nx.degree_centrality(network_graph).items(), key=lambda x: -x[1])[:5]]
                nx.draw_networkx_nodes(network_graph, slayout, nodelist=high_centrality_nodes, node_size=1000, node_color='orange')

                # Display the visualization in Streamlit
                plt.title(f"BioGRID PPI Network for {protein_id}")
                st.pyplot(plt)

        elif database_choice == "STRING":
            ppi_data = retrieve_ppi_string(protein_id)
            protein_column_a = "preferredName_A"
            protein_column_b = "preferredName_B"
            
            # Display STRING PPI Data and Network Info in the first column
            with st.expander("PPI Data Information (STRING)"):
                st.write("PPI Data (STRING):")
                st.dataframe(ppi_data.head())  # Display PPI Data

                # Create network graph and display info
                network_graph = generate_network(ppi_data, protein_column_a, protein_column_b)
                st.write("Number of edges:", network_graph.number_of_edges())
                st.write("Number of nodes:", network_graph.number_of_nodes())
                
                # STRING Network visualization with node labels, larger node size
                slayout = nx.spring_layout(network_graph, seed=123)
                plt.figure(figsize=(10, 10))
                
                # Draw STRING network with labels, larger node size
                nx.draw(network_graph, slayout, with_labels=True, node_size=1000, node_color='lightblue', font_size=8)

                # Highlight the nodes with high centrality
                high_centrality_nodes = [node for node, centrality in sorted(nx.degree_centrality(network_graph).items(), key=lambda x: -x[1])[:5]]
                nx.draw_networkx_nodes(network_graph, slayout, nodelist=high_centrality_nodes, node_size=1000, node_color='orange')

                # Display the visualization in Streamlit
                plt.title(f"STRING PPI Network for {protein_id}")
                st.pyplot(plt)

        # Compute centralities and display in the second column
        centralities_df = get_centralities(network_graph)
        with st.expander("Centrality Measures"):
            st.write("Top 5 Nodes by Degree Centrality and Other Measures")
            st.dataframe(centralities_df)  # Display centralities in table format
    else:
        st.error("Please enter a valid protein ID.")
