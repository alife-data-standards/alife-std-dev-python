import networkx as nx
import pandas as pd

def load_phylogeny_to_networkx(filename):
    phylogeny = nx.Graph()

    if filename.endswith(".csv"):
        data = pd.read_csv(filename)
    elif filename.endswith(".json"):    
        data = pd.read_json(filename)

    data.set_index("id")