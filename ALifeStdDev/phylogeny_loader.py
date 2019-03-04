import networkx as nx
import pandas as pd
import ast

def load_phylogeny_to_networkx(filename):
    """
    Loads a phylogeny in standards format (in csv or json) from the file 
    specified by the filename parameter. Returns the phylogeny as a 
    networkx digraph (with edges going from parent to child).

    Special ancestor values are encoded in an "origin" field within nodes.
    Example: `my_phylogeny.nodes["A"]["origin"]` would return the origin of 
    node "A"
    """
    phylogeny = nx.DiGraph()

    if filename.endswith(".csv"):
        data = pd.read_csv(filename)
    elif filename.endswith(".json"):    
        data = pd.read_json(filename)

    data.set_index("id", inplace=True)

    # Have to do this in two passes, (one to add nodes, one to
    # add eges) because order in file is not required/guaranteed
    for id in data.index:
        phylogeny.add_node(id)

    for id in data.index:
        ancestors = ast.literal_eval(data.loc[id, "ancestor_list"])

        for ancestor in ancestors:
            try:
                ancestor = int(ancestor)
                if phylogeny.has_node(ancestor):
                    phylogeny.add_edge(ancestor, id)        
                else:
                    raise Exception(f"{id}'s ancestor, {ancestor}, is not in this file.")

            except ValueError as e:
                phylogeny.nodes[id]["origin"] = ancestors[0]

    return phylogeny