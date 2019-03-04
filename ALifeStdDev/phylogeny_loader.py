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
        data.loc[:,"ancestor_list"] = data.loc[:,"ancestor_list"].apply(ast.literal_eval)
        data.set_index("id", inplace=True)
    elif filename.endswith(".json"):
        data = pd.read_json(filename, orient="index", convert_dates=False)

    # print(data)

    # Have to do this in two passes, (one to add nodes, one to
    # add eges) because order in file is not required/guaranteed
    for taxon_id in data.index:
        phylogeny.add_node(taxon_id)

    for taxon_id in data.index:

        for ancestor in data.loc[taxon_id, "ancestor_list"]:
            try:
                ancestor = int(ancestor)
                if phylogeny.has_node(ancestor):
                    phylogeny.add_edge(ancestor, taxon_id)        
                else:
                    raise Exception(f"{taxon_id}'s ancestor, {ancestor}, is not in this file.")

            except ValueError:
                phylogeny.nodes[taxon_id]["origin"] = ancestor

    return phylogeny