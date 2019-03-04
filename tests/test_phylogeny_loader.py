import ALifeStdDev as alsd

def test_load_phylogeny_to_networkx():
    phylo = alsd.load_phylogeny_to_networkx("example_data/asexual_phylogeny_test.csv")

def test_load_phylogeny_to_networkx_json_and_csv():
    phylo_json = alsd.load_phylogeny_to_networkx("example_data/example-standard-asexual-phylogeny.json")
    phylo_csv = alsd.load_phylogeny_to_networkx("example_data/example-standard-asexual-phylogeny.csv")

    assert phylo_csv.nodes == phylo_json.nodes
    assert phylo_csv.edges == phylo_json.edges

if __name__ == "__main__":
    test_load_phylogeny_to_networkx()