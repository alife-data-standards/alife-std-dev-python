import ALifeStdDev as alsd

def test_load_phylogeny_to_networkx_csv():
    alsd.load_phylogeny_to_networkx("example_data/asexual_phylogeny_test.csv")

def test_load_phylogeny_to_networkx_json():
    alsd.load_phylogeny_to_networkx("example_data/example-standard-asexual-phylogeny.json")


if __name__ == "__main__":
    test_load_phylogeny_to_networkx_json()