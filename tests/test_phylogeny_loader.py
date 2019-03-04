import ALifeStdDev as alsd

def test_load_phylogeny_to_networkx():
    alsd.load_phylogeny_to_networkx("example_data/asexual_phylogeny_test.csv")

if __name__ == "__main__":
    test_load_phylogeny_to_networkx()