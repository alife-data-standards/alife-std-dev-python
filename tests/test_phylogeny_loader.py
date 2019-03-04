import ALifeStdDev as alsd
import pytest


def test_load_phylogeny_to_networkx():
    phylo = alsd.load_phylogeny_to_networkx(
        "example_data/asexual_phylogeny_test.csv")

    assert phylo.has_node(36210211)
    assert phylo.has_node(36205850)
    assert phylo.has_edge(36205850, 36210211)
    assert phylo.nodes[1]["origin"] == "none"


def test_load_phylogeny_to_networkx_json_and_csv():
    phylo_json = alsd.load_phylogeny_to_networkx(
        "example_data/example-standard-asexual-phylogeny.json")
    phylo_csv = alsd.load_phylogeny_to_networkx(
        "example_data/example-standard-asexual-phylogeny.csv")

    assert phylo_csv.nodes == phylo_json.nodes
    assert phylo_csv.edges == phylo_json.edges


def test_failure():
    with pytest.raises(Exception):
        alsd.load_phylogeny_to_networkx("example_data/should_fail.csv")


if __name__ == "__main__":
    test_load_phylogeny_to_networkx()