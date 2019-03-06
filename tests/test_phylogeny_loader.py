import ALifeStdDev.phylogeny as phylodev
import pytest


def test_load_phylogeny_to_networkx():
    phylo = phylodev.load_phylogeny_to_networkx(
        "example_data/asexual_phylogeny_test.csv")

    assert phylo.has_node(36210211)
    assert phylo.has_node(36205850)
    assert phylo.has_edge(36205850, 36210211)
    assert phylo.nodes[1]["origin"] == "none"
    assert phylo.nodes[1]["origin_time"] == -1
    assert phylo.nodes[1]["src"] == "div:ext"


def test_load_phylogeny_to_networkx_json_and_csv():
    phylo_json = phylodev.load_phylogeny_to_networkx(
        "example_data/example-standard-asexual-phylogeny.json")
    phylo_csv = phylodev.load_phylogeny_to_networkx(
        "example_data/example-standard-asexual-phylogeny.csv")

    assert set(phylo_csv.nodes) == set(phylo_json.nodes)
    assert set(phylo_csv.edges) == set(phylo_json.edges)


def test_failure():
    with pytest.raises(Exception):
        phylodev.load_phylogeny_to_networkx("example_data/should_fail.csv")


if __name__ == "__main__":
    test_load_phylogeny_to_networkx()