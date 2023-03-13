import ALifeStdDev.phylogeny as phylodev
import pytest
import pandas as pd


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


def test_networkx_to_pandas():
    phylo = phylodev.load_phylogeny_to_networkx(
        "example_data/asexual_phylogeny_test.csv")
    df = phylodev.networkx_to_pandas_df(phylo, {"origin_time": "origin_time",
                                                "gestation": "gest_time"})
    df.set_index("id", inplace=True)
    assert df.loc[36205979, "origin_time"] == 199977
    assert df.loc[36205979, "gestation"] == 0
    assert df.loc[36205979, "ancestor_list"] == [36204695]

    df = phylodev.networkx_to_pandas_df(phylo)
    df.set_index("id", inplace=True)
    assert df.loc[36205979, "ancestor_list"] == [36204695]


def test_pandas_to_networkx():
    df = pd.read_csv("example_data/example-standard-toy-asexual-phylogeny.csv")
    g = phylodev.pandas_df_to_networkx(df)
    assert 1 in g[0]
    assert 2 in g[0]
    assert 3 in g[1]
    assert 4 in g[1]
    assert 5 in g[2]
    assert 2 not in g[5]
    

if __name__ == "__main__":
    # test_load_phylogeny_to_networkx()
    # test_networkx_to_pandas()
    test_pandas_to_networkx()