import ALifeStdDev as alsd
from ALifeStdDev import phylogeny_utils as phylo
import pytest

def test_all_taxa_have_attribute():
    fname = "example_data/example-standard-toy-asexual-phylogeny.csv"
    phylogeny = alsd.load_phylogeny_to_networkx(fname)
    # assert that all nodes have the 'trait1' atttribute
    assert phylo.all_taxa_have_attribute(phylogeny, "trait_a")
    assert phylo.all_taxa_have_attribute(phylogeny, "trait_b")
    assert phylo.all_taxa_have_attribute(phylogeny, "trait_c")
    # assert that all nodes do not have some garbage attributes (that they shouldn't
    # have)
    assert not phylo.all_taxa_have_attribute(phylogeny, "nothing_should_have_this_garbage_test_trait")

def test_is_asexual():
    asex_fname = "example_data/example-standard-toy-asexual-phylogeny.csv"
    sex_fname =  "example_data/example-standard-toy-sexual-phylogeny.csv"

    asex_phylogeny = alsd.load_phylogeny_to_networkx(asex_fname)
    sex_phylogeny = alsd.load_phylogeny_to_networkx(sex_fname)

    # assert that asexual phylogeny is asexual
    assert phylo.is_asexual(asex_phylogeny)
    assert not phylo.is_asexual(sex_phylogeny)


if __name__ == "__main__":
    test_all_taxa_have_attribute()