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

def test_has_single_root():
    single_root_fname = "example_data/example-standard-toy-asexual-phylogeny.csv"
    multi_root_fname = "example_data/example-standard-toy-asexual-phylogeny-multi-roots.csv"
    sroot = alsd.load_phylogeny_to_networkx(single_root_fname)
    mroot = alsd.load_phylogeny_to_networkx(multi_root_fname)

    # (1) assert true when there is a single root
    assert phylo.has_single_root(sroot)
    # (2) assert false when there are multiple roots
    assert not phylo.has_single_root(mroot)

def test_get_root_ids():
    single_root_fname = "example_data/example-standard-toy-asexual-phylogeny.csv"
    multi_root_fname = "example_data/example-standard-toy-asexual-phylogeny-multi-roots.csv"
    sex_fname =  "example_data/example-standard-toy-sexual-phylogeny.csv"
    sroot = alsd.load_phylogeny_to_networkx(single_root_fname)
    mroot = alsd.load_phylogeny_to_networkx(multi_root_fname)
    sexphylo = alsd.load_phylogeny_to_networkx(sex_fname)

    sroot_ids = phylo.get_root_ids(sroot)
    mroot_ids = phylo.get_root_ids(mroot)
    sexroot_ids = phylo.get_root_ids(sexphylo)

    assert len(sroot_ids) == 1
    assert 0 in sroot_ids

    assert len(mroot_ids) == 2
    assert 0 in mroot_ids and 6 in mroot_ids

    assert len(sexroot_ids)
    assert 0 in sexroot_ids and 100 in sexroot_ids


def test_get_roots():
    single_root_fname = "example_data/example-standard-toy-asexual-phylogeny.csv"
    multi_root_fname = "example_data/example-standard-toy-asexual-phylogeny-multi-roots.csv"
    sex_fname =  "example_data/example-standard-toy-sexual-phylogeny.csv"
    sroot = alsd.load_phylogeny_to_networkx(single_root_fname)
    mroot = alsd.load_phylogeny_to_networkx(multi_root_fname)
    sexphylo = alsd.load_phylogeny_to_networkx(sex_fname)

    sroots = phylo.get_roots(sroot)
    mroots = phylo.get_roots(mroot)
    sexroots = phylo.get_roots(sexphylo)

    assert len(sroots) == 1
    assert 0 in sroots
    assert sroots[0]["id"] == 0

    assert len(mroots) == 2
    assert 0 in mroots and 6 in mroots
    assert mroots[6]["id"] == 6

    assert len(sexroots)
    assert 0 in sexroots and 100 in sexroots
    assert sexroots[0]["id"] == 0

def test_get_num_roots():
    single_root_fname = "example_data/example-standard-toy-asexual-phylogeny.csv"
    multi_root_fname = "example_data/example-standard-toy-asexual-phylogeny-multi-roots.csv"
    sex_fname =  "example_data/example-standard-toy-sexual-phylogeny.csv"
    sroot = alsd.load_phylogeny_to_networkx(single_root_fname)
    mroot = alsd.load_phylogeny_to_networkx(multi_root_fname)
    sexphylo = alsd.load_phylogeny_to_networkx(sex_fname)

    assert phylo.get_num_roots(sroot) == 1
    assert phylo.get_num_roots(mroot) == 2
    assert phylo.get_num_roots(sexphylo) == 2

def test_get_num_independent_phylogenies():
    single_root_fname = "example_data/example-standard-toy-asexual-phylogeny.csv"
    multi_root_fname = "example_data/example-standard-toy-asexual-phylogeny-multi-roots.csv"
    sex_fname =  "example_data/example-standard-toy-sexual-phylogeny.csv"
    sroot = alsd.load_phylogeny_to_networkx(single_root_fname)
    mroot = alsd.load_phylogeny_to_networkx(multi_root_fname)
    sexphylo = alsd.load_phylogeny_to_networkx(sex_fname)

    assert phylo.get_num_independent_phylogenies(sroot) == 1
    assert phylo.get_num_independent_phylogenies(mroot) == 2
    assert phylo.get_num_independent_phylogenies(sexphylo) == 1

def test_get_independent_phylogenies():
    single_root_fname = "example_data/example-standard-toy-asexual-phylogeny.csv"
    multi_root_fname = "example_data/example-standard-toy-asexual-phylogeny-multi-roots.csv"
    sex_fname =  "example_data/example-standard-toy-sexual-phylogeny.csv"
    sroot = alsd.load_phylogeny_to_networkx(single_root_fname)
    mroot = alsd.load_phylogeny_to_networkx(multi_root_fname)
    sexphylo = alsd.load_phylogeny_to_networkx(sex_fname)

    sroot_indies = phylo.get_independent_phylogenies(sroot)
    mroot_indies = phylo.get_independent_phylogenies(mroot)
    sexphylo_indies = phylo.get_independent_phylogenies(sexphylo)

    assert len(sroot_indies) == 1
    assert sroot_indies[0].nodes == sroot.nodes

    assert len(mroot_indies) == 2
    assert (set(mroot_indies[0].nodes) == set([6,7,8])) or (set(mroot_indies[1].nodes) == set([6,7,8]))
    assert (set(mroot_indies[0].nodes) == set([0,1,2,3,4,5])) or (set(mroot_indies[1].nodes) == set([0,1,2,3,4,5]))

    assert len(sexphylo_indies) == 1
    assert sexphylo_indies[0].nodes == sexphylo.nodes

def test_get_leaf_taxa():
    single_root_fname = "example_data/example-standard-toy-asexual-phylogeny.csv"
    multi_root_fname = "example_data/example-standard-toy-asexual-phylogeny-multi-roots.csv"
    sex_fname =  "example_data/example-standard-toy-sexual-phylogeny.csv"
    sroot = alsd.load_phylogeny_to_networkx(single_root_fname)
    mroot = alsd.load_phylogeny_to_networkx(multi_root_fname)
    sexphylo = alsd.load_phylogeny_to_networkx(sex_fname)

    sroot_leafs = phylo.get_leaf_taxa(sroot)
    mroot_leafs = phylo.get_leaf_taxa(mroot)
    sexphylo_leafs = phylo.get_leaf_taxa(sexphylo)

    assert len(sroot_leafs) == 3
    assert set(sroot_leafs.keys()) == set([3,4,5])

    assert len(mroot_leafs) == 4
    assert set(mroot_leafs.keys()) == set([3,4,5,8])

    assert len(sexphylo_leafs) == 1
    assert set(sexphylo_leafs.keys()) == set([5])

def test_get_leaf_taxa_ids():
    single_root_fname = "example_data/example-standard-toy-asexual-phylogeny.csv"
    multi_root_fname = "example_data/example-standard-toy-asexual-phylogeny-multi-roots.csv"
    sex_fname =  "example_data/example-standard-toy-sexual-phylogeny.csv"
    sroot = alsd.load_phylogeny_to_networkx(single_root_fname)
    mroot = alsd.load_phylogeny_to_networkx(multi_root_fname)
    sexphylo = alsd.load_phylogeny_to_networkx(sex_fname)

    sroot_leafs = phylo.get_leaf_taxa_ids(sroot)
    mroot_leafs = phylo.get_leaf_taxa_ids(mroot)
    sexphylo_leafs = phylo.get_leaf_taxa_ids(sexphylo)

    assert len(sroot_leafs) == 3
    assert set(sroot_leafs) == set([3,4,5])

    assert len(mroot_leafs) == 4
    assert set(mroot_leafs) == set([3,4,5,8])

    assert len(sexphylo_leafs) == 1
    assert set(sexphylo_leafs) == set([5])

def test_get_extant_taxa_ids():
    pruned_fname = "example_data/example-standard-toy-asexual-phylogeny.csv"
    unpruned_fname = "example_data/example-standard-toy-asexual-phylogeny-not-pruned.csv"
    pruned_phylo = alsd.load_phylogeny_to_networkx(pruned_fname)
    unpruned_phylo = alsd.load_phylogeny_to_networkx(unpruned_fname)

    pruned_extant = phylo.get_extant_taxa_ids(pruned_phylo)
    unpruned_extant = phylo.get_extant_taxa_ids(unpruned_phylo)

    assert set(pruned_extant) == set([3,4,5])
    assert set(unpruned_extant) == set([3,4,5])

def test_get_extant_taxa():
    pruned_fname = "example_data/example-standard-toy-asexual-phylogeny.csv"
    unpruned_fname = "example_data/example-standard-toy-asexual-phylogeny-not-pruned.csv"
    pruned_phylo = alsd.load_phylogeny_to_networkx(pruned_fname)
    unpruned_phylo = alsd.load_phylogeny_to_networkx(unpruned_fname)

    pruned_extant = phylo.get_extant_taxa(pruned_phylo)
    unpruned_extant = phylo.get_extant_taxa(unpruned_phylo)

    assert set(pruned_extant.keys()) == set([3,4,5])
    assert set(unpruned_extant.keys()) == set([3,4,5])

def test_extract_asexual_lineage():
    multi_root_fname = "example_data/example-standard-toy-asexual-phylogeny-multi-roots.csv"
    sex_fname =  "example_data/example-standard-toy-sexual-phylogeny.csv"
    mroot = alsd.load_phylogeny_to_networkx(multi_root_fname)
    sexphylo = alsd.load_phylogeny_to_networkx(sex_fname)

    # check lineage of 5 8
    lineage_5 = phylo.extract_asexual_lineage(mroot, 5)
    lineage_8 = phylo.extract_asexual_lineage(mroot, 8)

    assert len(lineage_5.nodes) == 3
    assert set(lineage_5.nodes) == set([0,2,5])
    assert len(lineage_8.nodes) == 3
    assert set(lineage_8.nodes) == set([6,7,8])

    # make sure sexphylo fails
    with pytest.raises(Exception):
        phylo.extract_asexual_lineage(sexphylo, 5)

    # make sure that invalid taxa id fails
    with pytest.raises(Exception):
        phylo.extract_asexual_lineage(mroot, -99999)

if __name__ == "__main__":
    test_all_taxa_have_attribute()
    test_is_asexual()
    test_get_root_ids()
    test_get_roots()
    test_get_num_roots()
    test_get_num_independent_phylogenies()
    test_get_independent_phylogenies()
    test_get_leaf_taxa()
    test_get_leaf_taxa_ids()
    test_get_extant_taxa_ids()
    test_get_extant_taxa()
    test_extract_asexual_lineage()