import ALifeStdDev.phylogeny as phylodev
import pytest
import networkx as nx


def test_all_taxa_have_attribute():
    fname = "example_data/example-standard-toy-asexual-phylogeny.csv"
    phylogeny = phylodev.load_phylogeny_to_networkx(fname)
    # assert that all nodes have the 'trait1' atttribute
    assert phylodev.all_taxa_have_attribute(phylogeny, "trait_a")
    assert phylodev.all_taxa_have_attribute(phylogeny, "trait_b")
    assert phylodev.all_taxa_have_attribute(phylogeny, "trait_c")
    # assert that all nodes do not have some garbage attributes (that they
    # shouldn't have)
    assert not phylodev.all_taxa_have_attribute(phylogeny, "nothing_should_have_this_garbage_test_trait")


def test_is_asexual():
    asex_fname = "example_data/example-standard-toy-asexual-phylogeny.csv"
    sex_fname = "example_data/example-standard-toy-sexual-phylogeny.csv"

    asex_phylogeny = phylodev.load_phylogeny_to_networkx(asex_fname)
    sex_phylogeny = phylodev.load_phylogeny_to_networkx(sex_fname)

    # assert that asexual phylogeny is asexual
    assert phylodev.is_asexual(asex_phylogeny)
    assert not phylodev.is_asexual(sex_phylogeny)


def test_has_single_root():
    single_root_fname = "example_data/example-standard-toy-asexual-phylogeny.csv"
    multi_root_fname = "example_data/example-standard-toy-asexual-phylogeny-multi-roots.csv"
    sroot = phylodev.load_phylogeny_to_networkx(single_root_fname)
    mroot = phylodev.load_phylogeny_to_networkx(multi_root_fname)

    # (1) assert true when there is a single root
    assert phylodev.has_single_root(sroot)
    # (2) assert false when there are multiple roots
    assert not phylodev.has_single_root(mroot)


def test_get_root_ids():
    single_root_fname = "example_data/example-standard-toy-asexual-phylogeny.csv"
    multi_root_fname = "example_data/example-standard-toy-asexual-phylogeny-multi-roots.csv"
    sex_fname = "example_data/example-standard-toy-sexual-phylogeny.csv"
    sroot = phylodev.load_phylogeny_to_networkx(single_root_fname)
    mroot = phylodev.load_phylogeny_to_networkx(multi_root_fname)
    sexphylo = phylodev.load_phylogeny_to_networkx(sex_fname)

    sroot_ids = phylodev.get_root_ids(sroot)
    mroot_ids = phylodev.get_root_ids(mroot)
    sexroot_ids = phylodev.get_root_ids(sexphylo)

    assert len(sroot_ids) == 1
    assert 0 in sroot_ids

    assert len(mroot_ids) == 2
    assert 0 in mroot_ids and 6 in mroot_ids

    assert len(sexroot_ids)
    assert 0 in sexroot_ids and 100 in sexroot_ids


def test_get_roots():
    single_root_fname = "example_data/example-standard-toy-asexual-phylogeny.csv"
    multi_root_fname = "example_data/example-standard-toy-asexual-phylogeny-multi-roots.csv"
    sex_fname = "example_data/example-standard-toy-sexual-phylogeny.csv"
    sroot = phylodev.load_phylogeny_to_networkx(single_root_fname)
    mroot = phylodev.load_phylogeny_to_networkx(multi_root_fname)
    sexphylo = phylodev.load_phylogeny_to_networkx(sex_fname)

    sroots = phylodev.get_roots(sroot)
    mroots = phylodev.get_roots(mroot)
    sexroots = phylodev.get_roots(sexphylo)

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
    sex_fname = "example_data/example-standard-toy-sexual-phylogeny.csv"
    sroot = phylodev.load_phylogeny_to_networkx(single_root_fname)
    mroot = phylodev.load_phylogeny_to_networkx(multi_root_fname)
    sexphylo = phylodev.load_phylogeny_to_networkx(sex_fname)

    assert phylodev.get_num_roots(sroot) == 1
    assert phylodev.get_num_roots(mroot) == 2
    assert phylodev.get_num_roots(sexphylo) == 2


def test_get_num_independent_phylogenies():
    single_root_fname = "example_data/example-standard-toy-asexual-phylogeny.csv"
    multi_root_fname = "example_data/example-standard-toy-asexual-phylogeny-multi-roots.csv"
    sex_fname = "example_data/example-standard-toy-sexual-phylogeny.csv"
    sroot = phylodev.load_phylogeny_to_networkx(single_root_fname)
    mroot = phylodev.load_phylogeny_to_networkx(multi_root_fname)
    sexphylo = phylodev.load_phylogeny_to_networkx(sex_fname)

    assert phylodev.get_num_independent_phylogenies(sroot) == 1
    assert phylodev.get_num_independent_phylogenies(mroot) == 2
    assert phylodev.get_num_independent_phylogenies(sexphylo) == 1


def test_get_independent_phylogenies():
    single_root_fname = "example_data/example-standard-toy-asexual-phylogeny.csv"
    multi_root_fname = "example_data/example-standard-toy-asexual-phylogeny-multi-roots.csv"
    sex_fname = "example_data/example-standard-toy-sexual-phylogeny.csv"
    sroot = phylodev.load_phylogeny_to_networkx(single_root_fname)
    mroot = phylodev.load_phylogeny_to_networkx(multi_root_fname)
    sexphylo = phylodev.load_phylogeny_to_networkx(sex_fname)

    sroot_indies = phylodev.get_independent_phylogenies(sroot)
    mroot_indies = phylodev.get_independent_phylogenies(mroot)
    sexphylo_indies = phylodev.get_independent_phylogenies(sexphylo)

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
    sex_fname = "example_data/example-standard-toy-sexual-phylogeny.csv"
    sroot = phylodev.load_phylogeny_to_networkx(single_root_fname)
    mroot = phylodev.load_phylogeny_to_networkx(multi_root_fname)
    sexphylo = phylodev.load_phylogeny_to_networkx(sex_fname)

    sroot_leafs = phylodev.get_leaf_taxa(sroot)
    mroot_leafs = phylodev.get_leaf_taxa(mroot)
    sexphylo_leafs = phylodev.get_leaf_taxa(sexphylo)

    assert len(sroot_leafs) == 3
    assert set(sroot_leafs.keys()) == set([3,4,5])

    assert len(mroot_leafs) == 4
    assert set(mroot_leafs.keys()) == set([3,4,5,8])

    assert len(sexphylo_leafs) == 1
    assert set(sexphylo_leafs.keys()) == set([5])


def test_get_leaf_taxa_ids():
    single_root_fname = "example_data/example-standard-toy-asexual-phylogeny.csv"
    multi_root_fname = "example_data/example-standard-toy-asexual-phylogeny-multi-roots.csv"
    sex_fname = "example_data/example-standard-toy-sexual-phylogeny.csv"
    sroot = phylodev.load_phylogeny_to_networkx(single_root_fname)
    mroot = phylodev.load_phylogeny_to_networkx(multi_root_fname)
    sexphylo = phylodev.load_phylogeny_to_networkx(sex_fname)

    sroot_leafs = phylodev.get_leaf_taxa_ids(sroot)
    mroot_leafs = phylodev.get_leaf_taxa_ids(mroot)
    sexphylo_leafs = phylodev.get_leaf_taxa_ids(sexphylo)

    assert len(sroot_leafs) == 3
    assert set(sroot_leafs) == set([3,4,5])

    assert len(mroot_leafs) == 4
    assert set(mroot_leafs) == set([3,4,5,8])

    assert len(sexphylo_leafs) == 1
    assert set(sexphylo_leafs) == set([5])


def test_get_extant_taxa_ids():
    pruned_fname = "example_data/example-standard-toy-asexual-phylogeny.csv"
    unpruned_fname = "example_data/example-standard-toy-asexual-phylogeny-not-pruned.csv"
    pruned_phylo = phylodev.load_phylogeny_to_networkx(pruned_fname)
    unpruned_phylo = phylodev.load_phylogeny_to_networkx(unpruned_fname)

    pruned_extant = phylodev.get_extant_taxa_ids(pruned_phylo)
    unpruned_extant = phylodev.get_extant_taxa_ids(unpruned_phylo)

    assert set(pruned_extant) == set([3, 4, 5])
    assert set(unpruned_extant) == set([3, 4, 5])


def test_get_extant_taxa():
    pruned_fname = "example_data/example-standard-toy-asexual-phylogeny.csv"
    unpruned_fname = "example_data/example-standard-toy-asexual-phylogeny-not-pruned.csv"
    pruned_phylo = phylodev.load_phylogeny_to_networkx(pruned_fname)
    unpruned_phylo = phylodev.load_phylogeny_to_networkx(unpruned_fname)

    pruned_extant = phylodev.get_extant_taxa(pruned_phylo)
    unpruned_extant = phylodev.get_extant_taxa(unpruned_phylo)

    assert set(pruned_extant.keys()) == set([3, 4, 5])
    assert set(unpruned_extant.keys()) == set([3, 4, 5])


def test_get_extant_at_time():
    pruned_fname = "example_data/example-standard-toy-asexual-phylogeny.csv"

    pruned_phylo = phylodev.load_phylogeny_to_networkx(pruned_fname)

    # Should raise exception because no origin_time
    with pytest.raises(Exception):
        phylodev.get_extant_taxa_ids(pruned_phylo, time=1)
    with pytest.raises(Exception):
        phylodev.get_extant_taxa(pruned_phylo, time=1)

    extant = phylodev.get_extant_taxa_ids(pruned_phylo, time=0,
                                          origin_attribute="trait_a")
    assert set(extant) == set([0, 2, 5])
    extant = phylodev.get_extant_taxa_ids(pruned_phylo, time=1,
                                          origin_attribute="trait_a")
    assert set(extant) == set([1, 2, 3, 4, 5])
    extant = phylodev.get_extant_taxa_ids(pruned_phylo, time=2,
                                          origin_attribute="trait_a")
    assert set(extant) == set([3, 4, 5])


def test_validate_destruction_time():
    g = nx.DiGraph()
    g.add_node("A")

    with pytest.raises(Exception):
        phylodev.validate_destruction_time(g)


def test_extract_asexual_lineage():
    multi_root_fname = "example_data/example-standard-toy-asexual-phylogeny-multi-roots.csv"
    sex_fname = "example_data/example-standard-toy-sexual-phylogeny.csv"
    mroot = phylodev.load_phylogeny_to_networkx(multi_root_fname)
    sexphylo = phylodev.load_phylogeny_to_networkx(sex_fname)

    # check lineage of 5 8
    lineage_5 = phylodev.extract_asexual_lineage(mroot, 5)
    lineage_8 = phylodev.extract_asexual_lineage(mroot, 8)

    assert len(lineage_5.nodes) == 3
    assert set(lineage_5.nodes) == set([0, 2, 5])
    assert len(lineage_8.nodes) == 3
    assert set(lineage_8.nodes) == set([6, 7, 8])

    # make sure sexphylo fails
    with pytest.raises(Exception):
        phylodev.extract_asexual_lineage(sexphylo, 5)

    # make sure that invalid taxa id fails
    with pytest.raises(Exception):
        phylodev.extract_asexual_lineage(mroot, -99999)


def test_is_asexual_lineage():
    multi_root_fname = "example_data/example-standard-toy-asexual-phylogeny-multi-roots.csv"
    sex_fname = "example_data/example-standard-toy-sexual-phylogeny.csv"
    mroot = phylodev.load_phylogeny_to_networkx(multi_root_fname)
    sexphylo = phylodev.load_phylogeny_to_networkx(sex_fname)

    # check lineage of 5 8
    lineage_5 = phylodev.extract_asexual_lineage(mroot, 5)
    lineage_8 = phylodev.extract_asexual_lineage(mroot, 8)

    assert phylodev.is_asexual_lineage(lineage_5)
    assert phylodev.is_asexual_lineage(lineage_8)
    assert not phylodev.is_asexual_lineage(mroot)
    assert not phylodev.is_asexual_lineage(sexphylo)


def test_abstract_asexual_lineage():
    toy_lineage_fname = "example_data/example-standard-toy-asexual-lineage.csv"
    lineage = phylodev.load_phylogeny_to_networkx(toy_lineage_fname)

    # Abstract lineage by genotype
    abstract_lineage_genotype = phylodev.abstract_asexual_lineage(lineage,
                                                                  ["genotype"])
    assert len(abstract_lineage_genotype.nodes) == 4
    # Abstract lineage by phenotype
    abstract_lineage_phenotype = phylodev.abstract_asexual_lineage(lineage,
                                                                   ["trait_a",
                                                                    "trait_b"])
    assert len(abstract_lineage_phenotype) == 3
    # Abstract lineage by trait_b
    abstract_lineage_tb = phylodev.abstract_asexual_lineage(lineage, ["trait_b"])
    assert len(abstract_lineage_tb) == 2

def test_extract_asexual_lod():
    single_lineage_fname = "example_data/example-standard-toy-asexual-lineage.csv"
    single_root_fname = "example_data/example-standard-toy-asexual-phylogeny.csv"
    multi_root_fname = "example_data/example-standard-toy-asexual-phylogeny-multi-roots.csv"
    sex_fname = "example_data/example-standard-toy-sexual-phylogeny.csv"
    
    lineage = phylodev.load_phylogeny_to_networkx(single_lineage_fname)
    sroot = phylodev.load_phylogeny_to_networkx(single_root_fname)
    mroot = phylodev.load_phylogeny_to_networkx(multi_root_fname)
    sexphylo = phylodev.load_phylogeny_to_networkx(sex_fname)

    # check lod for a branching lineage
    sroot_lod = phylodev.extract_asexual_lod(sroot)
    assert len(sroot_lod.nodes) == 3
    assert set(sroot_lod.nodes) == set([5, 2, 0])

    # check lod for a multiroot lineage 
    mroot_lod = phylodev.extract_asexual_lod(mroot)
    assert len(mroot_lod.nodes) == 3
    assert set(mroot_lod.nodes) == set([8, 7, 6])

    # make sure something without living taxa fails
    with pytest.raises(Exception):
        phylodev.extract_asexual_lod(lineage)

    # make sure sexphylo fails
    with pytest.raises(Exception):
        phylodev.extract_asexual_lod(sexphylo)

def test_get_mrca_id_asexual():
    single_root_fname = "example_data/example-standard-toy-asexual-phylogeny.csv"
    multi_root_fname = "example_data/example-standard-toy-asexual-phylogeny-multi-roots.csv"

    sroot = phylodev.load_phylogeny_to_networkx(single_root_fname)
    mroot = phylodev.load_phylogeny_to_networkx(multi_root_fname)

    mrca_id = phylodev.get_mrca_id_asexual(sroot)
    assert mrca_id == 0

    mrca_id = phylodev.get_mrca_id_asexual(sroot, [3,4,5])
    assert mrca_id == 0

    mrca_id = phylodev.get_mrca_id_asexual(sroot, [3,4])
    assert mrca_id == 1

    mrca_id = phylodev.get_mrca_id_asexual(sroot, [5,0])
    assert mrca_id == 0

    mrca_id = phylodev.get_mrca_id_asexual(sroot, [0,1,2,3,4,5])
    assert mrca_id == 0

    mrca_id = phylodev.get_mrca_id_asexual(sroot, [2])
    assert mrca_id == 2

    mrca_id = phylodev.get_mrca_id_asexual(mroot)
    assert mrca_id == -1

    mrca_id = phylodev.get_mrca_id_asexual(mroot, [8, 8])
    assert mrca_id == 8

    mrca_id = phylodev.get_mrca_id_asexual(mroot, [7, 8])
    assert mrca_id == 7

    mrca_id = phylodev.get_mrca_id_asexual(mroot, [6,8])
    assert mrca_id == 6


def test_has_common_ancestor_asexual():
    single_root_fname = "example_data/example-standard-toy-asexual-phylogeny.csv"
    multi_root_fname = "example_data/example-standard-toy-asexual-phylogeny-multi-roots.csv"
    sroot = phylodev.load_phylogeny_to_networkx(single_root_fname)
    mroot = phylodev.load_phylogeny_to_networkx(multi_root_fname)
    assert phylodev.has_common_ancestor_asexual(sroot)
    assert not phylodev.has_common_ancestor_asexual(mroot)


def test_get_pairwise_distances():
    single_root_fname = "example_data/example-standard-toy-asexual-phylogeny.csv"
    multi_root_fname = "example_data/example-standard-toy-asexual-phylogeny-multi-roots.csv"
    sroot = phylodev.load_phylogeny_to_networkx(single_root_fname)
    mroot = phylodev.load_phylogeny_to_networkx(multi_root_fname)
    assert set([2, 4, 4]) == set(phylodev.get_pairwise_distances(sroot, [3, 4, 5]))
    with pytest.raises(nx.NetworkXNoPath):
        phylodev.get_pairwise_distances(mroot, [3, 4, 8])


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
    test_get_extant_at_time()
    test_extract_asexual_lineage()
    test_is_asexual_lineage()
    test_abstract_asexual_lineage()
    test_extract_asexual_lod()
    test_get_mrca_id_asexual()
    test_has_common_ancestor_asexual()
    test_get_pairwise_distances()