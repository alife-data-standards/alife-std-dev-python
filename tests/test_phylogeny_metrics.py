import ALifeStdDev.phylogeny as phylodev
import pytest

def test_get_asexual_lineage_length():
    toy_lineage_fname = "example_data/example-standard-toy-asexual-lineage.csv"
    lineage = phylodev.load_phylogeny_to_networkx(toy_lineage_fname)

    length = phylodev.get_asexual_lineage_length(lineage)
    assert length == 8

def test_get_asexual_lineage_num_discrete_state_changes():
    toy_lineage_fname = "example_data/example-standard-toy-asexual-lineage.csv"
    lineage = phylodev.load_phylogeny_to_networkx(toy_lineage_fname)

    assert 4 == phylodev.get_asexual_lineage_num_discrete_state_changes(lineage, attribute_list=["genotype"])
    assert 4 == phylodev.get_asexual_lineage_num_discrete_state_changes(lineage, attribute_list=["genotype","trait_a"])
    assert 3 == phylodev.get_asexual_lineage_num_discrete_state_changes(lineage, attribute_list=["trait_a"])
    assert 2 == phylodev.get_asexual_lineage_num_discrete_state_changes(lineage, attribute_list=["trait_b"])

    with pytest.raises(Exception):
        phylodev.get_asexual_lineage_num_discrete_state_changes(lineage, attribute_list=["garbage_attribute_that_nothing_should_have"])


def test_get_asexual_lineage_num_discrete_unique_states():
    toy_lineage_fname = "example_data/example-standard-toy-asexual-lineage.csv"
    lineage = phylodev.load_phylogeny_to_networkx(toy_lineage_fname)
    assert 4 == phylodev.get_asexual_lineage_num_discrete_unique_states(lineage, attribute_list=["genotype"])
    assert 2 == phylodev.get_asexual_lineage_num_discrete_unique_states(lineage, attribute_list=["trait_b"])

    with pytest.raises(Exception):
        phylodev.get_asexual_lineage_num_discrete_unique_states(lineage, attribute_list=["garbage_attribute_that_nothing_should_have"])

def test_get_asexual_lineage_mutation_accumulation():
    toy_lineage_fname = "example_data/example-standard-toy-asexual-lineage.csv"
    lineage = phylodev.load_phylogeny_to_networkx(toy_lineage_fname)

    mut_dist = phylodev.get_asexual_lineage_mutation_accumulation(lineage, mutation_attributes=["sub_mut_cnt","reverse_mut_cnt"])
    assert mut_dist["sub_mut_cnt"] == 2
    assert mut_dist["reverse_mut_cnt"] == 1

    with pytest.raises(Exception):
        mut_dist = phylodev.get_asexual_lineage_mutation_accumulation(lineage, mutation_attributes=["garbage_attribute_that_nothing_should_have"])

def test_get_mrca_tree_depth():
    single_root_fname = "example_data/example-standard-toy-asexual-phylogeny.csv"
    sroot = phylodev.load_phylogeny_to_networkx(single_root_fname)

    depth = phylodev.get_mrca_tree_depth_asexual(sroot)
    assert depth == 0

    depth = phylodev.get_mrca_tree_depth_asexual(sroot, [3,4,5])
    assert depth == 0

    depth = phylodev.get_mrca_tree_depth_asexual(sroot, [3,4])
    assert depth == 1

    depth = phylodev.get_mrca_tree_depth_asexual(sroot, [5,0])
    assert depth == 0

    depth = phylodev.get_mrca_tree_depth_asexual(sroot, [0,1,2,3,4,5])
    assert depth == 0

    depth = phylodev.get_mrca_tree_depth_asexual(sroot, [2])
    assert depth == 1

    depth = phylodev.get_mrca_tree_depth_asexual(sroot, [5])
    assert depth == 2

if __name__ == "__main__":
    test_get_asexual_lineage_length()
    test_get_asexual_lineage_num_discrete_state_changes()
    test_get_asexual_lineage_num_discrete_unique_states()
    test_get_asexual_lineage_mutation_accumulation()
    test_get_mrca_tree_depth()