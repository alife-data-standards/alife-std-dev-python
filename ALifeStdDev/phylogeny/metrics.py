import networkx as nx
from . import utils

# ===== asexual lineage metrics =====

def get_asexual_lineage_length(lineage):
    """Get asexual lineage length.

    Will check that given lineage is an asexual lineage.

    Args:
        lineage (networkx.DiGraph): an asexual lineage

    Returns:
        length (int) of given lineage
    """
    if not utils.is_asexual_lineage(lineage): raise Exception("the given lineage is not an asexual lineage")
    return len(lineage.nodes)

def get_asexual_lineage_num_discrete_state_changes(lineage, attribute_list):
    """Get the number of discrete state changes from an asexual lineage.

    State is described by the aggregation of all attributes give by attribute list.

    Args:
        lineage (networkx.DiGraph): an asexual lineage
        attribute_list (list): list of attributes (strings) to use when defining
            a state

    Returns:
        Returns the number of discrete states along the lineage.
    """
    # Check that lineage is an asexual lineage.
    if not utils.is_asexual_lineage(lineage): raise Exception("the given lineage is not an asexual lineage")
    # Check that all nodes have all given attributes in the attribute list
    if not utils.all_taxa_have_attributes(lineage, attribute_list): raise Exception("given attributes are not universal among all taxa along the lineage")
    # get the first state (root node)
    lineage_id = utils.get_root_ids(lineage)[0]
    num_states = 1
    cur_state = [lineage.nodes[lineage_id][attr] for attr in attribute_list]
    # count the number of state changes moving down the lineage
    while True:
        successor_ids = list(lineage.successors(lineage_id))
        if len(successor_ids) == 0: break # We've hit the last thing!
        lineage_id = successor_ids[0]
        state = [lineage.nodes[lineage_id][attr] for attr in attribute_list]
        if cur_state != state:
            cur_state = state
            num_states += 1
    return num_states

def get_asexual_lineage_num_discrete_unique_states(lineage, attribute_list):
    """Get the number of discrete unique states along a lineage where what it
    means to be a state is defined by attribute_list.

    Args:
        lineage (networkx.DiGraph): an asexual lineage
        attribute_list (list): list of attributes (strings) to use when defining
            a state

    Returns:
        The number of discrete unique states found along the lineage.
    """
    # Check that lineage is an asexual lineage.
    if not utils.is_asexual_lineage(lineage): raise Exception("the given lineage is not an asexual lineage")
    # Check that all nodes have all given attributes in the attribute list
    if not utils.all_taxa_have_attributes(lineage, attribute_list): raise Exception("given attributes are not universal among all taxa along the lineage")
    # get the first state (root node)
    lineage_id = utils.get_root_ids(lineage)[0]
    unique_states = set()
    unique_states.add(tuple([lineage.nodes[lineage_id][attr] for attr in attribute_list]))
    while True:
        successor_ids = list(lineage.successors(lineage_id))
        if len(successor_ids) == 0: break # We've hit the last thing!
        lineage_id = successor_ids[0]
        unique_states.add(tuple([lineage.nodes[lineage_id][attr] for attr in attribute_list]))
    return len(unique_states)

def get_asexual_lineage_mutation_accumulation(lineage, mutation_attributes, skip_root=False):
    """Get the distribution of mutation type accumulations over an asexual lineage.

    Args:
        lineage (networkx.DiGraph): an asexual lineage
        mutation_attributes (list of str): what are the mutation count attributes
            that we should accumulate over the lineage?
        skip_root (bool): Should we include root node mutation count values in
            our accumlation? Defaults to false.

    Returns:
        A dictionary indexed by mutation types (mutation_attributes) where each
        value in the dictionary is the sum of that type of mutation along the lineage.
    """
    # Check that lineage is an asexual lineage.
    if not utils.is_asexual_lineage(lineage): raise Exception("the given lineage is not an asexual lineage")
    # Check that all nodes have all given attributes in the attribute list
    if not utils.all_taxa_have_attributes(lineage, mutation_attributes): raise Exception("given mutation attributes are not universal among all taxa along the lineage")
    # initialize
    mut_accumulators = {mut_attr:0 for mut_attr in mutation_attributes}
    # get the root node
    lineage_id = utils.get_root_ids(lineage)[0]
    if not skip_root:
        for mut_attr in mutation_attributes:
            mut_accumulators[mut_attr] += lineage.nodes[lineage_id][mut_attr]
    while True:
        successor_ids = list(lineage.successors(lineage_id))
        if len(successor_ids) == 0: break # We've hit the last thing!
        # Is this a new state or a member of the current state?
        lineage_id = successor_ids[0]
        for mut_attr in mutation_attributes:
            mut_accumulators[mut_attr] += lineage.nodes[lineage_id][mut_attr]
    return mut_accumulators

# ===== asexual phylogeny metrics =====

def get_mrca_tree_depth_asexual(phylogeny, ids=None):
    """Get the tree depth of the most recent common ancestor shared by the specified
    taxa ids (ids) in an asexual phylogeny (phylogeny).
    """
    # Get the id of the most recent common ancestor
    mrca_id = utils.get_mrca_id_asexual(phylogeny, ids)
    if mrca_id == -1: raise Exception("phylogeny has no common ancestor")
    # Calculate distance from root to mrca
    cur_id = mrca_id
    depth = 0
    while True:
        ancestor_ids = list(phylogeny.predecessors(cur_id))
        if len(ancestor_ids) == 0: break
        depth+=1
        cur_id = ancestor_ids[0]
    return depth

# ===== phylogenetic richness =====

def calc_phylogenetic_diversity_asexual(phylogeny, ids=None):
    """Calculate phylogenetic diversity (i.e., the number of nodes in the minimum
    spanning tree from the MRCA to all extant taxa). Currently only for asexual
    phylogenies.

    (Faith, 1992)

    ids gives the set we want to calculate phylogenetic diversity on. i.e.,
    we'll get the mrca for those ids and compute the minimum spanning tree

    none defaults to including all leaf nodes
    """
    # if given no ids, default to leaf taxa; otherwise, validate given ids
    if ids == None:
        # Find MRCA on leaf nodes
        ids = utils.get_leaf_taxa_ids(phylogeny)
    # (1) get the mrca
    mrca_id = utils.get_mrca_id_asexual(phylogeny, ids)
    if mrca_id == -1: raise Exception("given ids have no common ancestor")
    # (2) collect paths from each id to mrca
    canopy = set([i for i in ids] + [mrca_id])
    for i in ids:
        cur_id = i
        while True:
            ancestor_ids = list(phylogeny.predecessors(cur_id))
            if len(ancestor_ids) == 0: break
            cur_id = ancestor_ids[0]
            # If we've encountered this path before, we can skip the rest because
            # we're guaranteed an asexual phylogeny.
            if cur_id in canopy: break
            canopy.add(cur_id)
    # Build a subgraph with only the canopy
    canopy_phylo = nx.subgraph(phylogeny, list(canopy))
    # Okay, now we can compute the minimum spanning tree.
    return len(nx.minimum_spanning_tree(canopy_phylo.to_undirected()).nodes)