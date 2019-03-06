from . import utils

# ===== asexual lineage metrics =====

def get_asexual_lineage_length(lineage):
    """Get asexual lineage length.

    Will check that given lineage is an asexual lineage.

    """
    if not utils.is_asexual_lineage(lineage): raise Exception("the given lineage is not an asexual lineage")
    return len(lineage.nodes)

def get_asexual_lineage_num_discrete_state_changes(lineage, attribute_list):
    """Get the number of discrete state changes from an asexual lineage.

    State is described by the aggregation of all attributes give by attribute list.

    Args: todo

    Returns: todo
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
    """
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
    """
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
