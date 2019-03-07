import networkx as nx

# ===== Verification =====
def all_taxa_have_attribute(phylogeny, attribute):
    """Do all taxa in the given phylogeny have the given attribute?

    Args:
        phylogeny (networkx.DiGraph): graph object that describes a phylogeny
        attribute (str): a possible attribute/descriptor for a taxa (node) in phylogeny

    Returns:
        True if all taxa (nodes) in the phylogeny have the given attribute and False
        otherwise.
    """
    for node in phylogeny.nodes:
        if not (attribute in phylogeny.nodes[node]): return False
    return True

def all_taxa_have_attributes(phylogeny, attribute_list):
    """Do all taxa in the given phylogeny have the given attributes?

    Args:
        phylogeny (networkx.DiGraph): graph object that describes a phylogeny
        attribute (str): a list of attributes to check for in given phylogeny

    Returns:
        True if all taxa (nodes) in the phylogeny have the all of the given attributes and False
        otherwise.
    """
    for node in phylogeny.nodes:
        for attribute in attribute_list:
            if not (attribute in phylogeny.nodes[node]): return False
    return True

def is_asexual(phylogeny):
    """Is this an asexual phylogeny?

    A phylogeny is considered to be asexual if all taxa (nodes) have a single direct
    ancestor (predecessor).

    Args:
        phylogeny (networkx.DiGraph): graph object that describes a phylogeny

    Returns:
        True if the phylogeny is asexual and False otherwise.
    """
    for node in phylogeny.nodes:
        if len(list(phylogeny.predecessors(node))) > 1: return False
    return True

def is_asexual_lineage(phylogeny):
    """Does this phylogeny give an asexual lineage?

    To be an asexual lineage, all internal nodes in the phylogeny have exactly
    one predecessor and exactly one successor.

    Args:
        phylogeny (networkx.DiGraph): graph object that describes a phylogeny

    Returns:
        True if the phylogeny is an asexual lineage and False otherwise.
    """
    lineage_ids = get_root_ids(phylogeny)
    # There should only be a single root if the given phylogeny is a single, asexual
    # lineage
    if len(lineage_ids) != 1: return False
    while True:
        successor_ids = list(phylogeny.successors(lineage_ids[-1]))
        if len(successor_ids) > 1: return False
        if len(successor_ids) == 0: break
        lineage_ids.append(successor_ids[0])
    return True

# ===== Rootedness-related utilities =====

def has_single_root(phylogeny):
    """Given phylogeny, return True if it has only a single root and False if it has
    mulitple roots.

    This function just wraps the networkx is_weekly_connected function.

    Args:
        phylogeny (networkx.DiGraph): graph object that describes a phylogeny

    Returns:
        True if it has only a single root and False if it has mulitple roots.
    """
    return nx.is_weakly_connected(phylogeny)

def get_root_ids(phylogeny):
    """Get ids of root nodes in phylogeny

    Args:
        phylogeny (networkx.DiGraph): graph object that describes a phylogeny

    Returns:
        For all nodes in phylogeny, return ids of nodes with no predecessors.
    """
    return [node for node in phylogeny.nodes if len(list(phylogeny.predecessors(node))) == 0]

def get_roots(phylogeny):
    """Get root nodes in phylogeny (does not assume that the given phylogeny has a single root).

    Args:
        phylogeny (networkx.DiGraph): graph object that describes a phylogeny

    Returns:
        For all nodes in phylogeny, return dictionary of root nodes (nodes with no predecessors).
        The returned dictionary is keyed by node ids.
        Each node in the returned list is a dictionary with all of the node's descriptors/attributes.
    """
    roots = {node:phylogeny.nodes[node] for node in phylogeny.nodes if len(list(phylogeny.predecessors(node))) == 0}
    for r in roots: roots[r]["id"] = r
    return roots

def get_num_roots(phylogeny):
    """Given a phylogeny (that may contain multiple roots), return number of roots
    where a root is a node with no predecessors.

    Args:
        phylogeny (networkx.DiGraph): graph object that describes a phylogeny

    Returns:
        Returns the number of independent trees (i.e., roots) in the given phylogeny.
    """
    return len(get_root_ids(phylogeny))

def get_num_independent_phylogenies(phylogeny):
    """Get number of the independently-rooted trees within the given phylogeny.

    This function wraps networkx's number_weakly_connected_components function.

    Args:
        phylogeny (networkx.DiGraph): graph object that describes a phylogeny

    Returns:
        Returns the number of weakly connected components (independent trees) in
        the given phylogeny.
    """
    return nx.number_weakly_connected_components(phylogeny)

def get_independent_phylogenies(phylogeny):
    """Get a list of the independently-rooted trees within the given phylogeny.

    Args:
        phylogeny (networkx.DiGraph): graph object that describes a phylogeny

    Returns:
        Returns a list of networkx.DiGraph objects.
        Each member of the returned list is an independent (not connected) subgraph
        of the given phylogeny. The returned list of networkx.DiGraph objects are
        copies.
    """
    components = [c for c in sorted(nx.weakly_connected_components(phylogeny), key=len, reverse=True)]
    phylogenies = [phylogeny.subgraph(comp).copy() for comp in components]
    return phylogenies

# ===== Extracting the extant taxa =====


def get_leaf_taxa(phylogeny):
    """Get the leaf taxa (taxa with no successors/descendants) of the given phylogeny.

    Args:
        phylogeny (networkx.DiGraph): graph object that describes a phylogeny

    Returns:
        Returns dictionary of leaf taxa nodes.
        The returned dictionary is keyed by node ids.
        Each node in the returned list is a dictionary with all of the node's descriptors/attributes.
    """
    extant = {node:phylogeny.nodes[node] for node in phylogeny.nodes if len(list(phylogeny.successors(node))) == 0}
    for e in extant: extant[e]["id"] = e
    return extant


def get_leaf_taxa_ids(phylogeny):
    """Given a phylogeny, return list of leaf taxa (taxa with no successors/descendants)

    Args:
        phylogeny (networkx.DiGraph): graph object that describes a phylogeny

    Returns:
        For all nodes in phylogeny, return ids of nodes with no successors (descendants).
    """
    extant_ids = [node for node in phylogeny.nodes if len(list(phylogeny.successors(node))) == 0]
    return extant_ids


def get_extant_taxa_ids(phylogeny, time="present", not_destroyed_value="none",
                        destruction_attribute="destruction_time",
                        origin_attribute="origin_time"):
    """
    Get ids of extant taxa from a phylogeny

    Args:
        phylogeny (networkx.DiGraph): graph object that describes a phylogeny
        not_destroyed_value (str): value of taxa[attribute] that indicates that
            the taxa is not destroyed (i.e., still exists)
        attribute (str): attribute to use to determine if taxa still exists

    Returns:
        List of extant taxa ids.
    """
    # Check if all taxa have destruction time attribute
    validate_destruction_time(phylogeny, destruction_attribute)

    # Check if all taxa have origin time attribute
    if (time != "present"):
        validate_origin_time(phylogeny, origin_attribute)

    extant_ids = [node for node in phylogeny.nodes
                  if taxon_is_alive(phylogeny.nodes[node],
                                    time=time,
                                    not_destroyed_value=not_destroyed_value,
                                    destruction_attribute=destruction_attribute,
                                    origin_attribute=origin_attribute)]
    return extant_ids


def get_extant_taxa(phylogeny, time="present", not_destroyed_value="none",
                    destruction_attribute="destruction_time",
                    origin_attribute="origin_time"):
    """
    Get extant taxa from a phylogeny

    Args:
        phylogeny (networkx.DiGraph): graph object that describes a phylogeny
        not_destroyed_value (str): value of taxa[attribute] that indicates that
            the taxa is not destroyed (i.e., still exists)
        attribute (str): attribute to use to determine if taxa still exists

    Returns:
        Returns dictionary of extant taxa.
        The returned dictionary is keyed by node ids.
        Each node in the returned list is a dictionary with all of the node's
        descriptors/attributes.
    """
    # Check if all taxa have destruction time attribute
    validate_destruction_time(phylogeny, destruction_attribute)

    # Check if all taxa have origin time attribute
    if (time != "present"):
        validate_origin_time(phylogeny, origin_attribute)

    extant = {node: phylogeny.nodes[node] for node in phylogeny.nodes
              if taxon_is_alive(phylogeny.nodes[node],
                                time=time,
                                not_destroyed_value=not_destroyed_value,
                                destruction_attribute=destruction_attribute,
                                origin_attribute=origin_attribute)}
    for e in extant:
        extant[e]["id"] = e
    return extant


def taxon_is_alive(node, time, not_destroyed_value="none",
                   destruction_attribute="destruction_time",
                   origin_attribute="origin_time"):
    return (node["destruction_time"] == not_destroyed_value  # not dead yet
            or node[destruction_attribute] > time) \
            and (time == "present" or
                 node[origin_attribute] < time)  # has been born


def validate_destruction_time(phylogeny, attribute="destruction_time"):
    if (not all_taxa_have_attribute(phylogeny, attribute)):
        raise Exception(f"Not all taxa have '{attribute}' data")


def validate_origin_time(phylogeny, attribute="origin_time"):
    if (not all_taxa_have_attribute(phylogeny, attribute)):
        raise Exception(f"Not all taxa have '{attribute}' data")


# ===== lineages-specific utilities =====

def extract_asexual_lineage(phylogeny, taxa_id):
    """Given a phylogeny, extract the ancestral lineage of the taxa specified by
    taxa_id. Only works for asexual phylogenies.

    Args:
        phylogeny (networkx.DiGraph): graph object that describes a phylogeny
        taxa_id (int): id of taxa to extract an ancestral lineage for (must be
            a valid node id in the given phylogeny).

    Returns:
        networkx.DiGraph that contains the ancestral lineage of the specified taxa.
    """
    return phylogeny.subgraph(extract_asexual_lineage_ids(phylogeny, taxa_id)).copy()

def extract_asexual_lineage_ids(phylogeny, taxa_id):
    """Given a phylogeny, extract the ids of the members of the ancestral lineage
        of the taxa specified by taxa_id. Only works for asexual phylogenies.

    Args:
        phylogeny (networkx.DiGraph): graph object that describes a phylogeny
        taxa_id (int): id of taxa to extract an ancestral lineage for (must be
            a valid node id in the given phylogeny).

    Returns:
        List of ids along specified taxa's lineage.
    """
    # Make sure taxa id is in the phylogeny
    if not taxa_id in phylogeny.nodes: raise Exception(f"Failed to find given taxa ({taxa_id}) in phylogeny")
    if not is_asexual(phylogeny): raise Exception("Given phylogeny is not asexual")
    # Get taxa ids on lineage
    ids_on_lineage = [taxa_id]
    while True:
        ancestor_ids = list(phylogeny.predecessors(ids_on_lineage[-1]))
        if len(ancestor_ids) == 0: break
        ids_on_lineage.append(ancestor_ids[0])
    return ids_on_lineage

def abstract_asexual_lineage(lineage, attribute_list, origin_time_attr="origin_time", destruction_time_attr="destruction_time"):
    """Given an asexual lineage, abstract as sequence of states where state-ness
    is described by attributes. I.e., compress the lineage into a sequence of states.

    Args:
        lineage (networkx.DiGraph): an asexual lineage
        attribute_list (list of strings): a list of attributes to use to define
            taxa state
        origin_time_attr (str): attribute key for origin_time of a taxa
        destruction_time_attr (str): attribute key for destruction time of a taxa

    Returns:
        networkx.DiGraph objects that describes an abstracted version of the given
        lineage.
    """
    # Check that lineage is an asexual lineage.
    if not is_asexual_lineage(lineage): raise Exception("the given lineage is not an asexual lineage")
    # Check that all nodes have all given attributes in the attribute list
    if not all_taxa_have_attributes(lineage, attribute_list): raise Exception("given attributes are not universal among all taxa along the lineage")
    # Make sure we prevent collisions between given attributes and attributes that are used to document states
    if "node_state" in attribute_list:
        raise Exception("'node_state' is a reserved attribute when using this function")
    if "members" in attribute_list:
        raise Exception("'members' is a reserved attribute when using this function")
    if "state_id" in attribute_list:
        raise Exception("'state_id' is a reserved attribute when using this function")

    track_origin = all_taxa_have_attribute(lineage, origin_time_attr)
    track_destruction = all_taxa_have_attribute(lineage, destruction_time_attr)

    abstract_lineage = nx.DiGraph() # Empty graph to hold our abstract lineage.

    # Start with the root node (to qualify as a lineage, graph must only have one root id)
    root_id = get_root_ids(lineage)[0]
    state_id = 0
    # Add the first lineage state to the abstract lineage
    abstract_lineage.add_node(state_id)
    abstract_lineage.nodes[state_id]["state_id"] = state_id
    abstract_lineage.nodes[state_id]["node_state"] = [lineage.nodes[root_id][attr] for attr in attribute_list]
    # Add attributes to state
    for attr in attribute_list:
        abstract_lineage.nodes[state_id][attr] = lineage.nodes[root_id][attr]
    if track_origin:
        abstract_lineage.nodes[state_id]["origin_time"] = lineage.nodes[root_id][origin_time_attr]
    if track_destruction: # (this might get updated as we go)
        abstract_lineage.nodes[state_id]["destruction_time"] = lineage.nodes[root_id][destruction_time_attr]
    # Add first member
    abstract_lineage.nodes[state_id]["members"] = {root_id:lineage.nodes[root_id]}

    lineage_id = root_id
    while True:
        successor_ids = list(lineage.successors(lineage_id))
        if len(successor_ids) == 0: break # We've hit the last thing!
        # Is this a new state or a member of the current state?
        lineage_id = successor_ids[0]
        state = [lineage.nodes[lineage_id][attr] for attr in attribute_list]
        if abstract_lineage.nodes[state_id]["node_state"] == state:
            # Add this taxa as member of current state
            # - update time of destruction etc
            abstract_lineage.nodes[state_id]["members"][lineage_id] = lineage.nodes[lineage_id]
            if track_destruction:
                abstract_lineage.nodes[state_id]["destruction_time"] = lineage.nodes[lineage_id][destruction_time_attr]
        else:
            # Add new state
            state_id += 1
            abstract_lineage.add_node(state_id)
            abstract_lineage.add_edge(state_id-1, state_id)
            # Document state information
            abstract_lineage.nodes[state_id]["state_id"] = state_id
            abstract_lineage.nodes[state_id]["node_state"] = [lineage.nodes[lineage_id][attr] for attr in attribute_list]
            if track_origin:
                abstract_lineage.nodes[state_id]["origin_time"] = lineage.nodes[lineage_id][origin_time_attr]
            for attr in attribute_list:
                abstract_lineage.nodes[state_id][attr] = lineage.nodes[lineage_id][attr]
            # Add first member
            abstract_lineage.nodes[state_id]["members"] = {lineage_id:lineage.nodes[lineage_id]}

    return abstract_lineage

# ===== mrca =====

def has_common_ancestor_asexual(phylogeny, ids=None):
    """Do the given set of ids share a common ancestor in the given phylogeny?
    """
    # check that phylogeny is asexual
    if not is_asexual(phylogeny): raise Exception("given phylogeny is not asexual")
    # if given no ids, default to leaf taxa; otherwise, validate given ids
    if ids == None:
        # Find MRCA on leaf nodes
        ids = get_leaf_taxa_ids(phylogeny)
    else:
        # Check validity of ids
        ids = list(set(ids)) # eliminate duplicates
        for i in ids:
            if not i in phylogeny.nodes: raise Exception(f"failed to find {i} in phylogeny")

    # Check for case where all ids are the same (return that id) and the case where we have no ids (return -1)
    if len(set(ids)) == 1: return True
    elif len(ids) == 0: return False

    # Get the lineages of each taxa
    lineages = [set(extract_asexual_lineage_ids(phylogeny,i)) for i in ids]
    common_ancestors = set.intersection(*lineages)

    if len(common_ancestors) > 0: return True
    else: return False # No common ancestors

def get_mrca_id_asexual(phylogeny, ids=None):
    """Get the id of the most recent common ancestor (mrca) shared among taxa specified
    by ids. Returns -1 if no mrca exists.
    """
    # check that phylogeny is asexual
    if not is_asexual(phylogeny): raise Exception("given phylogeny is not asexual")
    # if given no ids, default to leaf taxa; otherwise, validate given ids
    if ids == None:
        # Find MRCA on leaf nodes
        ids = get_leaf_taxa_ids(phylogeny)
    else:
        # Check validity of ids
        ids = list(set(ids)) # eliminate duplicates
        for i in ids:
            if not i in phylogeny.nodes: raise Exception(f"failed to find {i} in phylogeny")

    # Check for case where all ids are the same (return that id) and the case where we have no ids (return -1)
    if len(set(ids)) == 1: return ids[0]
    elif len(ids) == 0: return -1

    # Get the lineages of each taxa
    lineages = [set(extract_asexual_lineage_ids(phylogeny,i)) for i in ids]
    common_ancestors = set.intersection(*lineages)

    if len(common_ancestors) == 1: return list(common_ancestors)[0]
    elif len(common_ancestors) == 0: return -1 # No common ancestors

    # Multiple common ancestors, need to figure out which is most recent.
    # - We can assume that there _is_ a common ancestor, so we know we're working
    #   with a single tree.
    # The most recent common ancestor will be the ancestor most closely related
    # to any of the given ids.
    cur_taxa = ids[0]
    while True:
        # There are multiple common ancestors, which means each lineage length is > 1, so this [0] should never fail.
        if cur_taxa in common_ancestors: break
        cur_taxa = list(phylogeny.predecessors(cur_taxa))[0]
    return cur_taxa

def get_mrca_asexual(phylogeny, ids=None):
    """Get the most recent common ancestor (mrca) shared among taxa specified
    by ids. Returns None if no mrca exists.
    """
    if has_common_ancestor_asexual(phylogeny, ids):
        mrca_id = get_mrca_id_asexual(phylogeny, ids)
        mrca = phylogeny.nodes[mrca_id]
        mrca["id"] = mrca_id
        return mrca
    else:
        return None

# ===== miscellaneous =====

def get_pairwise_distances(phylogeny, ids):
    """
    given phylogeny and some ids to compute the pairwise distances between, return
    pairwise distances

    Raises:
        networkx.exception.NetworkXNoPath: if no path between any two of the given ids
    """
    dists = []
    undirected_phylo = phylogeny.to_undirected()
    for i in range(0, len(ids)):
        for j in range(i+1, len(ids)):
            dists.append(nx.shortest_path_length(undirected_phylo, source=ids[i], target=ids[j]))
            # print(f"Dist between: {ids[i]}-{ids[j]} = {dists[-1]}")
    return dists