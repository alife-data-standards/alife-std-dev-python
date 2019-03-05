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
        The returned dictionary is keyed by
        Each node in the returned list is a dictionary with all of the node's descriptors/attributes.
    """
    roots = {node:phylogeny.nodes[node] for node in phylogeny.nodes if len(list(phylogeny.predecessors(node))) == 0}
    for r in roots: roots[r]["id"] = r
    return roots

def get_num_roots(phylogeny):
    """Given a phylogeny (that may contain multiple roots), return number of roots.

    This function is a wrapper around networkx's number_weekly_connected_components
    function.

    Args:
        phylogeny (networkx.DiGraph): graph object that describes a phylogeny

    Returns:
        Returns the number of independent trees (i.e., roots) in the given phylogeny.
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

    """
    extant = {node:phylogeny.nodes[node] for node in phylogeny.nodes if len(list(phylogeny.successors(node))) == 0}
    for e in extant: extant[e]["id"] = e
    return extant


def get_leaf_taxa_ids(phylogeny):
    """
    given a phylogeny, return ids of extant taxa
    """
    extant_ids = [node for node in phylogeny.nodes if len(list(phylogeny.successors(node))) == 0]
    return extant_ids


def get_extant_taxa_ids(phylogeny, time="present", not_destroyed_value="none"):
    """
    given a phylogeny, return ids of extant taxa (i.e., taxa where 
    destruction_time is equal to given not_destroyed_value)
    """
    # Check if all taxa have destruction time attribute
    validate_destruction_time(phylogeny)

    # Check if all taxa have origin time attribute
    if (time != "present"):
        validate_origin_time(phylogeny)

    extant_ids = [node for node in phylogeny.nodes
                  if taxon_is_alive(phylogeny.nodes[node],
                                    time, not_destroyed_value)]
    return extant_ids


def get_extant_taxa(phylogeny, time="present", not_destroyed_value="none"):
    """
    given a phylogeny, return extant population details (as a dictionary)
    indexed by id
    """
    # Check if all taxa have destruction time attribute
    validate_destruction_time(phylogeny)

    # Check if all taxa have origin time attribute
    if (time != "present"):
        validate_origin_time(phylogeny)

    extant = {node: phylogeny.nodes[node] for node in phylogeny.nodes
              if taxon_is_alive(phylogeny.nodes[node],
                                time, not_destroyed_value)}
    for e in extant:
        extant[e]["id"] = e
    return extant


def taxon_is_alive(node, time, not_destroyed_value="none"):
    return (node["destruction_time"] == not_destroyed_value  # not dead yet
            or node["destruction_time"] > time) \
            and (time == "present" or
                 node["origin_time"] < time)  # has been born


def validate_destruction_time(phylogeny):
    if (not all_taxa_have_attribute(phylogeny, "destruction_time")):
        raise Exception(f"Not all taxa have 'destruction_time' attribute")


def validate_origin_time(phylogeny):
    if (not all_taxa_have_attribute(phylogeny, "origin_time")):
        raise Exception(f"Not all taxa have 'origin_time' attribute")


# ===== Extracting lineages =====

def extract_asexual_lineage(phylogeny, taxa_id):
    """
    Given a phylogeny, extract lineage of given id
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
    return phylogeny.subgraph(ids_on_lineage).copy()
