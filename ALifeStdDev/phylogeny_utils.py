import networkx as nx

# ===== Verification =====
def all_taxa_have_attribute(phylogeny, attribute):
    """
    Do all taxa in phylogeny have given attribute?
    """
    for node in phylogeny.nodes:
        if not (attribute in phylogeny.nodes[node]): return False
    return True

def is_asexual(phylogeny):
    """
    is this phylogeny asexual?
    """
    for node in phylogeny.nodes:
        if len(list(phylogeny.predecessors(node))) > 1: return False
    return True


# ===== Rootedness-related utilities =====

def has_single_root(phylogeny):
    """
    Given phylogeny, return True if it has only a single root and False if it has
    mulitple roots.
    (wrapper around nx.is_connected)
    """
    return nx.is_weakly_connected(phylogeny)

def get_root_ids(phylogeny):
    """
    get root ids of phylogeny
    """
    return [node for node in phylogeny.nodes if len(list(phylogeny.predecessors(node))) == 0]

def get_roots(phylogeny):
    roots = {node:phylogeny.nodes[node] for node in phylogeny.nodes if len(list(phylogeny.predecessors(node))) == 0}
    for r in roots: roots[r]["id"] = r
    return roots

def get_num_roots(phylogeny):
    """
    Given a phylogeny (that may contain multiple roots), return number of roots.
    """
    return nx.number_weakly_connected_components(phylogeny)

def get_independent_phylogenies(phylogeny):
    """
    Given phylogeny (that may contain multiple roots), return independently rooted phylogenies
    (deep copy)
    """
    components = [c for c in sorted(nx.weakly_connected_components(phylogeny), key=len, reverse=True)]
    # print("Components: {}".format(components))
    phylogenies = [phylogeny.subgraph(comp).copy() for comp in components]
    return phylogenies

# ===== Extracting the extant taxa =====

def get_extant_taxa_from_pruned(phylogeny):
    """
    given a phylogeny, return extant population details (as a dictionary) indexed
    by id
    """
    extant = {node:phylogeny.nodes[node] for node in phylogeny.nodes if len(list(phylogeny.successors(node))) == 0}
    for e in extant: extant[e]["id"] = e
    return extant

def get_extant_taxa_ids_from_pruned(phylogeny):
    """
    given a phylogeny, return ids of extant taxa
    """
    extant_ids = [node for node in phylogeny.nodes if len(list(phylogeny.successors(node))) == 0]
    return extant_ids

def get_extant_taxa_ids_by_destruction_time(phylogeny, not_destroyed_value="none"):
    """
    given a phylogeny, return ids of extant taxa (i.e., taxa where destruction_time is equal
    to given not_destroyed_value)
    """
    # Check if all taxa have destruction time attribute
    if (not all_taxa_have_attribute(phylogeny, "destruction_time")): raise Exception(f"Not all taxa have 'destruction_time' attribute")
    extant_ids = [node for node in phylogeny.nodes if phylogeny.nodes[node]["destruction_time"]==not_destroyed_value]
    return extant_ids

def get_extant_taxa_by_destruction_time(phylogeny, not_destroyed_value="none"):
    """
    given a phylogeny, return extant population details (as a dictionary) indexed
    by id
    """
    # Check if all taxa have destruction time attribute
    if (not all_taxa_have_attribute(phylogeny, "destruction_time")): raise Exception(f"Not all taxa have 'destruction_time' attribute")
    extant = {node:phylogeny.nodes[node] for node in phylogeny.nodes if phylogeny.nodes[node]["destruction_time"]==not_destroyed_value}
    for e in extant: extant[e]["id"] = e
    return extant

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
