"""
Substructure matching using VF2 algorithm
"""
from typing import Dict, List, Set, Tuple, Optional
from .tokenizer import tokenize_molecule

def vf2_match(
    pattern_graph: Dict[str, List[str]],
    target_graph: Dict[str, List[str]],
    pattern_atoms: List[Dict],
    target_atoms: List[Dict],
    pattern_bonds: List[Dict],
    target_bonds: List[Dict]
) -> Optional[List[Tuple[str, str]]]:
    """
    VF2 algorithm for subgraph isomorphism
    
    Returns:
        List of (pattern_node, target_node) mappings if match found, None otherwise
    """
    pattern_nodes = list(pattern_graph.keys())
    target_nodes = list(target_graph.keys())
    
    if len(pattern_nodes) > len(target_nodes):
        return None
    
    # Build atom element mapping for constraint checking
    pattern_elements = {a["id"]: a["element"] for a in pattern_atoms}
    target_elements = {a["id"]: a["element"] for a in target_atoms}
    
    # Build bond order mapping
    pattern_bond_orders = {}
    for b in pattern_bonds:
        key = tuple(sorted([b["atom1"], b["atom2"]]))
        pattern_bond_orders[key] = b.get("order", 1)
    
    target_bond_orders = {}
    for b in target_bonds:
        key = tuple(sorted([b["atom1"], b["atom2"]]))
        target_bond_orders[key] = b.get("order", 1)
    
    # Simple matching: try all permutations (for small patterns)
    if len(pattern_nodes) <= 5:
        return _brute_force_match(
            pattern_graph, target_graph, pattern_nodes, target_nodes,
            pattern_elements, target_elements, pattern_bond_orders, target_bond_orders
        )
    
    # For larger patterns, use recursive backtracking
    return _recursive_match(
        pattern_graph, target_graph, pattern_nodes, target_nodes,
        pattern_elements, target_elements, pattern_bond_orders, target_bond_orders,
        {}, set(), set()
    )

def _brute_force_match(
    pattern_graph: Dict[str, List[str]],
    target_graph: Dict[str, List[str]],
    pattern_nodes: List[str],
    target_nodes: List[str],
    pattern_elements: Dict[str, str],
    target_elements: Dict[str, str],
    pattern_bond_orders: Dict[Tuple[str, str], int],
    target_bond_orders: Dict[Tuple[str, str], int]
) -> Optional[List[Tuple[str, str]]]:
    """Brute force matching for small patterns"""
    from itertools import permutations
    
    for perm in permutations(target_nodes, len(pattern_nodes)):
        mapping = dict(zip(pattern_nodes, perm))
        
        # Check element constraints
        valid = True
        for p_node, t_node in mapping.items():
            if pattern_elements.get(p_node) != target_elements.get(t_node):
                valid = False
                break
        
        if not valid:
            continue
        
        # Check edge constraints
        for p_node, p_neighbors in pattern_graph.items():
            t_node = mapping[p_node]
            t_neighbors = set(target_graph.get(t_node, []))
            
            for p_neighbor in p_neighbors:
                t_neighbor = mapping.get(p_neighbor)
                if t_neighbor not in t_neighbors:
                    valid = False
                    break
                
                # Check bond order
                p_key = tuple(sorted([p_node, p_neighbor]))
                t_key = tuple(sorted([t_node, t_neighbor]))
                if pattern_bond_orders.get(p_key, 1) > target_bond_orders.get(t_key, 1):
                    valid = False
                    break
            
            if not valid:
                break
        
        if valid:
            return list(mapping.items())
    
    return None

def _recursive_match(
    pattern_graph: Dict[str, List[str]],
    target_graph: Dict[str, List[str]],
    pattern_nodes: List[str],
    target_nodes: List[str],
    pattern_elements: Dict[str, str],
    target_elements: Dict[str, str],
    pattern_bond_orders: Dict[Tuple[str, str], int],
    target_bond_orders: Dict[Tuple[str, str], int],
    mapping: Dict[str, str],
    used_pattern: Set[str],
    used_target: Set[str]
) -> Optional[List[Tuple[str, str]]]:
    """Recursive backtracking for VF2"""
    if len(mapping) == len(pattern_nodes):
        return list(mapping.items())
    
    # Find next pattern node to match
    next_pattern = None
    for node in pattern_nodes:
        if node not in used_pattern:
            next_pattern = node
            break
    
    if next_pattern is None:
        return None
    
    # Try matching with each target node
    for target_node in target_nodes:
        if target_node in used_target:
            continue
        
        # Check element constraint
        if pattern_elements.get(next_pattern) != target_elements.get(target_node):
            continue
        
        # Check if mapping is consistent with existing edges
        valid = True
        for mapped_pattern, mapped_target in mapping.items():
            # Check if edges are preserved
            if mapped_pattern in pattern_graph.get(next_pattern, []):
                if target_node not in target_graph.get(mapped_target, []):
                    valid = False
                    break
                
                # Check bond order
                p_key = tuple(sorted([next_pattern, mapped_pattern]))
                t_key = tuple(sorted([target_node, mapped_target]))
                if pattern_bond_orders.get(p_key, 1) > target_bond_orders.get(t_key, 1):
                    valid = False
                    break
        
        if not valid:
            continue
        
        # Try this mapping
        mapping[next_pattern] = target_node
        used_pattern.add(next_pattern)
        used_target.add(target_node)
        
        result = _recursive_match(
            pattern_graph, target_graph, pattern_nodes, target_nodes,
            pattern_elements, target_elements, pattern_bond_orders, target_bond_orders,
            mapping, used_pattern, used_target
        )
        
        if result:
            return result
        
        # Backtrack
        del mapping[next_pattern]
        used_pattern.remove(next_pattern)
        used_target.remove(target_node)
    
    return None

def find_substructure(query_molecule: Dict[str, Any], target_molecule: Dict[str, Any]) -> Tuple[bool, Optional[List[Tuple[str, str]]]]:
    """
    Find if query is a substructure of target
    
    Returns:
        (is_match, atom_mapping)
    """
    query_tokenized = tokenize_molecule(query_molecule)
    target_tokenized = tokenize_molecule(target_molecule)
    
    # Quick fingerprint check
    if query_tokenized["fingerprint"] not in target_tokenized["fingerprint"]:
        # Simple check - if query fingerprint is substring, might match
        pass
    
    mapping = vf2_match(
        query_tokenized["graph"],
        target_tokenized["graph"],
        query_tokenized["atoms"],
        target_tokenized["atoms"],
        query_tokenized["bonds"],
        target_tokenized["bonds"]
    )
    
    return (mapping is not None, mapping)

