"""
Feature extraction for molecules - matches frontend featurizer

Ensures consistent atom ordering between frontend and backend.
"""

from typing import Dict, List, Optional, Tuple
import numpy as np
import logging

# Optional PyTorch imports - delay import to avoid DLL issues at module load time
TORCH_AVAILABLE = None
torch = None
Data = None

def _ensure_torch():
    """Lazy import of PyTorch to avoid DLL issues at module load."""
    global TORCH_AVAILABLE, torch, Data
    if TORCH_AVAILABLE is not None:
        return TORCH_AVAILABLE
    
    try:
        import torch as _torch
        from torch_geometric.data import Data as _Data
        torch = _torch
        Data = _Data
        TORCH_AVAILABLE = True
        return True
    except (ImportError, OSError) as e:
        TORCH_AVAILABLE = False
        torch = None
        Data = None
        logging.warning(f"PyTorch/PyG not available: {e}. Featurization will use numpy only.")
        return False

from rdkit import Chem
from rdkit.Chem import Descriptors


# Element to atomic number mapping (matches frontend)
ELEMENT_TO_ATOMIC_NUMBER = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8,
    'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15,
    'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Fe': 26, 'Cu': 29,
    'Zn': 30, 'Br': 35, 'I': 53,
}

# Valence table (matches frontend)
VALENCE_TABLE = {
    'H': 1, 'C': 4, 'N': 3, 'O': 2, 'F': 1, 'S': 6, 'Cl': 1, 'Br': 1,
    'I': 1, 'P': 5, 'B': 3, 'Si': 4,
}


def compute_node_features(atom, mol, bonds) -> List[float]:
    """
    Compute node features for an atom.
    Matches frontend Featurizer.computeNodeFeatures()
    """
    features = []
    
    # Atomic number
    atomic_num = ELEMENT_TO_ATOMIC_NUMBER.get(atom.GetSymbol(), 0)
    features.append(float(atomic_num))
    
    # Degree (number of bonds)
    degree = atom.GetDegree()
    features.append(float(degree))
    
    # Bond order sum
    bond_order_sum = sum(bond.GetBondTypeAsDouble() for bond in bonds)
    features.append(float(bond_order_sum))
    
    # Formal charge
    formal_charge = atom.GetFormalCharge()
    features.append(float(formal_charge))
    
    # Hybridization estimate
    # 0 = sp, 1 = sp2, 2 = sp3
    hybridization = 2  # default sp3
    if bond_order_sum <= 2:
        hybridization = 0  # sp
    elif bond_order_sum <= 3:
        hybridization = 1  # sp2
    features.append(float(hybridization))
    
    # Valence electrons
    valence = VALENCE_TABLE.get(atom.GetSymbol(), 4)
    features.append(float(valence))
    
    # Position (normalized) - will be 0 if not available
    pos = mol.GetConformer().GetAtomPosition(atom.GetIdx()) if mol.GetNumConformers() > 0 else (0, 0, 0)
    features.append(pos[0] / 100.0)  # Normalize
    features.append(pos[1] / 100.0)
    features.append(pos[2] / 100.0)
    
    # Additional RDKit features
    features.append(float(atom.GetIsAromatic()))
    features.append(float(atom.GetHybridization()))
    features.append(float(atom.GetNumRadicalElectrons()))
    
    return features


def compute_edge_features(bond) -> List[float]:
    """
    Compute edge features for a bond.
    Matches frontend Featurizer.computeEdgeFeatures()
    """
    features = []
    
    # Bond order
    bond_order = bond.GetBondTypeAsDouble()
    features.append(float(bond_order))
    
    # Bond type encoding
    features.append(1.0 if bond_order == 1 else 0.0)  # Single
    features.append(1.0 if bond_order == 2 else 0.0)  # Double
    features.append(1.0 if bond_order == 3 else 0.0)  # Triple
    
    # Additional features
    features.append(float(bond.GetIsAromatic()))
    features.append(float(bond.GetIsConjugated()))
    features.append(float(bond.GetStereo()))
    
    return features


def featurize_smiles(
    smiles: str,
    atom_order: Optional[List[str]] = None,
    sanitize: bool = True
) -> Tuple[Data, Dict[int, str]]:
    """
    Featurize SMILES string to PyTorch Geometric Data object.
    
    Args:
        smiles: SMILES string
        atom_order: Optional list of frontend atom IDs to preserve order
        sanitize: Whether to sanitize molecule with RDKit
    
    Returns:
        Data: PyG Data object
        node_mapping: Dict mapping node index to atom ID (if atom_order provided)
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        
        if sanitize:
            try:
                Chem.SanitizeMol(mol)
            except:
                # Try to fix common issues
                mol = Chem.MolFromSmiles(smiles, sanitize=False)
                if mol is None:
                    raise ValueError("Molecule sanitization failed")
        
        # Add hydrogens
        mol = Chem.AddHs(mol)
        
        # Generate 2D coordinates if not present
        if mol.GetNumConformers() == 0:
            from rdkit.Chem import AllChem
            AllChem.Compute2DCoords(mol)
        
    except Exception as e:
        raise ValueError(f"Failed to parse SMILES: {str(e)}")
    
    # Node features
    node_features = []
    for atom in mol.GetAtoms():
        bonds = [mol.GetBondWithIdx(i) for i in range(mol.GetNumBonds())
                 if mol.GetBondWithIdx(i).GetBeginAtomIdx() == atom.GetIdx() or
                 mol.GetBondWithIdx(i).GetEndAtomIdx() == atom.GetIdx()]
        features = compute_node_features(atom, mol, bonds)
        node_features.append(features)
    
    # Create node mapping if atom_order provided
    node_mapping = {}
    if atom_order and len(atom_order) == mol.GetNumAtoms():
        for i, atom_id in enumerate(atom_order):
            node_mapping[i] = atom_id
    
    # Lazy import PyTorch - with fallback
    if not _ensure_torch():
        # Fallback: return mock data structure using numpy
        import numpy as np
        logging.warning("PyTorch not available, using numpy fallback for featurization")
        x = np.array(node_features, dtype=np.float32)
        
        # Create mock Data object structure
        class MockData:
            def __init__(self, x, edge_index, edge_attr):
                self.x = x
                self.edge_index = edge_index
                self.edge_attr = edge_attr
                self.num_nodes = len(node_features)
        
        # Create edge indices and features
        edge_indices = []
        edge_features = []
        for bond in mol.GetBonds():
            u = bond.GetBeginAtomIdx()
            v = bond.GetEndAtomIdx()
            edge_indices.append([u, v])
            edge_indices.append([v, u])
            edge_feat = compute_edge_features(bond)
            edge_features.append(edge_feat)
            edge_features.append(edge_feat)
        
        if len(edge_indices) > 0:
            edge_index = np.array(edge_indices, dtype=np.int64).T
            edge_attr = np.array(edge_features, dtype=np.float32)
        else:
            edge_index = np.empty((2, 0), dtype=np.int64)
            edge_attr = np.empty((0, 7), dtype=np.float32)
        
        data = MockData(x, edge_index, edge_attr)
        return data, node_mapping
    
    x = torch.tensor(node_features, dtype=torch.float)
    
    # Edge indices and features
    edge_indices = []
    edge_features = []
    
    for bond in mol.GetBonds():
        u = bond.GetBeginAtomIdx()
        v = bond.GetEndAtomIdx()
        edge_indices.append([u, v])
        edge_indices.append([v, u])  # Undirected graph
        
        edge_feat = compute_edge_features(bond)
        edge_features.append(edge_feat)
        edge_features.append(edge_feat)  # Same for both directions
    
    if len(edge_indices) > 0:
        edge_index = torch.tensor(edge_indices, dtype=torch.long).t().contiguous()
        edge_attr = torch.tensor(edge_features, dtype=torch.float)
    else:
        edge_index = torch.empty((2, 0), dtype=torch.long)
        if mol.GetNumBonds() > 0:
            edge_attr = torch.empty((0, len(compute_edge_features(mol.GetBondWithIdx(0)))), dtype=torch.float)
        else:
            edge_attr = torch.empty((0, 7), dtype=torch.float)  # Default edge feature dim
    
    # Create node mapping if atom_order provided
    node_mapping = {}
    if atom_order and len(atom_order) == mol.GetNumAtoms():
        for i, atom_id in enumerate(atom_order):
            node_mapping[i] = atom_id
    
    data = Data(
        x=x,
        edge_index=edge_index,
        edge_attr=edge_attr,
        num_nodes=mol.GetNumAtoms(),
    )
    
    return data, node_mapping


def featurize_json(payload: Dict) -> Tuple[Data, Dict[int, str]]:
    """
    Featurize molecule from frontend JSON payload.
    
    Expected payload format:
    {
        "atoms": [{"id": "a0", "element": "C", "x": 0, "y": 0, "z": 0, "charge": 0}, ...],
        "bonds": [{"id": "b0", "atoms": ["a0", "a1"], "order": 1}, ...],
        "node_order": ["a0", "a1", ...]  // optional
    }
    """
    atoms = payload.get('atoms', [])
    bonds = payload.get('bonds', [])
    node_order = payload.get('node_order', [])
    
    if not atoms:
        raise ValueError("No atoms in payload")
    
    # Build atom index to atom mapping
    atom_id_to_idx = {atom['id']: idx for idx, atom in enumerate(atoms)}
    
    # Node features
    node_features = []
    for atom in atoms:
        # Get bonds for this atom
        atom_bonds = [b for b in bonds if atom['id'] in b['atoms']]
        bond_order_sum = sum(b.get('order', 1) for b in atom_bonds)
        
        features = [
            float(ELEMENT_TO_ATOMIC_NUMBER.get(atom['element'], 0)),
            float(len(atom_bonds)),  # degree
            float(bond_order_sum),
            float(atom.get('charge', 0)),
            float(0 if bond_order_sum <= 2 else (1 if bond_order_sum <= 3 else 2)),  # hybridization
            float(VALENCE_TABLE.get(atom['element'], 4)),
            float(atom.get('x', 0) / 100.0),
            float(atom.get('y', 0) / 100.0),
            float(atom.get('z', 0) / 100.0),
            0.0,  # aromatic (would need RDKit to compute)
            0.0,  # hybridization enum
            0.0,  # radical electrons
        ]
        node_features.append(features)
    
    # Lazy import PyTorch
    if not _ensure_torch():
        # Fallback: return mock data structure using numpy
        import numpy as np
        logging.warning("PyTorch not available, using numpy fallback for featurization")
        x = np.array(node_features, dtype=np.float32)
        
        # Create mock Data object structure
        class MockData:
            def __init__(self, x, edge_index, edge_attr):
                self.x = x
                self.edge_index = edge_index
                self.edge_attr = edge_attr
                self.num_nodes = len(node_features)
        
        # Create edge indices and features
        edge_indices = []
        edge_features = []
        for bond in mol.GetBonds():
            u = bond.GetBeginAtomIdx()
            v = bond.GetEndAtomIdx()
            edge_indices.append([u, v])
            edge_indices.append([v, u])
            edge_feat = compute_edge_features(bond)
            edge_features.append(edge_feat)
            edge_features.append(edge_feat)
        
        if len(edge_indices) > 0:
            edge_index = np.array(edge_indices, dtype=np.int64).T
            edge_attr = np.array(edge_features, dtype=np.float32)
        else:
            edge_index = np.empty((2, 0), dtype=np.int64)
            edge_attr = np.empty((0, 7), dtype=np.float32)
        
        data = MockData(x, edge_index, edge_attr)
        return data, node_mapping
    
    x = torch.tensor(node_features, dtype=torch.float)
    
    # Edge indices and features
    edge_indices = []
    edge_features = []
    
    for bond in bonds:
        atom1_id = bond['atoms'][0]
        atom2_id = bond['atoms'][1]
        
        u = atom_id_to_idx.get(atom1_id)
        v = atom_id_to_idx.get(atom2_id)
        
        if u is None or v is None:
            continue
        
        edge_indices.append([u, v])
        edge_indices.append([v, u])
        
        bond_order = bond.get('order', 1)
        edge_feat = [
            float(bond_order),
            1.0 if bond_order == 1 else 0.0,
            1.0 if bond_order == 2 else 0.0,
            1.0 if bond_order == 3 else 0.0,
            0.0,  # aromatic
            0.0,  # conjugated
            0.0,  # stereo
        ]
        edge_features.append(edge_feat)
        edge_features.append(edge_feat)
    
    if len(edge_indices) > 0:
        edge_index = torch.tensor(edge_indices, dtype=torch.long).t().contiguous()
        edge_attr = torch.tensor(edge_features, dtype=torch.float)
    else:
        edge_index = torch.empty((2, 0), dtype=torch.long)
        edge_attr = torch.empty((0, 7), dtype=torch.float)
    
    # Node mapping
    node_mapping = {}
    if node_order:
        for i, atom_id in enumerate(node_order):
            if i < len(atoms):
                node_mapping[i] = atom_id
    else:
        # Use atom IDs from payload
        for i, atom in enumerate(atoms):
            node_mapping[i] = atom['id']
    
    data = Data(
        x=x,
        edge_index=edge_index,
        edge_attr=edge_attr,
        num_nodes=len(atoms),
    )
    
    return data, node_mapping


def canonicalize_smiles(smiles: str) -> str:
    """
    Canonicalize SMILES string.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return smiles
        return Chem.MolToSmiles(mol, canonical=True)
    except:
        return smiles

