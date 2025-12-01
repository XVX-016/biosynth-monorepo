"""
Molecule operations API endpoints

Phase 6: RDKit Backend Integration
Phase 8: 2D Layout Generation

Provides endpoints for:
- SMILES generation
- MolBlock generation
- Hydrogen normalization
- RDKit validation
- 2D coordinate generation
"""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import List, Dict, Any, Optional
import logging
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDepictor

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/molecule", tags=["molecule"])


class MoleculeRequest(BaseModel):
    """Request format for molecule operations."""
    molecule: Dict[str, Any]  # {atoms: [...], bonds: [...]}
    canonicalize: Optional[bool] = True


class ValidateRequest(BaseModel):
    """Request for molecule validation."""
    molecule: Dict[str, Any]
    smiles: Optional[str] = None  # Alternative: provide SMILES directly


class LayoutRequest(BaseModel):
    """Request for 2D layout generation."""
    molecule: Dict[str, Any]
    smiles: Optional[str] = None  # Alternative: provide SMILES directly
    method: Optional[str] = "coordgen"  # "coordgen" or "rdkit"
    spacing: Optional[float] = 1.5  # Bond length in Angstroms


@router.post("/to-smiles")
async def to_smiles(request: MoleculeRequest):
    """
    Convert molecule to SMILES string.
    
    Uses RDKit for accurate SMILES generation.
    """
    try:
        mol = molecule_dict_to_rdkit(request.molecule)
        if not mol:
            raise HTTPException(status_code=400, detail="Invalid molecule structure")
        
        if request.canonicalize:
            smiles = Chem.MolToSmiles(mol)
        else:
            smiles = Chem.MolToSmiles(mol, canonical=False)
        
        return {
            "smiles": smiles,
            "canonical": request.canonicalize,
        }
    except Exception as e:
        logger.error(f"Error generating SMILES: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/to-molblock")
async def to_molblock(request: MoleculeRequest):
    """
    Convert molecule to MolBlock format (V2000).
    """
    try:
        mol = molecule_dict_to_rdkit(request.molecule)
        if not mol:
            raise HTTPException(status_code=400, detail="Invalid molecule structure")
        
        molblock = Chem.MolToMolBlock(mol)
        
        return {
            "molblock": molblock,
        }
    except Exception as e:
        logger.error(f"Error generating MolBlock: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/normalize-hydrogens")
async def normalize_hydrogens(request: MoleculeRequest):
    """
    Normalize hydrogens in molecule (add implicit hydrogens).
    """
    try:
        mol = molecule_dict_to_rdkit(request.molecule)
        if not mol:
            raise HTTPException(status_code=400, detail="Invalid molecule structure")
        
        # Add hydrogens
        mol = Chem.AddHs(mol)
        
        # Convert back to molecule dict
        molecule_dict = rdkit_to_molecule_dict(mol)
        
        return {
            "molecule": molecule_dict,
        }
    except Exception as e:
        logger.error(f"Error normalizing hydrogens: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/generate-2d-layout")
async def generate_2d_layout(request: LayoutRequest):
    """
    Generate 2D coordinates for molecule.
    
    Uses RDKit's coordgen or standard 2D coordinate generation.
    Returns molecule with updated atom positions.
    """
    try:
        # Get molecule from dict or SMILES
        mol = None
        if request.smiles:
            mol = Chem.MolFromSmiles(request.smiles)
            if not mol:
                raise HTTPException(status_code=400, detail="Invalid SMILES string")
        else:
            mol = molecule_dict_to_rdkit(request.molecule)
            if not mol:
                raise HTTPException(status_code=400, detail="Invalid molecule structure")
        
        # Generate 2D coordinates
        if request.method == "coordgen":
            # Use CoordGen (better for complex molecules)
            try:
                rdDepictor.Compute2DCoords(mol)
            except:
                # Fallback to standard method
                AllChem.Compute2DCoords(mol)
        else:
            # Use standard RDKit 2D coordinate generation
            AllChem.Compute2DCoords(mol)
        
        # Optimize layout
        try:
            # Try to improve layout with additional optimization
            rdDepictor.Compute2DCoords(mol, clearConfs=True)
        except:
            pass
        
        # Convert back to molecule dict with 2D coordinates
        molecule_dict = rdkit_to_molecule_dict(mol)
        
        return {
            "molecule": molecule_dict,
            "method": request.method,
        }
    except Exception as e:
        logger.error(f"Error generating 2D layout: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/validate")
async def validate_molecule(request: ValidateRequest):
    """
    Validate molecule with RDKit.
    
    Returns:
    - valid: bool
    - sanitized_smiles: str (if valid)
    - molblock: str (if valid)
    - errors: List[str]
    """
    try:
        errors: List[str] = []
        
        # Try to create molecule from dict or SMILES
        mol = None
        if request.smiles:
            mol = Chem.MolFromSmiles(request.smiles)
            if not mol:
                errors.append(f"Invalid SMILES: {request.smiles}")
        else:
            mol = molecule_dict_to_rdkit(request.molecule)
            if not mol:
                errors.append("Invalid molecule structure")
        
        if not mol:
            return {
                "valid": False,
                "errors": errors,
            }
        
        # Try to sanitize
        try:
            Chem.SanitizeMol(mol)
            sanitized_smiles = Chem.MolToSmiles(mol)
            molblock = Chem.MolToMolBlock(mol)
            
            return {
                "valid": True,
                "sanitized_smiles": sanitized_smiles,
                "molblock": molblock,
                "errors": [],
            }
        except Exception as sanitize_error:
            errors.append(f"Sanitization failed: {str(sanitize_error)}")
            return {
                "valid": False,
                "errors": errors,
            }
            
    except Exception as e:
        logger.error(f"Error validating molecule: {e}", exc_info=True)
        return {
            "valid": False,
            "errors": [str(e)],
        }


def molecule_dict_to_rdkit(molecule_dict: Dict[str, Any]) -> Optional[Chem.Mol]:
    """
    Convert molecule dict to RDKit Mol object.
    
    Expected format:
    {
        "atoms": [{"id": "...", "element": "C", "position": [x, y, z], ...}, ...],
        "bonds": [{"id": "...", "atom1": "...", "atom2": "...", "order": 1}, ...]
    }
    """
    try:
        atoms = molecule_dict.get("atoms", [])
        bonds = molecule_dict.get("bonds", [])
        
        if not atoms:
            return None
        
        # Create RDKit molecule
        mol = Chem.RWMol()
        
        # Map atom IDs to RDKit atom indices
        atom_id_to_idx: Dict[str, int] = {}
        
        # Add atoms
        for atom_data in atoms:
            element = atom_data.get("element", "C")
            # Get atomic number from element symbol
            try:
                atomic_num = Chem.GetPeriodicTable().GetAtomicNumber(element)
            except:
                # Fallback: use common elements
                element_map = {
                    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8,
                    "F": 9, "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15,
                    "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Br": 35, "I": 53,
                }
                atomic_num = element_map.get(element, 6)  # Default to Carbon
            
            rdkit_atom = Chem.Atom(atomic_num)
            
            # Set charge if provided
            if "charge" in atom_data or "formalCharge" in atom_data:
                charge = atom_data.get("formalCharge") or atom_data.get("charge") or 0
                rdkit_atom.SetFormalCharge(int(charge))
            
            idx = mol.AddAtom(rdkit_atom)
            atom_id_to_idx[atom_data["id"]] = idx
        
        # Add bonds
        for bond_data in bonds:
            atom1_id = bond_data.get("atom1")
            atom2_id = bond_data.get("atom2")
            order = bond_data.get("order", 1)
            
            if atom1_id not in atom_id_to_idx or atom2_id not in atom_id_to_idx:
                continue
            
            idx1 = atom_id_to_idx[atom1_id]
            idx2 = atom_id_to_idx[atom2_id]
            
            # Convert bond order
            if order == 1.5:
                # Aromatic bond
                mol.AddBond(idx1, idx2, Chem.BondType.AROMATIC)
            elif order == 1:
                mol.AddBond(idx1, idx2, Chem.BondType.SINGLE)
            elif order == 2:
                mol.AddBond(idx1, idx2, Chem.BondType.DOUBLE)
            elif order == 3:
                mol.AddBond(idx1, idx2, Chem.BondType.TRIPLE)
            else:
                mol.AddBond(idx1, idx2, Chem.BondType.SINGLE)
        
        # Convert to regular Mol
        mol = mol.GetMol()
        
        # Try to sanitize
        try:
            Chem.SanitizeMol(mol)
        except:
            # If sanitization fails, return unsanitized mol
            # (caller can handle errors)
            pass
        
        return mol
        
    except Exception as e:
        logger.error(f"Error converting molecule dict to RDKit: {e}", exc_info=True)
        return None


def rdkit_to_molecule_dict(mol: Chem.Mol) -> Dict[str, Any]:
    """
    Convert RDKit Mol to molecule dict format.
    Extracts 2D or 3D coordinates from conformer.
    """
    atoms = []
    bonds = []
    
    # Get conformer (prefer 2D, fallback to 3D or default)
    conf = None
    if mol.GetNumConformers() > 0:
        # Prefer conformer 0 (usually 2D if generated)
        conf = mol.GetConformer(0)
    
    # Add atoms
    for i, atom in enumerate(mol.GetAtoms()):
        element = atom.GetSymbol()
        
        # Get position from conformer or default to [0, 0, 0]
        if conf:
            pos = conf.GetAtomPosition(i)
            position = [pos.x, pos.y, pos.z if pos.z else 0.0]
        else:
            position = [0.0, 0.0, 0.0]
        
        atoms.append({
            "id": f"atom_{i}",
            "element": element,
            "position": position,
            "charge": atom.GetFormalCharge(),
            "formalCharge": atom.GetFormalCharge(),
        })
    
    # Add bonds
    for i, bond in enumerate(mol.GetBonds()):
        order = bond.GetBondTypeAsDouble()
        bonds.append({
            "id": f"bond_{i}",
            "atom1": f"atom_{bond.GetBeginAtomIdx()}",
            "atom2": f"atom_{bond.GetEndAtomIdx()}",
            "order": int(order) if order == int(order) else 1.5,
        })
    
    return {
        "atoms": atoms,
        "bonds": bonds,
    }
