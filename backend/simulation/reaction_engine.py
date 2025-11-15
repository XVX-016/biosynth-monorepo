"""
Reaction engine for molecular transformations
"""
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
from rdkit import Chem
from rdkit.Chem import AllChem


@dataclass
class ReactionProduct:
    """Represents a reaction product"""
    smiles: str
    metadata: Dict[str, any]


@dataclass
class ReactionResult:
    """Result of a reaction operation"""
    success: bool
    products: List[ReactionProduct]
    error: Optional[str] = None
    metadata: Optional[Dict[str, any]] = None


class ReactionEngine:
    """Engine for performing molecular reactions"""
    
    def __init__(self):
        self.reaction_history: List[Dict] = []
    
    def bond_break(
        self,
        smiles: str,
        atom1_idx: int,
        atom2_idx: int,
        validate: bool = True
    ) -> ReactionResult:
        """
        Break a bond between two atoms
        
        Args:
            smiles: SMILES string of the molecule
            atom1_idx: Index of first atom (0-based)
            atom2_idx: Index of second atom (0-based)
            validate: Whether to validate the reaction
        
        Returns:
            ReactionResult with products
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return ReactionResult(
                    success=False,
                    products=[],
                    error="Invalid SMILES string"
                )
            
            # Convert to editable molecule
            mol_editable = Chem.EditableMol(mol)
            
            # Find bond between atoms
            bond = mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)
            if bond is None:
                return ReactionResult(
                    success=False,
                    products=[],
                    error=f"No bond found between atoms {atom1_idx} and {atom2_idx}"
                )
            
            # Remove bond
            mol_editable.RemoveBond(atom1_idx, atom2_idx)
            new_mol = mol_editable.GetMol()
            
            # Check if molecule split into fragments
            fragments = Chem.GetMolFrags(new_mol, asMols=True)
            
            products = []
            for frag in fragments:
                frag_smiles = Chem.MolToSmiles(frag)
                products.append(ReactionProduct(
                    smiles=frag_smiles,
                    metadata={
                        "reaction_type": "bond_break",
                        "original_smiles": smiles,
                        "broken_atoms": [atom1_idx, atom2_idx]
                    }
                ))
            
            result = ReactionResult(
                success=True,
                products=products,
                metadata={
                    "num_fragments": len(fragments),
                    "original_atoms": mol.GetNumAtoms()
                }
            )
            
            if validate:
                validation = self.validate_reaction(smiles, [p.smiles for p in products])
                if not validation["valid"]:
                    result.error = validation.get("error")
            
            return result
            
        except Exception as e:
            return ReactionResult(
                success=False,
                products=[],
                error=f"Error breaking bond: {str(e)}"
            )
    
    def bond_form(
        self,
        smiles1: str,
        smiles2: str,
        atom1_idx: Optional[int] = None,
        atom2_idx: Optional[int] = None,
        validate: bool = True
    ) -> ReactionResult:
        """
        Form a bond between two molecules
        
        Args:
            smiles1: SMILES string of first molecule
            smiles2: SMILES string of second molecule
            atom1_idx: Index of atom in first molecule (None for auto-select)
            atom2_idx: Index of atom in second molecule (None for auto-select)
            validate: Whether to validate the reaction
        
        Returns:
            ReactionResult with combined product
        """
        try:
            mol1 = Chem.MolFromSmiles(smiles1)
            mol2 = Chem.MolFromSmiles(smiles2)
            
            if mol1 is None or mol2 is None:
                return ReactionResult(
                    success=False,
                    products=[],
                    error="Invalid SMILES string(s)"
                )
            
            # Combine molecules
            combined = Chem.CombineMols(mol1, mol2)
            
            # If atom indices specified, form bond
            if atom1_idx is not None and atom2_idx is not None:
                # Adjust atom2_idx for combined molecule
                atom2_idx_adjusted = atom2_idx + mol1.GetNumAtoms()
                
                # Create editable molecule
                mol_editable = Chem.EditableMol(combined)
                mol_editable.AddBond(atom1_idx, atom2_idx_adjusted, Chem.BondType.SINGLE)
                combined = mol_editable.GetMol()
            
            # Generate SMILES
            combined_smiles = Chem.MolToSmiles(combined)
            
            products = [ReactionProduct(
                smiles=combined_smiles,
                metadata={
                    "reaction_type": "bond_form",
                    "reactant1": smiles1,
                    "reactant2": smiles2,
                    "bond_atoms": [atom1_idx, atom2_idx] if atom1_idx is not None else None
                }
            )]
            
            result = ReactionResult(
                success=True,
                products=products,
                metadata={
                    "num_atoms": combined.GetNumAtoms(),
                    "num_bonds": combined.GetNumBonds()
                }
            )
            
            if validate:
                validation = self.validate_reaction(
                    f"{smiles1}.{smiles2}",
                    [combined_smiles]
                )
                if not validation["valid"]:
                    result.error = validation.get("error")
            
            return result
            
        except Exception as e:
            return ReactionResult(
                success=False,
                products=[],
                error=f"Error forming bond: {str(e)}"
            )
    
    def ionize(
        self,
        smiles: str,
        atom_idx: int,
        charge: int = 1,
        validate: bool = True
    ) -> ReactionResult:
        """
        Ionize a molecule (add/remove charge)
        
        Args:
            smiles: SMILES string
            atom_idx: Index of atom to ionize
            charge: Charge to add (positive for cation, negative for anion)
            validate: Whether to validate the reaction
        
        Returns:
            ReactionResult with ionized product
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return ReactionResult(
                    success=False,
                    products=[],
                    error="Invalid SMILES string"
                )
            
            # Create copy with formal charge
            mol_copy = Chem.Mol(mol)
            atom = mol_copy.GetAtomWithIdx(atom_idx)
            atom.SetFormalCharge(atom.GetFormalCharge() + charge)
            
            # Regenerate SMILES
            ionized_smiles = Chem.MolToSmiles(mol_copy)
            
            products = [ReactionProduct(
                smiles=ionized_smiles,
                metadata={
                    "reaction_type": "ionize",
                    "original_smiles": smiles,
                    "ionized_atom": atom_idx,
                    "charge_added": charge
                }
            )]
            
            result = ReactionResult(
                success=True,
                products=products,
                metadata={
                    "total_charge": sum(atom.GetFormalCharge() for atom in mol_copy.GetAtoms())
                }
            )
            
            if validate:
                validation = self.validate_reaction(smiles, [ionized_smiles])
                if not validation["valid"]:
                    result.error = validation.get("error")
            
            return result
            
        except Exception as e:
            return ReactionResult(
                success=False,
                products=[],
                error=f"Error ionizing molecule: {str(e)}"
            )
    
    def recombine(
        self,
        fragments: List[str],
        validate: bool = True
    ) -> ReactionResult:
        """
        Recombine molecular fragments
        
        Args:
            fragments: List of SMILES strings for fragments
            validate: Whether to validate the reaction
        
        Returns:
            ReactionResult with recombined product
        """
        try:
            if len(fragments) < 2:
                return ReactionResult(
                    success=False,
                    products=[],
                    error="Need at least 2 fragments to recombine"
                )
            
            # Combine all fragments
            mols = []
            for frag in fragments:
                mol = Chem.MolFromSmiles(frag)
                if mol is None:
                    return ReactionResult(
                        success=False,
                        products=[],
                        error=f"Invalid SMILES in fragment: {frag}"
                    )
                mols.append(mol)
            
            # Combine molecules
            combined = mols[0]
            for mol in mols[1:]:
                combined = Chem.CombineMols(combined, mol)
            
            # Generate SMILES
            recombined_smiles = Chem.MolToSmiles(combined)
            
            products = [ReactionProduct(
                smiles=recombined_smiles,
                metadata={
                    "reaction_type": "recombine",
                    "fragments": fragments,
                    "num_fragments": len(fragments)
                }
            )]
            
            result = ReactionResult(
                success=True,
                products=products,
                metadata={
                    "num_atoms": combined.GetNumAtoms(),
                    "num_bonds": combined.GetNumBonds()
                }
            )
            
            if validate:
                validation = self.validate_reaction(
                    ".".join(fragments),
                    [recombined_smiles]
                )
                if not validation["valid"]:
                    result.error = validation.get("error")
            
            return result
            
        except Exception as e:
            return ReactionResult(
                success=False,
                products=[],
                error=f"Error recombining fragments: {str(e)}"
            )
    
    def validate_reaction(
        self,
        reactants: str,
        products: List[str]
    ) -> Dict[str, any]:
        """
        Validate a reaction
        
        Args:
            reactants: SMILES string of reactants (can be multiple separated by '.')
            products: List of product SMILES strings
        
        Returns:
            Dictionary with validation results
        """
        try:
            # Parse reactants
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants.split('.')]
            if any(m is None for m in reactant_mols):
                return {
                    "valid": False,
                    "error": "Invalid reactant SMILES"
                }
            
            # Parse products
            product_mols = [Chem.MolFromSmiles(p) for p in products]
            if any(m is None for m in product_mols):
                return {
                    "valid": False,
                    "error": "Invalid product SMILES"
                }
            
            # Count atoms (basic validation)
            reactant_atoms = sum(m.GetNumAtoms() for m in reactant_mols)
            product_atoms = sum(m.GetNumAtoms() for m in product_mols)
            
            # Allow some flexibility (atoms can be added/removed in some reactions)
            atom_balance = abs(reactant_atoms - product_atoms)
            
            # Basic validation: check if molecules are valid
            all_valid = all(
                mol is not None and mol.GetNumAtoms() > 0
                for mol in reactant_mols + product_mols
            )
            
            return {
                "valid": all_valid and atom_balance <= 10,  # Allow some flexibility
                "reactant_atoms": reactant_atoms,
                "product_atoms": product_atoms,
                "atom_balance": atom_balance,
                "error": None if all_valid else "Invalid molecule structure"
            }
            
        except Exception as e:
            return {
                "valid": False,
                "error": f"Validation error: {str(e)}"
            }


# Singleton instance
reaction_engine = ReactionEngine()

