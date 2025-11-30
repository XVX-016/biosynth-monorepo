"""
Validation endpoint - Uses RDKit to validate molecular structures
"""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import List, Optional
import sys
from pathlib import Path

# Add backend to path
ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

router = APIRouter(prefix="/api/validate", tags=["validation"])


class ValidateRequest(BaseModel):
    smiles: str


class ValidationIssue(BaseModel):
    type: str  # 'error' | 'warning' | 'info'
    message: str
    atomId: Optional[str] = None
    bondId: Optional[str] = None


class ValidateResponse(BaseModel):
    valid: bool
    issues: List[ValidationIssue]


@router.post("", response_model=ValidateResponse)
async def validate_molecule(req: ValidateRequest):
    """
    Validate a molecule structure using RDKit
    
    Returns validation results with any issues found
    """
    if not RDKIT_AVAILABLE:
        raise HTTPException(
            status_code=503,
            detail="RDKit is not available. Please install RDKit to use validation."
        )

    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(req.smiles)
        
        if mol is None:
            return ValidateResponse(
                valid=False,
                issues=[
                    ValidationIssue(
                        type="error",
                        message=f"Invalid SMILES string: {req.smiles}"
                    )
                ]
            )

        issues: List[ValidationIssue] = []

        # Sanitize molecule (this will catch many issues)
        try:
            Chem.SanitizeMol(mol)
        except Exception as e:
            issues.append(
                ValidationIssue(
                    type="error",
                    message=f"Sanitization failed: {str(e)}"
                )
            )

        # Check for valence issues
        for atom in mol.GetAtoms():
            try:
                # RDKit will have already checked valence during sanitization
                # But we can add custom checks here
                pass
            except Exception as e:
                issues.append(
                    ValidationIssue(
                        type="warning",
                        message=f"Atom {atom.GetIdx()}: {str(e)}",
                        atomId=str(atom.GetIdx())
                    )
                )

        # Check for disconnected fragments
        frags = Chem.GetMolFrags(mol, asMols=True)
        if len(frags) > 1:
            issues.append(
                ValidationIssue(
                    type="warning",
                    message=f"Molecule contains {len(frags)} disconnected fragments"
                )
            )

        # Check for unusual charges
        for atom in mol.GetAtoms():
            charge = atom.GetFormalCharge()
            if charge != 0:
                issues.append(
                    ValidationIssue(
                        type="info",
                        message=f"Atom {atom.GetSymbol()} has formal charge {charge}",
                        atomId=str(atom.GetIdx())
                    )
                )

        return ValidateResponse(
            valid=len([i for i in issues if i.type == "error"]) == 0,
            issues=issues
        )

    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Validation error: {str(e)}"
        )

