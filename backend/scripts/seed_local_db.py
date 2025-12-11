import sys
import traceback
from pathlib import Path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from sqlmodel import Session, select
from backend.db import engine
from backend.models.db.molecule import Molecule
from datetime import datetime

ITEMS = [
    {'name': 'Water', 'formula': 'H2O', 'smiles': 'O'},
    {'name': 'Ethanol', 'formula': 'C2H6O', 'smiles': 'CCO'},
    {'name': 'Methane', 'formula': 'CH4', 'smiles': 'C'},
    {'name': 'Benzene', 'formula': 'C6H6', 'smiles': 'c1ccccc1'},
    {'name': 'Acetone', 'formula': 'C3H6O', 'smiles': 'CC(=O)C'},
    {'name': 'Aspirin', 'formula': 'C9H8O4', 'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O'},
    {'name': 'Caffeine', 'formula': 'C8H10N4O2', 'smiles': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'},
    {'name': 'Glucose', 'formula': 'C6H12O6', 'smiles': 'C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O'},
    {'name': 'Ammonia', 'formula': 'NH3', 'smiles': 'N'},
]

def seed_only():
    print("Seeding...")
    try:
        with Session(engine) as session:
            for item in ITEMS:
                statement = select(Molecule).where(Molecule.name == item['name'])
                result = session.exec(statement).first()
                if not result:
                    # Attempt to set owner_id if possible, catch generic error if field unexpected
                    try:
                        mol = Molecule(
                            name=item['name'],
                            formula=item['formula'],
                            smiles=item['smiles'],
                            created_at=datetime.utcnow(),
                            updated_at=datetime.utcnow(),
                            properties={}, 
                            json_graph={},
                            owner_id=1 # Try providing owner_id
                        )
                    except TypeError:
                        # Fallback if owner_id not in __init__
                        mol = Molecule(
                            name=item['name'],
                            formula=item['formula'],
                            smiles=item['smiles'],
                            created_at=datetime.utcnow(),
                            updated_at=datetime.utcnow(),
                            properties={}, 
                            json_graph={}
                        )

                    session.add(mol)
                    print(f"Added {mol.name}")
                else:
                    print(f"Skipped {item['name']}")
            session.commit()
            print("Seeding Complete")
    except Exception:
        traceback.print_exc()

if __name__ == "__main__":
    seed_only()
