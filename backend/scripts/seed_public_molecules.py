#!/usr/bin/env python3
"""
Seed public_molecules table with 20 common molecules
Generates molfiles and thumbnails for each molecule
"""
import sys
from pathlib import Path
import json

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from backend.services.thumbnail_service import ThumbnailService

# 20 common molecules (reduced from 50 for initial seeding)
SAMPLE_MOLECULES = [
    {'name': 'Water', 'formula': 'H2O', 'smiles': 'O'},
    {'name': 'Ethanol', 'formula': 'C2H6O', 'smiles': 'CCO'},
    {'name': 'Methane', 'formula': 'CH4', 'smiles': 'C'},
    {'name': 'Carbon Dioxide', 'formula': 'CO2', 'smiles': 'O=C=O'},
    {'name': 'Oxygen', 'formula': 'O2', 'smiles': 'O=O'},
    {'name': 'Nitrogen', 'formula': 'N2', 'smiles': 'N#N'},
    {'name': 'Glucose', 'formula': 'C6H12O6', 'smiles': 'C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O'},
    {'name': 'Sucrose', 'formula': 'C12H22O11', 'smiles': 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)OC([C@@H]2[C@H]([C@@H]([C@H](C(O2)O)O)O)O)O)O)O)O)O'},
    {'name': 'Aspirin', 'formula': 'C9H8O4', 'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O'},
    {'name': 'Paracetamol', 'formula': 'C8H9NO2', 'smiles': 'CC(=O)NC1=CC=C(O)C=C1'},
    {'name': 'Caffeine', 'formula': 'C8H10N4O2', 'smiles': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'},
    {'name': 'Benzene', 'formula': 'C6H6', 'smiles': 'c1ccccc1'},
    {'name': 'Toluene', 'formula': 'C7H8', 'smiles': 'Cc1ccccc1'},
    {'name': 'Acetone', 'formula': 'C3H6O', 'smiles': 'CC(=O)C'},
    {'name': 'Acetic Acid', 'formula': 'C2H4O2', 'smiles': 'CC(=O)O'},
    {'name': 'Formic Acid', 'formula': 'CH2O2', 'smiles': 'C(=O)O'},
    {'name': 'Phenol', 'formula': 'C6H6O', 'smiles': 'c1ccc(cc1)O'},
    {'name': 'Aniline', 'formula': 'C6H7N', 'smiles': 'c1ccc(cc1)N'},
    {'name': 'Urea', 'formula': 'CH4N2O', 'smiles': 'C(=O)(N)N'},
    {'name': 'Citric Acid', 'formula': 'C6H8O7', 'smiles': 'C(C(=O)O)C(CC(=O)O)(C(=O)O)O'},
    {'name': 'Lactic Acid', 'formula': 'C3H6O3', 'smiles': 'CC(C(=O)O)O'},
    {'name': 'Methanol', 'formula': 'CH4O', 'smiles': 'CO'},
    {'name': 'Ammonia', 'formula': 'NH3', 'smiles': 'N'},
    # Total: 20 molecules for initial seeding
]

def generate_molecule_data(mol_data: dict) -> dict:
    """Generate complete molecule data with molfile and thumbnail"""
    print(f"Processing {mol_data['name']}...")
    
    try:
        # Generate 3D molfile
        molfile_3d = ThumbnailService.generate_3d_molfile(smiles=mol_data['smiles'])
        
        # Generate thumbnail
        thumbnail_b64 = ThumbnailService.generate_2d_thumbnail_base64(
            smiles=mol_data['smiles'],
            size=(600, 600)
        )
        
        return {
            'name': mol_data['name'],
            'formula': mol_data['formula'],
            'smiles': mol_data['smiles'],
            'molfile': molfile_3d,
            'thumbnail_b64': thumbnail_b64,
            'metadata': json.dumps({'description': f'Common molecule: {mol_data["name"]}'}),
        }
    except Exception as e:
        print(f"  ⚠️  Error processing {mol_data['name']}: {e}")
        return {
            'name': mol_data['name'],
            'formula': mol_data['formula'],
            'smiles': mol_data['smiles'],
            'molfile': None,
            'thumbnail_b64': None,
            'metadata': json.dumps({'description': f'Common molecule: {mol_data["name"]}'}),
        }

def generate_sql_inserts():
    """Generate SQL INSERT statements for Supabase"""
    print(f"Generating data for {len(SAMPLE_MOLECULES)} molecules...\n")
    
    molecules = []
    for mol_data in SAMPLE_MOLECULES:
        result = generate_molecule_data(mol_data)
        molecules.append(result)
        print(f"  ✓ {result['name']} - {result['formula']}")
    
    print("\n" + "="*60)
    print("SQL INSERT statements for public_molecules:")
    print("="*60 + "\n")
    
    for mol in molecules:
        name = mol['name'].replace("'", "''")
        formula = mol['formula'].replace("'", "''") if mol['formula'] else 'NULL'
        smiles = mol['smiles'].replace("'", "''") if mol['smiles'] else 'NULL'
        molfile = mol['molfile'].replace("'", "''") if mol['molfile'] and len(mol['molfile']) < 10000 else 'NULL'
        thumbnail = mol['thumbnail_b64'].replace("'", "''") if mol['thumbnail_b64'] and len(mol['thumbnail_b64']) < 50000 else 'NULL'
        metadata = mol['metadata'].replace("'", "''") if mol['metadata'] else 'NULL'
        
        sql = f"""INSERT INTO public_molecules (name, formula, smiles, molfile, thumbnail_b64, metadata, created_at, updated_at)
VALUES ('{name}', '{formula}', '{smiles}', 
        {'NULL' if molfile == 'NULL' else f"'{molfile[:10000]}'"},
        {'NULL' if thumbnail == 'NULL' else f"'{thumbnail[:50000]}'"},
        '{metadata}', NOW(), NOW());"""
        
        print(sql)
        print()

if __name__ == "__main__":
    generate_sql_inserts()

