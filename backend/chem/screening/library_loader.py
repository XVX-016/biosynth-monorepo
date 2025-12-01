"""
Library Loader - Load molecular libraries and build search indices

Supports loading from various formats (SDF, SMILES, CSV) and indexing.
"""

import csv
import json
from typing import List, Dict, Optional, Iterator
from pathlib import Path
import logging

from .fingerprint_index import FingerprintIndex, get_fingerprint_index, compute_ecfp_fingerprint

logger = logging.getLogger(__name__)


class LibraryLoader:
    """
    Load molecular libraries from files and build search indices.
    """
    
    def __init__(self, index: Optional[FingerprintIndex] = None):
        """
        Args:
            index: FingerprintIndex instance (default: global instance)
        """
        self.index = index or get_fingerprint_index()
    
    def load_from_smiles_file(
        self,
        filepath: str,
        delimiter: str = '\t',
        smiles_column: int = 0,
        name_column: Optional[int] = None,
        id_column: Optional[int] = None,
        skip_header: bool = True
    ) -> int:
        """
        Load molecules from SMILES file (TSV/CSV format).
        
        Args:
            filepath: Path to file
            delimiter: Column delimiter
            smiles_column: Column index for SMILES
            name_column: Optional column index for name
            id_column: Optional column index for ID
            skip_header: Skip first line if header
        
        Returns:
            Number of molecules loaded
        """
        path = Path(filepath)
        if not path.exists():
            raise FileNotFoundError(f"File not found: {filepath}")
        
        count = 0
        
        with open(path, 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter=delimiter)
            
            if skip_header:
                next(reader, None)
            
            for row_num, row in enumerate(reader, start=2 if skip_header else 1):
                if len(row) <= smiles_column:
                    logger.warning(f"Row {row_num}: insufficient columns")
                    continue
                
                smiles = row[smiles_column].strip()
                if not smiles:
                    continue
                
                # Generate molecule ID
                molecule_id = None
                if id_column is not None and id_column < len(row):
                    molecule_id = row[id_column].strip()
                
                if not molecule_id:
                    molecule_id = f"mol_{count + 1}"
                
                # Get name
                name = None
                if name_column is not None and name_column < len(row):
                    name = row[name_column].strip()
                
                # Compute fingerprint
                fingerprint = compute_ecfp_fingerprint(smiles)
                if not fingerprint:
                    logger.warning(f"Row {row_num}: failed to compute fingerprint for {smiles}")
                    continue
                
                # Add to index
                metadata = {
                    'smiles': smiles,
                    'name': name or molecule_id,
                    'source_file': str(path),
                    'row_number': row_num,
                }
                
                self.index.add_molecule(molecule_id, fingerprint, metadata)
                count += 1
        
        logger.info(f"Loaded {count} molecules from {filepath}")
        return count
    
    def load_from_json(
        self,
        filepath: str,
        smiles_key: str = 'smiles',
        id_key: str = 'id',
        name_key: str = 'name'
    ) -> int:
        """
        Load molecules from JSON file.
        
        Args:
            filepath: Path to JSON file
            smiles_key: Key for SMILES in each object
            id_key: Key for molecule ID
            name_key: Key for molecule name
        
        Returns:
            Number of molecules loaded
        """
        path = Path(filepath)
        if not path.exists():
            raise FileNotFoundError(f"File not found: {filepath}")
        
        with open(path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        if not isinstance(data, list):
            raise ValueError("JSON file must contain a list of molecules")
        
        count = 0
        
        for idx, mol_data in enumerate(data):
            if not isinstance(mol_data, dict):
                continue
            
            smiles = mol_data.get(smiles_key)
            if not smiles:
                continue
            
            molecule_id = mol_data.get(id_key) or f"mol_{idx + 1}"
            name = mol_data.get(name_key) or molecule_id
            
            # Compute fingerprint
            fingerprint = compute_ecfp_fingerprint(smiles)
            if not fingerprint:
                logger.warning(f"Failed to compute fingerprint for {smiles}")
                continue
            
            # Add to index
            metadata = {
                'smiles': smiles,
                'name': name,
                'source_file': str(path),
                **{k: v for k, v in mol_data.items() if k not in [smiles_key, id_key, name_key]},
            }
            
            self.index.add_molecule(molecule_id, fingerprint, metadata)
            count += 1
        
        logger.info(f"Loaded {count} molecules from {filepath}")
        return count
    
    def load_from_sdf(
        self,
        filepath: str,
        use_rdkit: bool = True
    ) -> int:
        """
        Load molecules from SDF file.
        
        Args:
            filepath: Path to SDF file
            use_rdkit: Use RDKit for parsing (if available)
        
        Returns:
            Number of molecules loaded
        """
        # TODO: Implement SDF parsing
        # For now, return 0 and log warning
        logger.warning("SDF loading not yet implemented")
        return 0
    
    def load_from_molecule_dicts(
        self,
        molecules: List[Dict],
        id_key: str = 'id',
        smiles_key: str = 'smiles',
        name_key: str = 'name'
    ) -> int:
        """
        Load molecules from list of dicts (e.g., from database).
        
        Args:
            molecules: List of molecule dicts
            id_key: Key for molecule ID
            smiles_key: Key for SMILES
            name_key: Key for name
        
        Returns:
            Number of molecules loaded
        """
        count = 0
        
        for mol_data in molecules:
            if not isinstance(mol_data, dict):
                continue
            
            smiles = mol_data.get(smiles_key)
            if not smiles:
                continue
            
            molecule_id = mol_data.get(id_key) or f"mol_{count + 1}"
            name = mol_data.get(name_key) or molecule_id
            
            # Compute fingerprint
            fingerprint = compute_ecfp_fingerprint(smiles)
            if not fingerprint:
                continue
            
            # Add to index
            metadata = {
                'smiles': smiles,
                'name': name,
                **{k: v for k, v in mol_data.items() if k not in [id_key, smiles_key, name_key]},
            }
            
            self.index.add_molecule(molecule_id, fingerprint, metadata)
            count += 1
        
        return count

