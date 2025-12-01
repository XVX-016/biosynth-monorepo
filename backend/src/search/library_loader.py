"""
Library Loader - Lazy streaming loader for .smi files

Loads molecules from /data/libraries/*.smi files with lazy evaluation.
"""

from typing import Iterator, Dict, Optional, List
from pathlib import Path
import logging

from .fingerprint_index import FingerprintIndex
from .rdkit_index import compute_fingerprint, validate_smiles

logger = logging.getLogger(__name__)


class LibraryLoader:
    """
    Lazy loader for molecular libraries from .smi files.
    
    Supports streaming loading to avoid loading entire libraries into memory.
    """
    
    def __init__(self, index: Optional[FingerprintIndex] = None):
        """
        Args:
            index: FingerprintIndex instance to populate
        """
        self.index = index or FingerprintIndex()
        self.data_dir = Path("data/libraries")
    
    def stream_smi_file(self, filepath: str) -> Iterator[Dict]:
        """
        Stream molecules from .smi file (lazy loading).
        
        .smi format: SMILES [optional_name] [optional_id]
        
        Args:
            filepath: Path to .smi file
        
        Yields:
            Dict with 'smiles', 'name', 'id' keys
        """
        path = Path(filepath)
        if not path.exists():
            logger.error(f"File not found: {filepath}")
            return
        
        with open(path, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, start=1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                # Parse line: SMILES [name] [id]
                parts = line.split()
                if not parts:
                    continue
                
                smiles = parts[0]
                name = parts[1] if len(parts) > 1 else None
                mol_id = parts[2] if len(parts) > 2 else None
                
                # Generate ID if not provided
                if not mol_id:
                    mol_id = f"{path.stem}_{line_num}"
                
                yield {
                    'smiles': smiles,
                    'name': name or mol_id,
                    'id': mol_id,
                    'source_file': str(path),
                    'line_number': line_num,
                }
    
    def load_from_file(
        self,
        filepath: str,
        max_molecules: Optional[int] = None,
        validate: bool = True
    ) -> int:
        """
        Load molecules from .smi file into index.
        
        Args:
            filepath: Path to .smi file
            max_molecules: Maximum molecules to load (None = all)
            validate: Validate SMILES before adding
        
        Returns:
            Number of molecules loaded
        """
        count = 0
        
        for mol_data in self.stream_smi_file(filepath):
            if max_molecules and count >= max_molecules:
                break
            
            smiles = mol_data['smiles']
            
            # Validate if requested
            if validate and not validate_smiles(smiles):
                logger.warning(f"Invalid SMILES at line {mol_data.get('line_number')}: {smiles}")
                continue
            
            # Compute fingerprint
            fingerprint = compute_fingerprint(smiles)
            if not fingerprint:
                logger.warning(f"Failed to compute fingerprint for: {smiles}")
                continue
            
            # Add to index
            molecule_id = mol_data['id']
            metadata = {
                'smiles': smiles,
                'name': mol_data.get('name', molecule_id),
                'source_file': mol_data.get('source_file'),
                'line_number': mol_data.get('line_number'),
            }
            
            self.index.add_molecule(molecule_id, fingerprint, metadata)
            count += 1
        
        logger.info(f"Loaded {count} molecules from {filepath}")
        return count
    
    def load_from_directory(
        self,
        directory: Optional[str] = None,
        pattern: str = "*.smi",
        max_molecules_per_file: Optional[int] = None
    ) -> int:
        """
        Load all .smi files from directory.
        
        Args:
            directory: Directory path (default: data/libraries)
            pattern: File pattern (default: *.smi)
            max_molecules_per_file: Max molecules per file (None = all)
        
        Returns:
            Total number of molecules loaded
        """
        if directory is None:
            directory = str(self.data_dir)
        
        dir_path = Path(directory)
        if not dir_path.exists():
            logger.warning(f"Directory not found: {directory}")
            return 0
        
        total_count = 0
        
        for smi_file in dir_path.glob(pattern):
            logger.info(f"Loading library: {smi_file}")
            count = self.load_from_file(
                str(smi_file),
                max_molecules=max_molecules_per_file
            )
            total_count += count
        
        logger.info(f"Loaded {total_count} total molecules from {directory}")
        return total_count
    
    def get_index(self) -> FingerprintIndex:
        """Get the fingerprint index."""
        return self.index

