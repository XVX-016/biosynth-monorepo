"""
Library Loader - Lazy streaming loader for .smi files

Loads molecules from .smi files with validation.
"""

from typing import Iterator, Dict, Optional
from pathlib import Path
import logging

from .fingerprint_index import FingerprintIndex
from .rdkit_index import compute_ecfp
from backend.chem.utils.validators import validate_smiles

logger = logging.getLogger(__name__)


class LibraryLoader:
    """
    Lazy loader for molecular libraries from .smi files.
    """
    
    def __init__(self, index: Optional[FingerprintIndex] = None):
        """
        Args:
            index: FingerprintIndex instance to populate
        """
        from .fingerprint_index import FingerprintIndex
        self.index = index or FingerprintIndex()
    
    def stream_smi_file(self, filepath: str) -> Iterator[Dict]:
        """
        Stream molecules from .smi file (lazy loading).
        
        .smi format: SMILES [optional_name]
        
        Args:
            filepath: Path to .smi file
        
        Yields:
            Dict with 'smiles', 'name' keys
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
                
                # Parse line: SMILES [name]
                parts = line.split()
                if not parts:
                    continue
                
                smiles = parts[0]
                name = parts[1] if len(parts) > 1 else None
                
                yield {
                    'smiles': smiles,
                    'name': name or smiles,
                    'line_number': line_num,
                }
    
    def load_from_directory(
        self,
        directory: str = "data/libraries",
        pattern: str = "*.smi"
    ) -> int:
        """
        Load all .smi files from directory.
        
        Args:
            directory: Directory path
            pattern: File pattern (default: *.smi)
        
        Returns:
            Number of molecules loaded
        """
        dir_path = Path(directory)
        if not dir_path.exists():
            logger.warning(f"Directory not found: {directory}")
            return 0
        
        total_count = 0
        
        for smi_file in dir_path.glob(pattern):
            logger.info(f"Loading library: {smi_file}")
            count = self._load_file(str(smi_file))
            total_count += count
        
        logger.info(f"Loaded {total_count} total molecules from {directory}")
        return total_count
    
    def _load_file(self, filepath: str) -> int:
        """
        Load molecules from single .smi file.
        
        Args:
            filepath: Path to .smi file
        
        Returns:
            Number of molecules loaded
        """
        count = 0
        
        for mol_data in self.stream_smi_file(filepath):
            smiles = mol_data['smiles']
            
            # Validate SMILES
            if not validate_smiles(smiles):
                logger.warning(f"Invalid SMILES at line {mol_data.get('line_number')}: {smiles}")
                continue
            
            # Compute fingerprint
            fp = compute_ecfp(smiles)
            if not fp:
                logger.warning(f"Failed to compute fingerprint for: {smiles}")
                continue
            
            # Add to index with metadata
            metadata = {
                'name': mol_data.get('name', smiles),
                'smiles': smiles,
            }
            
            self.index.add(smiles, fp, metadata)
            count += 1
        
        return count
    
    def get_index(self) -> FingerprintIndex:
        """Get the fingerprint index."""
        return self.index
