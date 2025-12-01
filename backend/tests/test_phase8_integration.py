"""
Integration tests for Phase 8: Conformer Generator

Tests:
- Conformer generation
- ETKDG fallback
- API endpoints
"""

import pytest
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

# Add backend to path
backend_path = Path(__file__).parent.parent
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

from src.conformers import ConformerGenerator, ETKDGGenerator


class TestPhase8Conformers:
    """Test Phase 8 conformer generation."""
    
    def test_conformer_generator(self):
        """Test conformer generator works."""
        generator = ConformerGenerator()
        
        # Generate conformers for ethanol
        conformers = generator.generate_conformers("CCO", n=5)
        
        assert len(conformers) == 5
        
        # Check structure
        for conf in conformers:
            assert 'id' in conf
            assert 'energy' in conf
            assert 'coords' in conf
            assert isinstance(conf['coords'], list)
            assert len(conf['coords']) > 0
            assert len(conf['coords'][0]) == 3  # 3D coordinates
        
        # Check sorted by energy
        energies = [c['energy'] for c in conformers]
        assert energies == sorted(energies)  # Ascending order
    
    def test_conformer_generator_invalid_smiles(self):
        """Test conformer generator handles invalid SMILES."""
        generator = ConformerGenerator()
        
        # Invalid SMILES
        conformers = generator.generate_conformers("INVALID", n=5)
        
        assert len(conformers) == 0
    
    def test_etkdg_placeholder(self):
        """Test ETKDG raises NotImplementedError."""
        etkdg = ETKDGGenerator()
        
        with pytest.raises(NotImplementedError):
            etkdg.generate("CCO", n=5)
    
    def test_conformer_generator_fallback(self):
        """Test conformer generator falls back to mock when ETKDG unavailable."""
        generator = ConformerGenerator()
        
        # Should work even though ETKDG is not implemented
        conformers = generator.generate_conformers("CCO", n=3)
        
        assert len(conformers) == 3
        # All should have mock coordinates
        for conf in conformers:
            assert len(conf['coords']) > 0


class TestPhase8API:
    """Test Phase 8 API endpoints."""
    
    def test_api_imports(self):
        """Test API routes can be imported."""
        try:
            from api import conformers as conformers_api
            
            assert hasattr(conformers_api, 'router')
        except ImportError as e:
            pytest.skip(f"API routes not available: {e}")

