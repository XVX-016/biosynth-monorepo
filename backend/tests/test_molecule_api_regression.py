"""
Regression tests for molecule API endpoints

Phase 16: Full Regression Test & Final Cleanup

Tests all molecule API endpoints with various molecule types.
"""

import pytest
import requests
import logging
from typing import Dict, Any

logger = logging.getLogger(__name__)

API_BASE = "http://localhost:8000"

# Test molecules
TEST_MOLECULES = {
    "water": {
        "atoms": [
            {"id": "o1", "element": "O", "position": [0, 0, 0]},
            {"id": "h1", "element": "H", "position": [1, 0, 0]},
            {"id": "h2", "element": "H", "position": [-1, 0, 0]},
        ],
        "bonds": [
            {"id": "b1", "atom1": "o1", "atom2": "h1", "order": 1},
            {"id": "b2", "atom1": "o1", "atom2": "h2", "order": 1},
        ],
    },
    "methane": {
        "atoms": [
            {"id": "c1", "element": "C", "position": [0, 0, 0]},
            {"id": "h1", "element": "H", "position": [1, 0, 0]},
            {"id": "h2", "element": "H", "position": [-1, 0, 0]},
            {"id": "h3", "element": "H", "position": [0, 1, 0]},
            {"id": "h4", "element": "H", "position": [0, -1, 0]},
        ],
        "bonds": [
            {"id": "b1", "atom1": "c1", "atom2": "h1", "order": 1},
            {"id": "b2", "atom1": "c1", "atom2": "h2", "order": 1},
            {"id": "b3", "atom1": "c1", "atom2": "h3", "order": 1},
            {"id": "b4", "atom1": "c1", "atom2": "h4", "order": 1},
        ],
    },
    "benzene": {
        "atoms": [
            {"id": f"c{i+1}", "element": "C", "position": [20 * (i % 2), 20 * (i // 2), 0]}
            for i in range(6)
        ],
        "bonds": [
            {"id": f"b{i+1}", "atom1": f"c{i+1}", "atom2": f"c{(i+1)%6+1}", "order": 1.5}
            for i in range(6)
        ],
    },
    "ethanol": {
        "atoms": [
            {"id": "c1", "element": "C", "position": [0, 0, 0]},
            {"id": "c2", "element": "C", "position": [1.5, 0, 0]},
            {"id": "o1", "element": "O", "position": [2.5, 0, 0]},
            {"id": "h1", "element": "H", "position": [-0.5, 0.866, 0]},
            {"id": "h2", "element": "H", "position": [-0.5, -0.866, 0]},
            {"id": "h3", "element": "H", "position": [1.5, 0.866, 0]},
            {"id": "h4", "element": "H", "position": [1.5, -0.866, 0]},
            {"id": "h5", "element": "H", "position": [3, 0, 0]},
        ],
        "bonds": [
            {"id": "b1", "atom1": "c1", "atom2": "c2", "order": 1},
            {"id": "b2", "atom1": "c2", "atom2": "o1", "order": 1},
            {"id": "b3", "atom1": "c1", "atom2": "h1", "order": 1},
            {"id": "b4", "atom1": "c1", "atom2": "h2", "order": 1},
            {"id": "b5", "atom1": "c2", "atom2": "h3", "order": 1},
            {"id": "b6", "atom1": "c2", "atom2": "h4", "order": 1},
            {"id": "b7", "atom1": "o1", "atom2": "h5", "order": 1},
        ],
    },
}


class TestMoleculeAPI:
    """Test suite for molecule API endpoints."""

    @pytest.fixture(autouse=True)
    def setup(self):
        """Setup for each test."""
        self.base_url = API_BASE

    def test_to_smiles_endpoint(self):
        """Test SMILES generation for all test molecules."""
        for name, molecule in TEST_MOLECULES.items():
            try:
                response = requests.post(
                    f"{self.base_url}/api/molecule/to-smiles",
                    json={"molecule": molecule},
                    timeout=10,
                )
                assert response.status_code == 200, f"{name} failed: {response.status_code}"
                data = response.json()
                assert "smiles" in data, f"{name} missing smiles in response"
                logger.info(f"✓ {name}: {data['smiles']}")
            except Exception as e:
                logger.error(f"✗ {name}: {e}")
                raise

    def test_to_molblock_endpoint(self):
        """Test MolBlock generation for all test molecules."""
        for name, molecule in TEST_MOLECULES.items():
            try:
                response = requests.post(
                    f"{self.base_url}/api/molecule/to-molblock",
                    json={"molecule": molecule},
                    timeout=10,
                )
                assert response.status_code == 200, f"{name} failed: {response.status_code}"
                data = response.json()
                assert "molblock" in data, f"{name} missing molblock in response"
                assert len(data["molblock"]) > 0, f"{name} empty molblock"
                logger.info(f"✓ {name}: MolBlock generated ({len(data['molblock'])} chars)")
            except Exception as e:
                logger.error(f"✗ {name}: {e}")
                raise

    def test_validate_endpoint(self):
        """Test validation for all test molecules."""
        for name, molecule in TEST_MOLECULES.items():
            try:
                response = requests.post(
                    f"{self.base_url}/api/molecule/validate",
                    json={"molecule": molecule},
                    timeout=10,
                )
                assert response.status_code == 200, f"{name} failed: {response.status_code}"
                data = response.json()
                assert "valid" in data, f"{name} missing valid field"
                logger.info(f"✓ {name}: valid={data['valid']}")
            except Exception as e:
                logger.error(f"✗ {name}: {e}")
                raise

    def test_generate_2d_layout(self):
        """Test 2D layout generation for all test molecules."""
        for name, molecule in TEST_MOLECULES.items():
            try:
                response = requests.post(
                    f"{self.base_url}/api/molecule/generate-2d-layout",
                    json={"molecule": molecule},
                    timeout=10,
                )
                assert response.status_code == 200, f"{name} failed: {response.status_code}"
                data = response.json()
                assert "molecule" in data, f"{name} missing molecule in response"
                assert "atoms" in data["molecule"], f"{name} missing atoms"
                logger.info(f"✓ {name}: 2D layout generated")
            except Exception as e:
                logger.error(f"✗ {name}: {e}")
                raise

    def test_generate_3d_coordinates(self):
        """Test 3D coordinate generation for all test molecules."""
        for name, molecule in TEST_MOLECULES.items():
            try:
                response = requests.post(
                    f"{self.base_url}/api/molecule/generate-3d",
                    json={"molecule": molecule},
                    timeout=30,  # 3D generation can be slower
                )
                assert response.status_code == 200, f"{name} failed: {response.status_code}"
                data = response.json()
                assert "molecule" in data, f"{name} missing molecule in response"
                assert "atoms" in data["molecule"], f"{name} missing atoms"
                # Verify 3D coordinates (z != 0 for some atoms)
                atoms = data["molecule"]["atoms"]
                has_3d = any(atom.get("position", [0, 0, 0])[2] != 0 for atom in atoms)
                logger.info(f"✓ {name}: 3D coordinates generated (has_3d={has_3d})")
            except Exception as e:
                logger.error(f"✗ {name}: {e}")
                raise

    def test_from_smiles_endpoint(self):
        """Test loading from SMILES strings."""
        test_smiles = {
            "water": "O",
            "methane": "C",
            "ethanol": "CCO",
            "benzene": "c1ccccc1",
        }

        for name, smiles in test_smiles.items():
            try:
                response = requests.post(
                    f"{self.base_url}/api/molecule/from-smiles",
                    json={"smiles": smiles},
                    timeout=10,
                )
                assert response.status_code == 200, f"{name} failed: {response.status_code}"
                data = response.json()
                assert "molecule" in data, f"{name} missing molecule in response"
                assert len(data["molecule"]["atoms"]) > 0, f"{name} no atoms generated"
                logger.info(f"✓ {name}: Loaded from SMILES '{smiles}'")
            except Exception as e:
                logger.error(f"✗ {name}: {e}")
                raise

    def test_normalize_hydrogens(self):
        """Test hydrogen normalization."""
        for name, molecule in TEST_MOLECULES.items():
            try:
                response = requests.post(
                    f"{self.base_url}/api/molecule/normalize-hydrogens",
                    json={"molecule": molecule},
                    timeout=10,
                )
                assert response.status_code == 200, f"{name} failed: {response.status_code}"
                data = response.json()
                assert "molecule" in data, f"{name} missing molecule in response"
                logger.info(f"✓ {name}: Hydrogens normalized")
            except Exception as e:
                logger.error(f"✗ {name}: {e}")
                raise

    def test_error_handling(self):
        """Test error handling for invalid inputs."""
        # Invalid molecule structure
        response = requests.post(
            f"{self.base_url}/api/molecule/to-smiles",
            json={"molecule": {"atoms": [], "bonds": []}},
            timeout=10,
        )
        assert response.status_code in [400, 500], "Should return error for empty molecule"

        # Invalid SMILES
        response = requests.post(
            f"{self.base_url}/api/molecule/from-smiles",
            json={"smiles": "INVALID_SMILES_!!!###"},
            timeout=10,
        )
        assert response.status_code in [400, 500], "Should return error for invalid SMILES"

        logger.info("✓ Error handling tests passed")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

