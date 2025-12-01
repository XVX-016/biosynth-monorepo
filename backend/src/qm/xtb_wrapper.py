"""
XTB Wrapper - Interface for xTB semi-empirical quantum chemistry

This is a placeholder for future xTB integration.
Currently returns mock results.
"""

from typing import List, Optional
import logging

from .qm_interfaces import QMEngineProtocol, QMResult
from .qm_engine import QMEngine

logger = logging.getLogger(__name__)


class XTBWrapper(QMEngine):
    """
    Wrapper for xTB semi-empirical quantum chemistry package.
    
    Currently uses mock engine. In the future, this would:
    1. Check if xTB is installed
    2. Write input files
    3. Run xTB executable
    4. Parse output files
    """
    
    def __init__(self):
        super().__init__()
        self.xtb_available = False
        self._check_xtb()
    
    def _check_xtb(self):
        """Check if xTB is available."""
        # TODO: Check for xTB executable
        # import shutil
        # self.xtb_available = shutil.which("xtb") is not None
        self.xtb_available = False
        logger.warning("xTB not available, using mock engine")
    
    async def compute_energy(
        self,
        smiles: str,
        coordinates: Optional[List[List[float]]] = None,
        method: str = "GFN2-xTB",  # xTB default method
        basis: str = ""  # xTB doesn't use basis sets
    ) -> QMResult:
        """
        Compute energy using xTB (or mock if unavailable).
        
        TODO: Implement real xTB integration:
        1. Convert SMILES to xyz coordinates
        2. Write xyz file
        3. Run: xtb molecule.xyz --gfn 2
        4. Parse energy from output
        """
        if not self.xtb_available:
            return await super().compute_energy(smiles, coordinates, method, basis)
        
        # TODO: Real xTB implementation
        logger.warning("Real xTB integration not yet implemented, using mock")
        return await super().compute_energy(smiles, coordinates, method, basis)
    
    async def optimize_geometry(
        self,
        smiles: str,
        initial_coordinates: Optional[List[List[float]]] = None,
        method: str = "GFN2-xTB",
        basis: str = ""
    ) -> QMResult:
        """Optimize geometry using xTB (or mock if unavailable)."""
        if not self.xtb_available:
            return await super().optimize_geometry(smiles, initial_coordinates, method, basis)
        
        # TODO: Real xTB implementation
        logger.warning("Real xTB integration not yet implemented, using mock")
        return await super().optimize_geometry(smiles, initial_coordinates, method, basis)

