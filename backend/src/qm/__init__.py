"""
Phase 8: Quantum Chemistry Engine

Provides interfaces and mock implementations for QM calculations.
"""

from .qm_interfaces import QMEngineProtocol, QMResult
from .qm_engine import QMEngine
from .psi4_wrapper import Psi4Wrapper
from .xtb_wrapper import XTBWrapper

__all__ = [
    'QMEngineProtocol',
    'QMResult',
    'QMEngine',
    'Psi4Wrapper',
    'XTBWrapper',
]

