"""
Phase 8: Molecular Dynamics Engine

Provides interfaces and mock implementations for MD simulations.
"""

from .md_interfaces import MDEngineProtocol, MDResult, TrajectoryFrame
from .md_engine import MDEngine
from .forcefields import ForceField, SimpleForceField
from .integrators import Integrator, VerletIntegrator

__all__ = [
    'MDEngineProtocol',
    'MDResult',
    'TrajectoryFrame',
    'MDEngine',
    'ForceField',
    'SimpleForceField',
    'Integrator',
    'VerletIntegrator',
]

