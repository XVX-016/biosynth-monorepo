# backend/core/exceptions.py
"""
Custom exceptions for MolForge backend
"""
from fastapi import HTTPException, status


class MolForgeException(Exception):
    """Base exception for MolForge"""
    pass


class InvalidSMILESError(MolForgeException):
    """Raised when SMILES string is invalid"""
    pass


class ModelLoadError(MolForgeException):
    """Raised when model fails to load"""
    pass


class PredictionError(MolForgeException):
    """Raised when prediction fails"""
    pass


def raise_http_exception(e: Exception, status_code: int = 500, detail: str = None):
    """Convert exception to HTTPException"""
    if detail is None:
        detail = str(e)
    raise HTTPException(status_code=status_code, detail=detail)

