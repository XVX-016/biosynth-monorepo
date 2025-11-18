"""
Main entrypoint for MolForge backend
Can be used with: uvicorn main:app --reload
"""
from app import app

__all__ = ["app"]

