# backend/core/dependencies.py
"""
FastAPI dependencies
"""
from sqlmodel import Session
from backend.db import engine


def get_db():
    """
    Database session dependency
    Yields a database session and closes it after use
    """
    with Session(engine) as session:
        yield session

