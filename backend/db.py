# backend/db.py
from sqlmodel import SQLModel, create_engine, Session
from backend.config import settings
from backend.models.db.molecule import Molecule, Item

DATABASE_URL = settings.DATABASE_URL
engine = create_engine(DATABASE_URL, echo=False, connect_args={"check_same_thread": False})

def init_db():
    """Initialize database tables"""
    SQLModel.metadata.create_all(engine)

def get_session():
    """Get database session (deprecated - use core.dependencies.get_db instead)"""
    with Session(engine) as session:
        yield session

