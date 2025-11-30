# backend/db.py
from sqlmodel import SQLModel, create_engine, Session
from sqlalchemy import text
from backend.config import settings
from backend.models.db.molecule import Molecule, Item

DATABASE_URL = settings.DATABASE_URL
engine = create_engine(DATABASE_URL, echo=False, connect_args={"check_same_thread": False})

def init_db():
    """Initialize database tables and run lightweight migrations"""
    SQLModel.metadata.create_all(engine)
    _run_migrations()

def get_session():
    """Get database session (deprecated - use core.dependencies.get_db instead)"""
    with Session(engine) as session:
        yield session


def _run_migrations():
    """
    Minimal migrations to keep existing SQLite databases compatible
    without requiring Alembic.
    """
    # Only needed for SQLite (current default)
    if engine.url.get_backend_name() != "sqlite":
        return

    with engine.begin() as connection:
        columns = connection.execute(text("PRAGMA table_info(molecule);")).fetchall()
        column_names = {row[1] for row in columns}

        if "molfile" not in column_names:
            connection.execute(text("ALTER TABLE molecule ADD COLUMN molfile TEXT"))

