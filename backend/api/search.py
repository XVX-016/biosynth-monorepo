"""
Phase 7: Search and Screening API routes

Provides endpoints for similarity search, substructure search, and screening.
"""

from fastapi import APIRouter, HTTPException, Query
from typing import Optional
import logging
import sys
from pathlib import Path

# Add backend to path
backend_path = Path(__file__).parent.parent
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

from src.search import (
    SearchEngine,
    LibraryLoader,
    ScreeningPipeline,
    SimilaritySearchResponse,
    SubstructureSearchResponse,
    ScreeningRequest,
    ScreeningResponse,
)

logger = logging.getLogger(__name__)

router = APIRouter()

# Global instances
_search_engine = SearchEngine()
_screening_pipeline = ScreeningPipeline(_search_engine.get_index())


@router.get("/similarity", response_model=SimilaritySearchResponse)
async def similarity_search(
    smiles: str = Query(..., description="Query SMILES string"),
    k: int = Query(10, ge=1, le=1000, description="Number of results"),
    threshold: Optional[float] = Query(None, ge=0.0, le=1.0, description="Minimum similarity")
):
    """
    Search for similar molecules by SMILES.
    """
    try:
        results = _search_engine.similarity_search(
            smiles=smiles,
            k=k,
            threshold=threshold
        )
        
        return SimilaritySearchResponse(
            query_smiles=smiles,
            results=results,
            count=len(results)
        )
    except Exception as e:
        logger.error(f"Similarity search error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/substructure", response_model=SubstructureSearchResponse)
async def substructure_search(
    smarts: str = Query(..., description="SMARTS pattern"),
    max_results: int = Query(100, ge=1, le=10000, description="Maximum results")
):
    """
    Search for molecules containing SMARTS pattern.
    """
    try:
        results = _search_engine.substructure_search(
            smarts=smarts,
            max_results=max_results
        )
        
        return SubstructureSearchResponse(
            query_smarts=smarts,
            results=results,
            count=len(results)
        )
    except Exception as e:
        logger.error(f"Substructure search error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/library/load")
async def load_library(
    directory: str = Query("data/libraries", description="Directory path"),
    pattern: str = Query("*.smi", description="File pattern")
):
    """
    Load molecules from library files into search index.
    """
    try:
        loader = LibraryLoader(_search_engine.get_index())
        count = loader.load_from_directory(directory=directory, pattern=pattern)
        
        return {
            "loaded": count,
            "total_in_index": _search_engine.get_index().size()
        }
    except Exception as e:
        logger.error(f"Library load error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/library/stats")
async def get_library_stats():
    """Get statistics about loaded library."""
    try:
        index = _search_engine.get_index()
        return {
            "num_molecules": index.size(),
        }
    except Exception as e:
        logger.error(f"Stats error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))

