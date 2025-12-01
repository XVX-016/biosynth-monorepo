"""
Phase 7: Search and Screening API routes

Provides endpoints for similarity search, substructure search, and screening.
"""

from fastapi import APIRouter, HTTPException, Query
from typing import Optional
import logging

import sys
from pathlib import Path

# Add backend to path for imports
backend_path = Path(__file__).parent.parent
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

from src.search import SearchEngine, LibraryLoader, ScreeningPipeline
from src.search.schemas import (
    SimilaritySearchRequest,
    SimilaritySearchResponse,
    SimilaritySearchResult,
    SubstructureSearchRequest,
    SubstructureSearchResponse,
    SubstructureSearchResult,
    ScreeningRequest,
    ScreeningResponse,
    ScreeningResult,
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
    
    Uses ECFP fingerprints and Tanimoto similarity.
    """
    try:
        results = _search_engine.similarity_search(
            smiles=smiles,
            k=k,
            threshold=threshold
        )
        
        return SimilaritySearchResponse(
            query_smiles=smiles,
            results=[
                SimilaritySearchResult(**result) for result in results
            ],
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
            results=[
                SubstructureSearchResult(**result) for result in results
            ],
            count=len(results)
        )
    except Exception as e:
        logger.error(f"Substructure search error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/screening/run", response_model=ScreeningResponse)
async def run_screening(request: ScreeningRequest):
    """
    Run screening pipeline with predicate.
    
    Supports property-threshold and custom predicates.
    """
    try:
        # Create predicate from config
        predicate = _screening_pipeline.create_custom_predicate(
            request.predicate_config
        )
        
        # Get library to screen
        if request.library_path:
            # Load from file
            loader = LibraryLoader(_screening_pipeline.index)
            loader.load_from_file(request.library_path)
            library = list(_screening_pipeline.index.fingerprints.keys())
        else:
            # Screen existing index
            library = list(_screening_pipeline.index.fingerprints.keys())
        
        # Run screening
        results = _screening_pipeline.batch_screen(
            library=library,
            predicate_fn=predicate,
            max_results=request.max_results
        )
        
        return ScreeningResponse(
            results=[
                ScreeningResult(**result) for result in results
            ],
            count=len(results),
            total_screened=len(library)
        )
    except Exception as e:
        logger.error(f"Screening error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/library/load")
async def load_library(
    directory: Optional[str] = None,
    pattern: str = "*.smi"
):
    """
    Load molecules from library files into search index.
    """
    try:
        loader = LibraryLoader(_search_engine.get_index())
        count = loader.load_from_directory(
            directory=directory,
            pattern=pattern
        )
        
        return {
            "loaded": count,
            "total_in_index": len(_search_engine.get_index().fingerprints)
        }
    except Exception as e:
        logger.error(f"Library load error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/library/stats")
async def get_library_stats():
    """Get statistics about loaded library."""
    try:
        stats = _search_engine.get_index().get_stats()
        return stats
    except Exception as e:
        logger.error(f"Stats error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))

