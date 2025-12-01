"""
Phase 7: Screening API endpoints

Provides endpoints for similarity search, library loading, and screening pipelines.
"""

from fastapi import APIRouter, HTTPException, UploadFile, File
from pydantic import BaseModel
from typing import List, Dict, Optional, Any
import logging

from chem.screening import (
    SimilaritySearchEngine,
    LibraryLoader,
    ScreeningPipeline,
    get_fingerprint_index,
)

logger = logging.getLogger(__name__)

router = APIRouter()


# Request models
class SimilaritySearchRequest(BaseModel):
    query_smiles: str
    threshold: float = 0.5
    max_results: int = 100
    radius: int = 2


class BatchSearchRequest(BaseModel):
    query_smiles_list: List[str]
    threshold: float = 0.5
    max_results_per_query: int = 10


class ScreeningRequest(BaseModel):
    query_smiles: str
    similarity_threshold: float = 0.5
    max_results: int = 1000
    filters: Optional[Dict[str, Any]] = None


class LoadLibraryRequest(BaseModel):
    molecules: List[Dict[str, Any]]
    id_key: str = 'id'
    smiles_key: str = 'smiles'
    name_key: str = 'name'


# Response models
class SearchResult(BaseModel):
    molecule_id: str
    similarity: float
    smiles: str
    name: str
    metadata: Dict[str, Any]


class SearchResponse(BaseModel):
    results: List[SearchResult]
    count: int
    query_smiles: str
    threshold: float


@router.post("/similarity", response_model=SearchResponse)
async def similarity_search(req: SimilaritySearchRequest):
    """
    Search for similar molecules by SMILES.
    """
    try:
        engine = SimilaritySearchEngine()
        results = engine.search_by_smiles(
            req.query_smiles,
            threshold=req.threshold,
            max_results=req.max_results,
            radius=req.radius
        )
        
        return SearchResponse(
            results=[
                SearchResult(**result) for result in results
            ],
            count=len(results),
            query_smiles=req.query_smiles,
            threshold=req.threshold,
        )
    except Exception as e:
        logger.error(f"Similarity search error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/similarity/batch")
async def batch_similarity_search(req: BatchSearchRequest):
    """
    Batch similarity search for multiple queries.
    """
    try:
        engine = SimilaritySearchEngine()
        results = engine.batch_search(
            req.query_smiles_list,
            threshold=req.threshold,
            max_results_per_query=req.max_results_per_query
        )
        
        return {
            "results": {
                query: [
                    SearchResult(**result).dict() for result in result_list
                ]
                for query, result_list in results.items()
            },
            "count": sum(len(r) for r in results.values()),
        }
    except Exception as e:
        logger.error(f"Batch search error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/similarity/knn")
async def k_nearest_neighbors(
    query_smiles: str,
    k: int = 10,
    radius: int = 2
):
    """
    Get k nearest neighbors.
    """
    try:
        engine = SimilaritySearchEngine()
        results = engine.get_k_nearest(query_smiles, k=k, radius=radius)
        
        return {
            "results": [SearchResult(**result).dict() for result in results],
            "count": len(results),
            "query_smiles": query_smiles,
            "k": k,
        }
    except Exception as e:
        logger.error(f"KNN search error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/screening")
async def screen_library(req: ScreeningRequest):
    """
    Run screening pipeline with filters.
    """
    try:
        pipeline = ScreeningPipeline()
        
        # Add filters if provided
        if req.filters:
            for filter_name, filter_params in req.filters.items():
                if filter_name == 'property':
                    pipeline.add_property_filter(
                        filter_params.get('name'),
                        min_value=filter_params.get('min'),
                        max_value=filter_params.get('max'),
                    )
                elif filter_name == 'smiles':
                    pipeline.add_smiles_filter(
                        filter_params.get('pattern'),
                        mode=filter_params.get('mode', 'contains'),
                    )
        
        results = pipeline.screen_by_similarity(
            req.query_smiles,
            similarity_threshold=req.similarity_threshold,
            max_results=req.max_results,
            apply_filters=True
        )
        
        return {
            "results": [SearchResult(**result).dict() for result in results],
            "count": len(results),
            "query_smiles": req.query_smiles,
        }
    except Exception as e:
        logger.error(f"Screening error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/library/load")
async def load_library(req: LoadLibraryRequest):
    """
    Load molecules into search index from JSON.
    """
    try:
        loader = LibraryLoader()
        count = loader.load_from_molecule_dicts(
            req.molecules,
            id_key=req.id_key,
            smiles_key=req.smiles_key,
            name_key=req.name_key
        )
        
        return {
            "loaded": count,
            "total_in_index": len(get_fingerprint_index().fingerprints),
        }
    except Exception as e:
        logger.error(f"Library load error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/library/load-file")
async def load_library_file(
    file: UploadFile = File(...),
    file_type: str = "smiles",  # "smiles", "json", "sdf"
    delimiter: str = "\t",
    smiles_column: int = 0,
    name_column: Optional[int] = None,
):
    """
    Load library from uploaded file.
    """
    try:
        # Save uploaded file temporarily
        import tempfile
        import os
        
        with tempfile.NamedTemporaryFile(delete=False, suffix=f".{file_type}") as tmp:
            content = await file.read()
            tmp.write(content)
            tmp_path = tmp.name
        
        try:
            loader = LibraryLoader()
            
            if file_type == "smiles":
                count = loader.load_from_smiles_file(
                    tmp_path,
                    delimiter=delimiter,
                    smiles_column=smiles_column,
                    name_column=name_column,
                )
            elif file_type == "json":
                count = loader.load_from_json(tmp_path)
            elif file_type == "sdf":
                count = loader.load_from_sdf(tmp_path)
            else:
                raise HTTPException(status_code=400, detail=f"Unsupported file type: {file_type}")
            
            return {
                "loaded": count,
                "total_in_index": len(get_fingerprint_index().fingerprints),
                "file_type": file_type,
            }
        finally:
            # Clean up temp file
            os.unlink(tmp_path)
    
    except Exception as e:
        logger.error(f"File load error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/library/stats")
async def get_library_stats():
    """
    Get statistics about loaded library.
    """
    try:
        index = get_fingerprint_index()
        stats = index.get_stats()
        
        return {
            "index_stats": stats,
            "num_molecules": stats["num_molecules"],
        }
    except Exception as e:
        logger.error(f"Stats error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.delete("/library/clear")
async def clear_library():
    """
    Clear all molecules from search index.
    """
    try:
        index = get_fingerprint_index()
        index.fingerprints.clear()
        index.inverted_index.clear()
        index.metadata.clear()
        
        return {"cleared": True}
    except Exception as e:
        logger.error(f"Clear error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))

