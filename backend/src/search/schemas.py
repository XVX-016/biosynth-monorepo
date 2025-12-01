"""
Pydantic schemas for search and screening API
"""

from pydantic import BaseModel, Field
from typing import List, Optional, Dict, Any


class SimilaritySearchResponse(BaseModel):
    """Response for similarity search"""
    query_smiles: str
    results: List[Dict[str, Any]]
    count: int


class SubstructureSearchResponse(BaseModel):
    """Response for substructure search"""
    query_smarts: str
    results: List[Dict[str, Any]]
    count: int


class ScreeningRequest(BaseModel):
    """Request for screening pipeline"""
    predicate_type: str = Field(..., description="Type: 'property_threshold'")
    predicate_config: Dict[str, Any] = Field(..., description="Predicate configuration")
    max_results: int = Field(1000, ge=1, le=100000, description="Maximum results")


class ScreeningResponse(BaseModel):
    """Response for screening"""
    results: List[Dict[str, Any]]
    count: int
    total_screened: int
