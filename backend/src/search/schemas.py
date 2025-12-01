"""
Pydantic schemas for search and screening API
"""

from pydantic import BaseModel, Field
from typing import List, Optional, Dict, Any


class SimilaritySearchRequest(BaseModel):
    """Request for similarity search"""
    query_smiles: str = Field(..., description="Query SMILES string")
    k: int = Field(10, ge=1, le=1000, description="Number of results")
    threshold: float = Field(0.0, ge=0.0, le=1.0, description="Minimum similarity threshold")


class SimilaritySearchResult(BaseModel):
    """Single similarity search result"""
    smiles: str
    similarity: float
    index: Optional[int] = None


class SimilaritySearchResponse(BaseModel):
    """Response for similarity search"""
    query_smiles: str
    results: List[Dict[str, Any]]
    count: int


class SubstructureSearchRequest(BaseModel):
    """Request for substructure search"""
    query_smarts: str = Field(..., description="Query SMARTS pattern")
    max_results: int = Field(1000, ge=1, le=100000, description="Maximum results")


class SubstructureSearchResult(BaseModel):
    """Single substructure search result"""
    smiles: str
    index: Optional[int] = None


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


class ScreeningResult(BaseModel):
    """Single screening result"""
    smiles: str
    properties: Dict[str, Any]
    index: Optional[int] = None


class ScreeningResponse(BaseModel):
    """Response for screening"""
    results: List[Dict[str, Any]]
    count: int
    total_screened: int
