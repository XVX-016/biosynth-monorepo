"""
Pydantic schemas for search and screening API
"""

from pydantic import BaseModel, Field
from typing import List, Optional, Dict, Any, Callable


class SimilaritySearchRequest(BaseModel):
    """Request for similarity search"""
    smiles: str = Field(..., description="Query SMILES string")
    k: int = Field(10, ge=1, le=1000, description="Number of results to return")
    threshold: Optional[float] = Field(None, ge=0.0, le=1.0, description="Minimum similarity threshold")


class SimilaritySearchResult(BaseModel):
    """Single similarity search result"""
    molecule_id: str
    smiles: str
    similarity: float
    metadata: Dict[str, Any] = Field(default_factory=dict)


class SimilaritySearchResponse(BaseModel):
    """Response for similarity search"""
    query_smiles: str
    results: List[SimilaritySearchResult]
    count: int


class SubstructureSearchRequest(BaseModel):
    """Request for substructure search"""
    smarts: str = Field(..., description="SMARTS pattern")
    max_results: int = Field(100, ge=1, le=10000, description="Maximum number of results")


class SubstructureSearchResult(BaseModel):
    """Single substructure search result"""
    molecule_id: str
    smiles: str
    matched_atoms: List[int] = Field(default_factory=list)
    metadata: Dict[str, Any] = Field(default_factory=dict)


class SubstructureSearchResponse(BaseModel):
    """Response for substructure search"""
    query_smarts: str
    results: List[SubstructureSearchResult]
    count: int


class ScreeningRequest(BaseModel):
    """Request for screening pipeline"""
    library_path: Optional[str] = Field(None, description="Path to library file")
    predicate_type: str = Field(..., description="Type of predicate: 'property_threshold' or 'custom'")
    predicate_config: Dict[str, Any] = Field(..., description="Predicate configuration")
    max_results: int = Field(1000, ge=1, le=100000, description="Maximum results to return")


class ScreeningResult(BaseModel):
    """Single screening result"""
    molecule_id: str
    smiles: str
    passed: bool
    score: Optional[float] = None
    metadata: Dict[str, Any] = Field(default_factory=dict)


class ScreeningResponse(BaseModel):
    """Response for screening"""
    results: List[ScreeningResult]
    count: int
    total_screened: int

