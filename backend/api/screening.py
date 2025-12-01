"""
Phase 7: Screening API routes

Provides endpoint for screening pipeline.
"""

from fastapi import APIRouter, HTTPException
import logging
import sys
from pathlib import Path

# Add backend to path
backend_path = Path(__file__).parent.parent
if str(backend_path) not in sys.path:
    sys.path.insert(0, str(backend_path))

from src.search import ScreeningPipeline, ScreeningRequest, ScreeningResponse

logger = logging.getLogger(__name__)

router = APIRouter()

# Global instance
_screening_pipeline = ScreeningPipeline()


@router.post("/run", response_model=ScreeningResponse)
async def run_screening(request: ScreeningRequest):
    """
    Run screening pipeline with predicate.
    """
    try:
        # Create predicate from config
        if request.predicate_type == 'property_threshold':
            config = request.predicate_config
            predicate = _screening_pipeline.property_threshold_predicate(
                name=config.get('name'),
                min_value=config.get('min'),
                max_value=config.get('max'),
            )
        else:
            raise HTTPException(
                status_code=400,
                detail=f"Unknown predicate type: {request.predicate_type}"
            )
        
        # Run screening
        results = _screening_pipeline.batch_screen(
            pred=predicate,
            max_results=request.max_results
        )
        
        total_screened = _screening_pipeline.get_index().size()
        
        return ScreeningResponse(
            results=results,
            count=len(results),
            total_screened=total_screened
        )
    except Exception as e:
        logger.error(f"Screening error: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))

