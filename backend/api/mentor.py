"""
Mentor API Routes
Chemistry mentor chat endpoints powered by Gemini
"""
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import Optional
import logging

from backend.services.mentor_service import get_mentor_service

logger = logging.getLogger(__name__)

router = APIRouter()


class MentorChatRequest(BaseModel):
    """Request to chat with chemistry mentor"""
    prompt: str
    mentor_skill: str  # "molecule-builder", "property-analyst", "reaction-planner", "reaction-simulator"
    conversation_id: Optional[str] = None


class MentorChatResponse(BaseModel):
    """Response from chemistry mentor"""
    response: str
    mentor_skill: str
    blocked: bool = False


@router.post("/chat", response_model=MentorChatResponse)
async def mentor_chat(request: MentorChatRequest):
    """
    Chat with chemistry mentor
    
    Enforces chemistry-only scope through:
    - Frontend validation (should be done before calling)
    - Backend keyword check
    - Gemini system prompt
    - Output filtering
    """
    try:
        mentor_service = get_mentor_service()
        result = await mentor_service.chat(request.prompt, request.mentor_skill)
        
        if result.get("blocked"):
            raise HTTPException(
                status_code=400,
                detail=result.get("error", "Request blocked: not chemistry-related")
            )
        
        if "error" in result:
            raise HTTPException(
                status_code=500,
                detail=result["error"]
            )
        
        return MentorChatResponse(
            response=result["response"],
            mentor_skill=result["mentor_skill"],
            blocked=result.get("blocked", False)
        )
    except ValueError as e:
        # Service not available (missing package or API key)
        raise HTTPException(
            status_code=503,
            detail=f"Mentor service unavailable: {str(e)}"
        )
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Mentor chat endpoint error: {e}", exc_info=True)
        raise HTTPException(
            status_code=500,
            detail=f"Internal server error: {str(e)}"
        )


@router.get("/health")
async def mentor_health():
    """Check if mentor service is available"""
    try:
        mentor_service = get_mentor_service()
        return {
            "status": "available",
            "service": "chemistry-mentor",
            "model": "gemini-1.5-pro"
        }
    except Exception as e:
        return {
            "status": "unavailable",
            "error": str(e)
        }

