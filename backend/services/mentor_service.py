"""
Mentor Service
Chemistry mentor powered by Gemini API with strict chemistry-only enforcement
"""
import os
import logging
from typing import Dict, Optional

logger = logging.getLogger(__name__)

# Try to import Gemini, but handle gracefully if not installed
try:
    import google.generativeai as genai
    GEMINI_AVAILABLE = True
except ImportError:
    GEMINI_AVAILABLE = False
    logger.warning("google-generativeai not installed. Mentor service will not work.")

CHEMISTRY_SYSTEM_PROMPT = """
You are MolForge Chemistry Mentor.

Scope:
- Chemistry only
- Molecules, reactions, synthesis, properties
- Organic, inorganic, physical chemistry
- Drug-like molecules and reasoning

Rules:
- If a question is not chemistry-related, refuse politely
- Explain step-by-step
- Use chemical terminology
- Never discuss politics, programming, or general topics

If asked something outside chemistry, reply:
"I'm here to help only with chemistry and molecular science."
"""


class MentorService:
    """Chemistry mentor service using Gemini API"""
    
    def __init__(self):
        if not GEMINI_AVAILABLE:
            raise ValueError("google-generativeai package not installed")
        
        api_key = os.getenv("GEMINI_KEY")
        if not api_key:
            raise ValueError("GEMINI_KEY not found in environment variables")
        
        genai.configure(api_key=api_key)
        self.model = genai.GenerativeModel(
            model_name="gemini-1.5-pro",
            system_instruction=CHEMISTRY_SYSTEM_PROMPT
        )
        logger.info("MentorService initialized with Gemini API")
    
    async def chat(self, prompt: str, mentor_skill: str) -> Dict:
        """
        Chat with chemistry mentor
        
        Args:
            prompt: User's question/prompt
            mentor_skill: Mentor skill type (molecule-builder, property-analyst, etc.)
        
        Returns:
            Dict with response or error
        """
        # Validate chemistry-only
        if not self._is_chemistry_prompt(prompt):
            return {
                "error": "I can help only with chemistry-related topics.",
                "blocked": True
            }
        
        # Add mentor skill context
        skill_context = self._get_skill_context(mentor_skill)
        full_prompt = f"{skill_context}\n\nUser: {prompt}"
        
        try:
            response = self.model.generate_content(full_prompt)
            text = response.text
            
            # Post-filter output
            if self._contains_blocked_content(text):
                return {
                    "error": "I can assist only with chemistry-related topics.",
                    "blocked": True
                }
            
            return {
                "response": text,
                "mentor_skill": mentor_skill,
                "blocked": False
            }
        except Exception as e:
            logger.error(f"Mentor chat error: {e}", exc_info=True)
            return {"error": str(e), "blocked": False}
    
    def _is_chemistry_prompt(self, prompt: str) -> bool:
        """Check if prompt is chemistry-related"""
        keywords = [
            "molecule", "reaction", "bond", "atom", "smiles",
            "synthesis", "chemistry", "pKa", "solubility",
            "mechanism", "functional group", "organic", "inorganic",
            "compound", "chemical", "molecular", "reagent", "catalyst",
            "substrate", "product", "reactant", "equilibrium", "thermodynamics"
        ]
        prompt_lower = prompt.lower()
        return any(kw in prompt_lower for kw in keywords)
    
    def _contains_blocked_content(self, text: str) -> bool:
        """Check if response contains blocked non-chemistry content"""
        blocked = ["politics", "finance", "dating", "religion", "programming", "code"]
        text_lower = text.lower()
        return any(b in text_lower for b in blocked)
    
    def _get_skill_context(self, skill: str) -> str:
        """Get context prompt for specific mentor skill"""
        contexts = {
            "molecule-builder": "You are helping with molecule generation and 3D structure building. Focus on molecular structure, SMILES notation, and 3D coordinate generation.",
            "property-analyst": "You are helping with molecular property prediction and QSAR analysis. Focus on solubility, pKa, toxicity, drug-likeness, and other molecular properties.",
            "reaction-planner": "You are helping with synthesis route planning and retrosynthesis. Focus on reaction pathways, synthetic strategies, and reaction mechanisms.",
            "reaction-simulator": "You are helping with reaction simulation and outcome prediction. Focus on reaction mechanisms, product prediction, and reaction conditions."
        }
        return contexts.get(skill, "You are a chemistry mentor helping with chemistry questions.")


# Global service instance (lazy initialization)
_mentor_service: Optional[MentorService] = None


def get_mentor_service() -> MentorService:
    """Get or create mentor service instance"""
    global _mentor_service
    if _mentor_service is None:
        _mentor_service = MentorService()
    return _mentor_service

