"""
Pathway scoring and ranking functions
"""
from typing import Dict, List, Any
from .utils import get_reagent_availability
from chem.energy.forcefield import calculate_total_energy

def score_pathway(pathway: Dict[str, Any]) -> Dict[str, Any]:
    """
    Score a retrosynthesis pathway based on multiple criteria
    
    Returns pathway with score and breakdown
    """
    steps = pathway.get("steps", [])
    total_steps = pathway.get("total_steps", len(steps))
    
    # Scoring factors
    step_score = _score_steps(total_steps)
    reagent_score = _score_reagents(pathway)
    energy_score = _score_energy_feasibility(pathway)
    yield_score = _score_yield(pathway)
    
    # Weighted total score (higher is better)
    total_score = (
        0.3 * step_score +
        0.3 * reagent_score +
        0.2 * energy_score +
        0.2 * yield_score
    )
    
    return {
        **pathway,
        "score": total_score,
        "score_breakdown": {
            "step_score": step_score,
            "reagent_score": reagent_score,
            "energy_score": energy_score,
            "yield_score": yield_score
        }
    }

def _score_steps(num_steps: int) -> float:
    """
    Score based on number of steps (fewer steps = higher score)
    """
    if num_steps == 0:
        return 0.0
    # Normalize: 1 step = 1.0, 5 steps = 0.2, 10+ steps = 0.0
    if num_steps <= 1:
        return 1.0
    elif num_steps <= 5:
        return 1.0 - (num_steps - 1) * 0.2
    else:
        return max(0.0, 1.0 - (num_steps - 5) * 0.1)

def _score_reagents(pathway: Dict[str, Any]) -> float:
    """
    Score based on reagent availability and cost
    """
    steps = pathway.get("steps", [])
    if not steps:
        return 0.0
    
    total_cost = 0.0
    unavailable_count = 0
    
    for step in steps:
        reaction = step.get("reaction")
        if not reaction:
            continue
        
        precursors = step.get("precursors", [])
        for precursor in precursors:
            reagent = precursor.get("reagent", "")
            if reagent:
                availability = get_reagent_availability(reagent)
                if not availability.get("available", False):
                    unavailable_count += 1
                total_cost += availability.get("cost", 100.0)
    
    # Score: lower cost and fewer unavailable = higher score
    cost_score = max(0.0, 1.0 - (total_cost / 100.0))  # Normalize to 0-1
    availability_score = max(0.0, 1.0 - (unavailable_count / len(steps))) if steps else 0.0
    
    return (cost_score + availability_score) / 2.0

def _score_energy_feasibility(pathway: Dict[str, Any]) -> float:
    """
    Score based on energy feasibility of reactions
    """
    steps = pathway.get("steps", [])
    if not steps:
        return 0.0
    
    feasible_count = 0
    
    for step in steps:
        molecule = step.get("molecule")
        if not molecule:
            continue
        
        try:
            energy = calculate_total_energy(molecule)["total_energy"]
            # Simple heuristic: reasonable energy range
            if -1000 < energy < 1000:
                feasible_count += 1
        except:
            pass
    
    return feasible_count / len(steps) if steps else 0.0

def _score_yield(pathway: Dict[str, Any]) -> float:
    """
    Score based on estimated overall yield
    """
    steps = pathway.get("steps", [])
    if not steps:
        return 0.0
    
    overall_yield = 1.0
    
    for step in steps:
        reaction = step.get("reaction")
        if reaction:
            yield_val = reaction.get("yield", 0.5)
            overall_yield *= yield_val
    
    return overall_yield

def rank_pathways(pathways: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Rank pathways by score (highest first)
    """
    scored = [score_pathway(p) for p in pathways]
    return sorted(scored, key=lambda x: x.get("score", 0.0), reverse=True)

