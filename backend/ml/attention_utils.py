"""
Attention weight normalization and processing utilities.
"""

import numpy as np
from typing import List, Union


def minmax_normalize(arr: Union[List[float], np.ndarray]) -> List[float]:
    """
    Min-max normalize attention weights to [0, 1].
    
    Args:
        arr: Attention weights array
    
    Returns:
        Normalized array
    """
    a = np.array(arr, dtype=float)
    lo, hi = a.min(), a.max()
    
    if hi - lo < 1e-12:
        return (a * 0.0).tolist()  # All equal -> zeros
    
    return ((a - lo) / (hi - lo)).tolist()


def softmax_per_edge(arr: Union[List[float], np.ndarray]) -> List[float]:
    """
    Softmax normalization for attention weights.
    
    Args:
        arr: Attention weights array
    
    Returns:
        Softmax-normalized array (sums to 1)
    """
    a = np.array(arr, dtype=float)
    
    # Subtract max for numerical stability
    a_shifted = a - a.max()
    e = np.exp(a_shifted)
    
    return (e / e.sum()).tolist()


def aggregate_heads(attentions: np.ndarray, method: str = "mean") -> np.ndarray:
    """
    Aggregate attention across heads.
    
    Args:
        attentions: Array of shape [E, heads] or [heads, E]
        method: "mean" or "max"
    
    Returns:
        Aggregated attention [E]
    """
    if attentions.ndim == 1:
        return attentions
    
    if attentions.ndim == 2:
        if attentions.shape[0] < attentions.shape[1]:
            # [heads, E] -> [E, heads]
            attentions = attentions.T
        
        if method == "mean":
            return attentions.mean(axis=1)
        elif method == "max":
            return attentions.max(axis=1)
        else:
            return attentions.mean(axis=1)
    
    return attentions.flatten()


def aggregate_layers(
    attentions_list: List[np.ndarray],
    method: str = "last"
) -> np.ndarray:
    """
    Aggregate attention across layers.
    
    Args:
        attentions_list: List of attention arrays (one per layer)
        method: "last", "mean", "max", or "weighted"
    
    Returns:
        Aggregated attention
    """
    if not attentions_list:
        return np.array([])
    
    if method == "last":
        return attentions_list[-1]
    elif method == "mean":
        return np.mean(attentions_list, axis=0)
    elif method == "max":
        return np.max(attentions_list, axis=0)
    elif method == "weighted":
        # Weight later layers more
        weights = np.linspace(0.1, 1.0, len(attentions_list))
        weights = weights / weights.sum()
        return np.average(attentions_list, axis=0, weights=weights)
    else:
        return attentions_list[-1]


def threshold_attentions(
    attentions: np.ndarray,
    threshold: float = 0.1
) -> np.ndarray:
    """
    Threshold low attention values to zero.
    
    Args:
        attentions: Attention weights
        threshold: Minimum value to keep
    
    Returns:
        Thresholded attentions
    """
    a = np.array(attentions)
    a[a < threshold] = 0.0
    return a

