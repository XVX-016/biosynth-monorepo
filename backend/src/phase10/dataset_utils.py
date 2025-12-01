"""
Dataset Utilities - Helper functions for managing generated molecules

Utilities for storing, loading, and managing datasets of generated molecules
and their associated rewards and properties.
"""

from typing import List, Dict, Optional, Any
from pathlib import Path
import json
import logging
from dataclasses import dataclass, asdict

logger = logging.getLogger(__name__)


@dataclass
class MoleculeRecord:
    """Record for a generated molecule with metadata."""
    smiles: str
    reward: float
    properties: Dict[str, float]
    generation_iteration: int
    generation_method: str  # "rl", "generative", etc.
    metadata: Dict[str, Any]


class DatasetUtils:
    """
    Utilities for managing generated molecule datasets.
    """
    
    def __init__(self, storage_path: Optional[str] = None):
        """
        Initialize dataset utilities.
        
        Args:
            storage_path: Optional path to store datasets
        """
        self.storage_path = Path(storage_path) if storage_path else Path("data/phase10")
        self.storage_path.mkdir(parents=True, exist_ok=True)
        self.records: List[MoleculeRecord] = []
        logger.info(f"DatasetUtils initialized (storage: {self.storage_path})")
    
    def add_record(
        self,
        smiles: str,
        reward: float,
        properties: Dict[str, float],
        generation_iteration: int,
        generation_method: str,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> MoleculeRecord:
        """
        Add a molecule record to the dataset.
        
        Args:
            smiles: SMILES string
            reward: Reward score
            properties: Property predictions
            generation_iteration: Iteration number
            generation_method: Generation method used
            metadata: Optional metadata
        
        Returns:
            Created record
        """
        record = MoleculeRecord(
            smiles=smiles,
            reward=reward,
            properties=properties,
            generation_iteration=generation_iteration,
            generation_method=generation_method,
            metadata=metadata or {},
        )
        self.records.append(record)
        return record
    
    def get_top_candidates(self, n: int = 10) -> List[MoleculeRecord]:
        """
        Get top N candidates by reward.
        
        Args:
            n: Number of top candidates
        
        Returns:
            List of top records
        """
        sorted_records = sorted(self.records, key=lambda r: r.reward, reverse=True)
        return sorted_records[:n]
    
    def filter_by_reward(self, min_reward: float) -> List[MoleculeRecord]:
        """
        Filter records by minimum reward.
        
        Args:
            min_reward: Minimum reward threshold
        
        Returns:
            Filtered records
        """
        return [r for r in self.records if r.reward >= min_reward]
    
    def filter_by_property(
        self,
        property_name: str,
        min_value: Optional[float] = None,
        max_value: Optional[float] = None,
    ) -> List[MoleculeRecord]:
        """
        Filter records by property value.
        
        Args:
            property_name: Property to filter on
            min_value: Minimum value (optional)
            max_value: Maximum value (optional)
        
        Returns:
            Filtered records
        """
        filtered = []
        for record in self.records:
            if property_name in record.properties:
                value = record.properties[property_name]
                if min_value is not None and value < min_value:
                    continue
                if max_value is not None and value > max_value:
                    continue
                filtered.append(record)
        return filtered
    
    def save_dataset(self, filename: str = "generated_molecules.json"):
        """
        Save dataset to file.
        
        Args:
            filename: Output filename
        """
        filepath = self.storage_path / filename
        data = [asdict(record) for record in self.records]
        with open(filepath, "w") as f:
            json.dump(data, f, indent=2)
        logger.info(f"Saved {len(self.records)} records to {filepath}")
    
    def load_dataset(self, filename: str = "generated_molecules.json") -> int:
        """
        Load dataset from file.
        
        Args:
            filename: Input filename
        
        Returns:
            Number of records loaded
        """
        filepath = self.storage_path / filename
        if not filepath.exists():
            logger.warning(f"Dataset file not found: {filepath}")
            return 0
        
        with open(filepath, "r") as f:
            data = json.load(f)
        
        self.records = [MoleculeRecord(**record) for record in data]
        logger.info(f"Loaded {len(self.records)} records from {filepath}")
        return len(self.records)
    
    def get_statistics(self) -> Dict[str, Any]:
        """
        Get dataset statistics.
        
        Returns:
            Statistics dictionary
        """
        if not self.records:
            return {
                "total": 0,
                "avg_reward": 0.0,
                "max_reward": 0.0,
                "min_reward": 0.0,
            }
        
        rewards = [r.reward for r in self.records]
        return {
            "total": len(self.records),
            "avg_reward": sum(rewards) / len(rewards),
            "max_reward": max(rewards),
            "min_reward": min(rewards),
            "by_method": self._count_by_method(),
        }
    
    def _count_by_method(self) -> Dict[str, int]:
        """Count records by generation method."""
        counts = {}
        for record in self.records:
            method = record.generation_method
            counts[method] = counts.get(method, 0) + 1
        return counts
    
    def clear(self):
        """Clear all records."""
        self.records = []
        logger.info("Dataset cleared")

