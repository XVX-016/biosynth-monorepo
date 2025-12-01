"""
Model Pool - LRU-based model loading and unloading

Manages model instances in memory with automatic unloading of least recently used models.
"""

import time
from typing import Dict, Optional, Any
from collections import OrderedDict
import threading
import logging

logger = logging.getLogger(__name__)


class ModelPool:
    """
    Pool for managing loaded ML models with LRU eviction.
    """
    
    def __init__(
        self,
        max_models: int = 3,
        max_idle_time: int = 300,  # 5 minutes
        unload_on_idle: bool = True,
    ):
        """
        Args:
            max_models: Maximum number of models to keep in memory
            max_idle_time: Seconds before idle model is unloaded
            unload_on_idle: Whether to automatically unload idle models
        """
        self.max_models = max_models
        self.max_idle_time = max_idle_time
        self.unload_on_idle = unload_on_idle
        
        # LRU cache: OrderedDict maintains insertion order
        # Most recently used at end, least recently used at beginning
        self.models: OrderedDict[str, Dict[str, Any]] = OrderedDict()
        self.lock = threading.Lock()
        
        # Background cleanup thread
        if self.unload_on_idle:
            self._start_cleanup_thread()
    
    def get(self, model_id: str) -> Optional[Any]:
        """
        Get model instance, moving it to end (most recently used).
        
        Returns:
            Model instance or None if not loaded
        """
        with self.lock:
            if model_id in self.models:
                # Move to end (most recently used)
                model_info = self.models.pop(model_id)
                self.models[model_id] = model_info
                model_info['last_used'] = time.time()
                return model_info['model']
            return None
    
    def put(self, model_id: str, model: Any, metadata: Optional[Dict] = None):
        """
        Add model to pool.
        
        Args:
            model_id: Model identifier
            model: Model instance
            metadata: Optional metadata dict
        """
        with self.lock:
            # If already exists, update
            if model_id in self.models:
                self.models[model_id]['model'] = model
                self.models[model_id]['last_used'] = time.time()
                if metadata:
                    self.models[model_id]['metadata'] = metadata
                # Move to end
                self.models.move_to_end(model_id)
                return
            
            # Check if we need to evict
            if len(self.models) >= self.max_models:
                self._evict_lru()
            
            # Add new model
            self.models[model_id] = {
                'model': model,
                'last_used': time.time(),
                'loaded_at': time.time(),
                'metadata': metadata or {},
            }
    
    def remove(self, model_id: str) -> bool:
        """
        Remove model from pool.
        
        Returns:
            True if removed, False if not found
        """
        with self.lock:
            if model_id in self.models:
                model_info = self.models.pop(model_id)
                # Cleanup if model has cleanup method
                if hasattr(model_info['model'], 'cleanup'):
                    try:
                        model_info['model'].cleanup()
                    except Exception as e:
                        logger.error(f"Error cleaning up model {model_id}: {e}")
                return True
            return False
    
    def clear(self):
        """Clear all models from pool."""
        with self.lock:
            for model_id, model_info in self.models.items():
                if hasattr(model_info['model'], 'cleanup'):
                    try:
                        model_info['model'].cleanup()
                    except Exception as e:
                        logger.error(f"Error cleaning up model {model_id}: {e}")
            self.models.clear()
    
    def _evict_lru(self):
        """Evict least recently used model."""
        if not self.models:
            return
        
        # Remove first item (least recently used)
        model_id, model_info = self.models.popitem(last=False)
        logger.info(f"Evicting model {model_id} from pool")
        
        # Cleanup if model has cleanup method
        if hasattr(model_info['model'], 'cleanup'):
            try:
                model_info['model'].cleanup()
            except Exception as e:
                logger.error(f"Error cleaning up evicted model {model_id}: {e}")
    
    def _cleanup_idle(self):
        """Remove models that have been idle too long."""
        if not self.unload_on_idle:
            return
        
        now = time.time()
        to_remove = []
        
        with self.lock:
            for model_id, model_info in self.models.items():
                idle_time = now - model_info['last_used']
                if idle_time > self.max_idle_time:
                    to_remove.append(model_id)
        
        for model_id in to_remove:
            logger.info(f"Unloading idle model {model_id}")
            self.remove(model_id)
    
    def _start_cleanup_thread(self):
        """Start background thread for cleanup."""
        def cleanup_loop():
            while True:
                time.sleep(60)  # Check every minute
                try:
                    self._cleanup_idle()
                except Exception as e:
                    logger.error(f"Error in cleanup thread: {e}")
        
        thread = threading.Thread(target=cleanup_loop, daemon=True)
        thread.start()
    
    def get_stats(self) -> Dict[str, Any]:
        """Get pool statistics."""
        with self.lock:
            return {
                'loaded_models': len(self.models),
                'max_models': self.max_models,
                'models': [
                    {
                        'id': model_id,
                        'last_used': model_info['last_used'],
                        'loaded_at': model_info['loaded_at'],
                        'idle_time': time.time() - model_info['last_used'],
                    }
                    for model_id, model_info in self.models.items()
                ],
            }
    
    def __len__(self) -> int:
        """Get number of loaded models."""
        return len(self.models)
    
    def __contains__(self, model_id: str) -> bool:
        """Check if model is in pool."""
        return model_id in self.models

