"""
Prediction Cache - In-memory and Redis-backed caching for predictions

Caches predictions by canonical SMILES to avoid redundant computation.
"""

import hashlib
import json
import time
from typing import Optional, Dict, Any
import logging

logger = logging.getLogger(__name__)

try:
    import redis
    REDIS_AVAILABLE = True
except ImportError:
    REDIS_AVAILABLE = False
    logger.warning("Redis not available, using in-memory cache only")


class PredictionCache:
    """
    Cache for molecular property predictions.
    
    Supports both in-memory and Redis backends.
    """
    
    def __init__(
        self,
        use_redis: bool = False,
        redis_host: str = 'localhost',
        redis_port: int = 6379,
        redis_db: int = 0,
        default_ttl: int = 3600,  # 1 hour
        max_memory_size: int = 1000,  # Max items in memory cache
    ):
        self.use_redis = use_redis and REDIS_AVAILABLE
        self.default_ttl = default_ttl
        self.max_memory_size = max_memory_size
        
        # In-memory cache
        self.memory_cache: Dict[str, Dict[str, Any]] = {}
        self.memory_timestamps: Dict[str, float] = {}
        
        # Redis client
        self.redis_client = None
        if self.use_redis:
            try:
                self.redis_client = redis.Redis(
                    host=redis_host,
                    port=redis_port,
                    db=redis_db,
                    decode_responses=True,
                )
                # Test connection
                self.redis_client.ping()
                logger.info("Redis cache connected")
            except Exception as e:
                logger.warning(f"Failed to connect to Redis: {e}, falling back to memory cache")
                self.use_redis = False
    
    def _get_cache_key(self, smiles: str, model_id: Optional[str] = None) -> str:
        """Generate cache key from SMILES and model ID."""
        # Canonicalize SMILES for consistent keys
        from .featurize import canonicalize_smiles
        canonical = canonicalize_smiles(smiles)
        
        key_data = {
            'smiles': canonical,
            'model_id': model_id or 'default',
        }
        
        key_str = json.dumps(key_data, sort_keys=True)
        return hashlib.sha256(key_str.encode()).hexdigest()
    
    def get(
        self,
        smiles: str,
        model_id: Optional[str] = None
    ) -> Optional[Dict[str, Any]]:
        """
        Get cached prediction.
        
        Returns:
            Cached prediction dict or None if not found/expired
        """
        cache_key = self._get_cache_key(smiles, model_id)
        
        if self.use_redis and self.redis_client:
            try:
                cached = self.redis_client.get(cache_key)
                if cached:
                    return json.loads(cached)
            except Exception as e:
                logger.error(f"Redis get error: {e}")
        
        # Fallback to memory cache
        if cache_key in self.memory_cache:
            # Check TTL
            if time.time() - self.memory_timestamps[cache_key] < self.default_ttl:
                return self.memory_cache[cache_key]
            else:
                # Expired, remove
                del self.memory_cache[cache_key]
                del self.memory_timestamps[cache_key]
        
        return None
    
    def set(
        self,
        smiles: str,
        prediction: Dict[str, Any],
        model_id: Optional[str] = None,
        ttl: Optional[int] = None
    ):
        """
        Cache a prediction.
        
        Args:
            smiles: SMILES string
            prediction: Prediction result dict
            model_id: Model identifier
            ttl: Time to live in seconds (default: self.default_ttl)
        """
        cache_key = self._get_cache_key(smiles, model_id)
        ttl = ttl or self.default_ttl
        
        if self.use_redis and self.redis_client:
            try:
                self.redis_client.setex(
                    cache_key,
                    ttl,
                    json.dumps(prediction)
                )
            except Exception as e:
                logger.error(f"Redis set error: {e}")
        
        # Also store in memory cache
        self.memory_cache[cache_key] = prediction
        self.memory_timestamps[cache_key] = time.time()
        
        # Trim memory cache if too large
        if len(self.memory_cache) > self.max_memory_size:
            # Remove oldest entries
            sorted_items = sorted(
                self.memory_timestamps.items(),
                key=lambda x: x[1]
            )
            to_remove = len(self.memory_cache) - self.max_memory_size
            for i in range(to_remove):
                key = sorted_items[i][0]
                del self.memory_cache[key]
                del self.memory_timestamps[key]
    
    def invalidate(self, smiles: str, model_id: Optional[str] = None):
        """Invalidate cache entry."""
        cache_key = self._get_cache_key(smiles, model_id)
        
        if self.use_redis and self.redis_client:
            try:
                self.redis_client.delete(cache_key)
            except Exception as e:
                logger.error(f"Redis delete error: {e}")
        
        if cache_key in self.memory_cache:
            del self.memory_cache[cache_key]
            del self.memory_timestamps[cache_key]
    
    def clear(self):
        """Clear all cache entries."""
        if self.use_redis and self.redis_client:
            try:
                # Clear all keys matching pattern (if using namespaced keys)
                # For now, just clear memory cache
                pass
            except Exception as e:
                logger.error(f"Redis clear error: {e}")
        
        self.memory_cache.clear()
        self.memory_timestamps.clear()
    
    def get_stats(self) -> Dict[str, Any]:
        """Get cache statistics."""
        stats = {
            'backend': 'redis' if self.use_redis else 'memory',
            'memory_size': len(self.memory_cache),
            'max_memory_size': self.max_memory_size,
        }
        
        if self.use_redis and self.redis_client:
            try:
                info = self.redis_client.info('memory')
                stats['redis_memory'] = info.get('used_memory_human', 'N/A')
            except:
                pass
        
        return stats

