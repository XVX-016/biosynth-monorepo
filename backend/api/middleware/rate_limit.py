"""
Rate Limiting Middleware - Token bucket algorithm for API endpoints

Limits request rate per IP and per user to prevent abuse.
"""

import time
from typing import Dict, Optional
from collections import defaultdict
from fastapi import Request, HTTPException, status
from starlette.middleware.base import BaseHTTPMiddleware
import logging

logger = logging.getLogger(__name__)


class TokenBucket:
    """
    Token bucket for rate limiting.
    
    Refills tokens at a constant rate up to capacity.
    """
    
    def __init__(self, capacity: int, refill_rate: float):
        """
        Args:
            capacity: Maximum tokens
            refill_rate: Tokens per second
        """
        self.capacity = capacity
        self.refill_rate = refill_rate
        self.tokens = float(capacity)
        self.last_refill = time.time()
    
    def consume(self, tokens: int = 1) -> bool:
        """
        Try to consume tokens.
        
        Returns:
            True if tokens were consumed, False if insufficient
        """
        self._refill()
        
        if self.tokens >= tokens:
            self.tokens -= tokens
            return True
        
        return False
    
    def _refill(self):
        """Refill tokens based on elapsed time."""
        now = time.time()
        elapsed = now - self.last_refill
        self.tokens = min(
            self.capacity,
            self.tokens + elapsed * self.refill_rate
        )
        self.last_refill = now
    
    def get_remaining(self) -> float:
        """Get remaining tokens."""
        self._refill()
        return self.tokens
    
    def get_reset_time(self) -> float:
        """Get time until bucket is full."""
        if self.tokens >= self.capacity:
            return 0.0
        
        needed = self.capacity - self.tokens
        return needed / self.refill_rate


class RateLimiter:
    """
    Rate limiter using token bucket per identifier (IP, user, etc.).
    """
    
    def __init__(
        self,
        capacity: int = 100,
        refill_rate: float = 10.0,  # 10 tokens per second
        default_limit: int = 100,
        default_window: int = 60,  # 60 seconds
    ):
        self.capacity = capacity
        self.refill_rate = refill_rate
        self.default_limit = default_limit
        self.default_window = default_window
        
        # Per-identifier buckets
        self.buckets: Dict[str, TokenBucket] = {}
        self.window_limits: Dict[str, Dict[str, Any]] = defaultdict(dict)
    
    def get_identifier(self, request: Request) -> str:
        """
        Get identifier for rate limiting (IP address or user ID).
        """
        # Try to get user ID from request (if authenticated)
        user_id = getattr(request.state, 'user_id', None)
        if user_id:
            return f"user:{user_id}"
        
        # Fall back to IP address
        client_host = request.client.host if request.client else "unknown"
        return f"ip:{client_host}"
    
    def is_allowed(
        self,
        identifier: str,
        tokens: int = 1,
        limit: Optional[int] = None,
        window: Optional[int] = None
    ) -> tuple[bool, Dict[str, Any]]:
        """
        Check if request is allowed.
        
        Returns:
            (allowed, info_dict)
        """
        # Use token bucket for general rate limiting
        if identifier not in self.buckets:
            limit = limit or self.default_limit
            window = window or self.default_window
            refill_rate = limit / window
            self.buckets[identifier] = TokenBucket(limit, refill_rate)
        
        bucket = self.buckets[identifier]
        allowed = bucket.consume(tokens)
        
        info = {
            'allowed': allowed,
            'remaining': bucket.get_remaining(),
            'reset_time': bucket.get_reset_time(),
        }
        
        return allowed, info
    
    def cleanup_old_buckets(self, max_age: int = 3600):
        """Remove buckets that haven't been used recently."""
        now = time.time()
        to_remove = []
        
        for identifier, bucket in self.buckets.items():
            if now - bucket.last_refill > max_age:
                to_remove.append(identifier)
        
        for identifier in to_remove:
            del self.buckets[identifier]


# Global rate limiter instance
_rate_limiter = RateLimiter()


class RateLimitMiddleware(BaseHTTPMiddleware):
    """
    FastAPI middleware for rate limiting.
    """
    
    def __init__(
        self,
        app,
        limit: int = 100,
        window: int = 60,
        exempt_paths: Optional[list] = None,
    ):
        super().__init__(app)
        self.limit = limit
        self.window = window
        self.exempt_paths = exempt_paths or []
    
    async def dispatch(self, request: Request, call_next):
        # Check if path is exempt
        if any(request.url.path.startswith(path) for path in self.exempt_paths):
            return await call_next(request)
        
        # Get identifier
        identifier = _rate_limiter.get_identifier(request)
        
        # Check rate limit
        allowed, info = _rate_limiter.is_allowed(
            identifier,
            tokens=1,
            limit=self.limit,
            window=self.window,
        )
        
        if not allowed:
            raise HTTPException(
                status_code=status.HTTP_429_TOO_MANY_REQUESTS,
                detail={
                    'error': 'Rate limit exceeded',
                    'reset_time': info['reset_time'],
                },
                headers={
                    'X-RateLimit-Limit': str(self.limit),
                    'X-RateLimit-Remaining': str(int(info['remaining'])),
                    'X-RateLimit-Reset': str(int(time.time() + info['reset_time'])),
                },
            )
        
        # Add rate limit headers to response
        response = await call_next(request)
        response.headers['X-RateLimit-Limit'] = str(self.limit)
        response.headers['X-RateLimit-Remaining'] = str(int(info['remaining']))
        response.headers['X-RateLimit-Reset'] = str(int(time.time() + info['reset_time']))
        
        return response


def get_rate_limiter() -> RateLimiter:
    """Get global rate limiter instance."""
    return _rate_limiter

