# backend/core/security.py
"""
Security utilities (JWT, API keys, etc.)
Future implementation for authentication
"""
from typing import Optional
from datetime import datetime, timedelta
from jose import JWTError, jwt
from backend.config import settings


def create_access_token(data: dict, expires_delta: Optional[timedelta] = None):
    """
    Create JWT access token
    
    Args:
        data: Data to encode in token
        expires_delta: Optional expiration time
        
    Returns:
        Encoded JWT token
    """
    to_encode = data.copy()
    if expires_delta:
        expire = datetime.utcnow() + expires_delta
    else:
        expire = datetime.utcnow() + timedelta(minutes=settings.ACCESS_TOKEN_EXPIRE_MINUTES)
    
    to_encode.update({"exp": expire})
    encoded_jwt = jwt.encode(to_encode, settings.SECRET_KEY or "secret", algorithm=settings.ALGORITHM)
    return encoded_jwt


def verify_token(token: str) -> Optional[dict]:
    """
    Verify JWT token
    
    Args:
        token: JWT token to verify
        
    Returns:
        Decoded token data or None if invalid
    """
    try:
        payload = jwt.decode(token, settings.SECRET_KEY or "secret", algorithms=[settings.ALGORITHM])
        return payload
    except JWTError:
        return None

