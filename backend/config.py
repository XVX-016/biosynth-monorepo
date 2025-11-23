# backend/config.py
"""
Centralized configuration management
"""
import os
from dotenv import load_dotenv
from typing import Optional, List

# Load environment variables
load_dotenv()


class Settings:
    """Application settings"""
    
    # Database
    DATABASE_URL: str = os.getenv(
        "DATABASE_URL",
        "sqlite:///./biosynth.db"
    )
    
    # API
    API_TITLE: str = "MolForge Backend API"
    API_VERSION: str = "0.1.0"
    API_DESCRIPTION: str = "Backend API for MolForge molecular design platform"
    
    # CORS
    # Allow CORS origins from environment variable or default to localhost
    # Format in env: "http://localhost:5173,https://your-app.web.app"
    _cors_origins_env = os.getenv("CORS_ORIGINS", "")
    _is_dev = os.getenv("ENVIRONMENT", "development").lower() in ["development", "dev"]
    
    if _cors_origins_env:
        CORS_ORIGINS: List[str] = [origin.strip() for origin in _cors_origins_env.split(",") if origin.strip()]
    elif _is_dev:
        # In development, allow all localhost origins (more permissive)
        CORS_ORIGINS: List[str] = [
            "http://localhost:5173",
            "http://localhost:5174",
            "http://localhost:3000",
            "http://localhost:5175",
            "http://localhost:8080",
            "http://127.0.0.1:5173",
            "http://127.0.0.1:5174",
            "http://127.0.0.1:3000",
            "http://127.0.0.1:5175",
            "http://127.0.0.1:8080",
        ]
    else:
        # Production: only specific origins
        CORS_ORIGINS: List[str] = [
            "http://localhost:5173",
            "http://localhost:5174",
            "http://127.0.0.1:5173",
            "http://127.0.0.1:5174",
        ]
    
    # Model paths
    MODEL_WEIGHTS_PATH: str = os.getenv(
        "MODEL_WEIGHTS_PATH",
        "backend/weights/property_predictor.pt"
    )
    
    ONNX_MODEL_PATH: str = os.getenv(
        "ONNX_MODEL_PATH",
        "backend/weights/property_predictor.onnx"
    )
    
    # Security (future)
    SECRET_KEY: Optional[str] = os.getenv("SECRET_KEY")
    ALGORITHM: str = "HS256"
    ACCESS_TOKEN_EXPIRE_MINUTES: int = 30
    
    # Logging
    LOG_LEVEL: str = os.getenv("LOG_LEVEL", "INFO")


# Global settings instance
settings = Settings()

