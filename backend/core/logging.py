# backend/core/logging.py
"""
Logging configuration
"""
import logging
import sys
from backend.config import settings

# Configure root logger
logging.basicConfig(
    level=getattr(logging, settings.LOG_LEVEL.upper()),
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
    ]
)

# Create logger for backend
logger = logging.getLogger("molforge")

