from setuptools import setup, find_packages

setup(
    name="molforge-backend",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "fastapi",
        "uvicorn",
        "pydantic",
    ],
    extras_require={
        "test": [
            "pytest>=7.0.0",
        ],
    },
)
