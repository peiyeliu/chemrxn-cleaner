# setup.py
from pathlib import Path
import re

from setuptools import setup, find_packages


def read_version():
    """Parse __version__ from chemrxn_cleaner/__init__.py without importing."""
    init_path = Path(__file__).parent / "chemrxn_cleaner" / "__init__.py"
    version_match = re.search(r'__version__\s*=\s*[\'"]([^\'"]+)[\'"]', init_path.read_text())
    if not version_match:
        raise RuntimeError("Unable to find __version__ in chemrxn_cleaner/__init__.py")
    return version_match.group(1)


setup(
    name="chemrxn-cleaner",
    version=read_version(),
    description="A lightweight toolkit for cleaning and standardizing organic reaction datasets.",
    author="Peiye Liu",
    author_email="peiyeliu.work@outlook.com",
    url="https://github.com/peiyeliu/chemrxn-cleaner",
    packages=find_packages(exclude=("tests", "examples")),
    python_requires=">=3.9",
    install_requires=[
        "rdkit",
        "pandas>=1.5.0",
        "tqdm>=4.64.0",
        "ord-schema>=0.3.0"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
    ]
)
