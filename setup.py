# setup.py
from setuptools import setup, find_packages

setup(
    name="chemrxn-cleaner",
    version="0.0.2",
    description="A lightweight toolkit for cleaning and standardizing organic reaction datasets.",
    author="Peiye Liu",
    author_email="peiyeliu.work@outlook.com",
    url="https://github.com/peiyeliu/chemrxn-cleaner",
    packages=find_packages(exclude=("tests", "examples")),
    python_requires=">=3.9",
    install_requires=[
        "rdkit-pypi",
        "pandas>=1.5.0",
        "tqdm>=4.64.0",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
    ],
    extras_require={
        "ord": ["ord-schema>=0.3.0"],
    }
)
