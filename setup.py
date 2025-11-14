# setup.py
from setuptools import setup, find_packages

setup(
    name="chemrxn-cleaner",  # 如果以后要发 PyPI，可以改成唯一名称
    version="0.0.1",
    description="A lightweight toolkit for cleaning and standardizing organic reaction datasets.",
    author="Your Name",
    author_email="you@example.com",
    url="https://github.com/your-github-username/chemrxn-cleaner",
    packages=find_packages(exclude=("tests", "examples")),
    python_requires=">=3.9",
    install_requires=[
        "rdkit-pypi>=2022.9.5",
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
)
