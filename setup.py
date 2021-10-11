#!/usr/bin/env python

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ribonorma",
    version="1.1",
    author="ChocoParrot (Hilbert L.), Kroppeb (Robbe P.)",
    author_email="lachocoparrot@gmail.com",
    description="Implements additional normalisation methods for RNA Sequencing.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Chokyotager/Ribonorma",
    project_urls={
        "Bug Tracker": "https://github.com/Chokyotager/Ribonorma/issues",
    },
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
    install_requires=[],
    scripts=[],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6"
)
