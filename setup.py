# rnamotifs2, https://github.com/grexor/rnamotifs2

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_desc = fh.read()
    long_desc = "\n".join(long_desc.split("\n")[1:])

setup(
    name="rnamotifs2",
    license="GNU General Public License v3.0",
    version=open("rnamotifs2/version", "rt").readlines()[0].strip(),
    packages=find_packages(include=["rnamotifs2", "rnamotifs2.*"]),
    description="cluster motif analysis for alternative splicing and polyadenylation",
    long_description=long_desc,
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    zip_safe=False,
    python_requires=">=3.8",
    author="Gregor Rot",
    author_email="gregor.rot@gmail.com",
    url="https://github.com/grexor/rnamotifs2",
    keywords=["rnamotifs2", "bioinformatics", "splicing", "RNA maps", "motifs"],
    include_package_data=True,
    package_data={
        "rnamotifs2": ["version"],
    },
    install_requires=["pybio", "numpy", "scipy", "matplotlib", "fisher"],
    extras_require={
        "lift": ["pyliftover"],  # only for comps/paper.bh/lift_hg19_to_ensembl115.py
    },
    entry_points={
        "console_scripts": [
            "rnamotifs2 = rnamotifs2.cli:main",
        ],
    },
)
