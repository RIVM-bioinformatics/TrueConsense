import sys

if sys.version_info.major != 3 or sys.version_info.minor < 8:
    print("Error: you must execute setup.py using Python 3.8 or later")
    sys.exit(1)

from setuptools import find_packages, setup

from TrueConsense.version import __version__

with open("README.md", "rb") as readme:
    DESCR = readme.read().decode()


setup(
    name="TrueConsense",
    version=__version__,
    author="Florian Zwagemaker, Dennis Schmitz",
    author_email="ids-bioinformatics@rivm.nl",
    license="AGPLv3",
    packages=find_packages(),
    install_requires=[
        "pysam=0.19",
        "pandas>=1.2.3",
        "gffpandas>=1.2.0",
        "tqdm>=4.59.0",
        "biopython>=1.78",
    ],
    entry_points={
        "console_scripts": [
            "trueconsense = TrueConsense.TrueConsense:main",
            "TrueConsense = TrueConsense.TrueConsense:main",
        ]
    },
    keywords=[],
    zip_safe=False,
)
