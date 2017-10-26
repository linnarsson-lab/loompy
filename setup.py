from setuptools import setup, find_packages

# pylint: disable=exec-used
__version__ = '0.0.0'
exec(open('loompy/_version.py').read())

setup(
    name="loompy",
    version=__version__,
    packages=find_packages(),
    install_requires=['h5py', 'numpy', 'scipy', "typing", "setuptools"],
    # metadata for upload to PyPI
    author="Linnarsson Lab",
    author_email="sten.linnarsson@ki.se",
    description="Work with .loom files for single-cell RNA-seq data",
    license="BSD",
    keywords="loom omics transcriptomics bioinformatics",
    url="https://github.com/linnarsson-lab/loompy"
)
