from setuptools import setup, find_packages

setup(
    name = "loom",
    version = "0.1",
    packages = find_packages(),
    install_requires = ['h5py>=2.6.0', 'numexpr>=2.6.1', 'pandas>=0.18.1', 'numpy>=1.11.1', 'scipy>=0.18.0', 'scikit-learn>=0.17.1']

    # metadata for upload to PyPI
    author = "Sten Linnarsson",
    author_email = "sten.linnarsson@ki.se",
    description = "Create and use .loom files from Python",
    license = "BSD",
    keywords = "loom omics transcriptomics bioinformatics",
    url = "https://github.com/linnarsson-lab/loom.py",
)