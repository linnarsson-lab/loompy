from setuptools import setup, find_packages

setup(
    name = "loompy",
    version = "0.1",
    packages = find_packages(),
    install_requires = [
        'scikit-learn>=0.17.1',
        'h5py>=2.6.0', 
        'numexpr>=2.6.1', 
        'pandas>=0.18.1', 
        'scipy>=0.18.0',
        'numpy>=1.11.1', 
    ],

    # metadata for upload to PyPI
    author = "Sten Linnarsson",
    author_email = "sten.linnarsson@ki.se",
    description = "Create and use .loom files from Python",
    license = "BSD",
    keywords = "loom omics transcriptomics bioinformatics",
    url = "https://github.com/linnarsson-lab/loompy",
)