from setuptools import setup, find_packages

setup(
    name = "loompy",
    version = "0.6.1",
    packages = find_packages(),
    install_requires = [
        'scikit-learn',
        'h5py', 
        'pandas', 
        'scipy',
        'numpy', 
    ],

    # metadata for upload to PyPI
    author = "Sten Linnarsson",
    author_email = "sten.linnarsson@ki.se",
    description = "Create and use .loom files from Python",
    license = "BSD",
    keywords = "loom omics transcriptomics bioinformatics",
    url = "https://github.com/linnarsson-lab/loompy",
)