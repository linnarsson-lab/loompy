from setuptools import setup, find_packages

# First update the version in loompy/_version.py, then:

# cd loompy  (the root loompy folder, not the one inside!)
# rm -r dist   (otherwise twine will upload the oldest build!)
# python setup.py sdist
# twine upload dist/*

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
	url="https://github.com/linnarsson-lab/loompy",
	download_url=f"https://github.com/linnarsson-lab/loompy/archive/{__version__}.tar.gz",
)
