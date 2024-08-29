from pathlib import Path
from setuptools import find_packages, setup

# First update the version in loompy/_version.py, then:

# cd ~/code/loompy  (the directory where loompy resides)
# rm -r dist   (otherwise twine will upload the oldest build!)
# python setup.py sdist
# twine upload dist/*

# NOTE: Don't forget to update the release version at loompy.github.io (index.html)!

# pylint: disable=exec-used
__version__ = '0.0.0'
exec(open('loompy/_version.py').read())

setup(
	name="loompy",
	version=__version__,
	packages=find_packages(),
	python_requires='>=3.6,<3.12',
	install_requires=['h5py', 'numpy<1.24', 'scipy', 'setuptools', 'numba', 'click', "numpy-groupies"],
	entry_points='''
		[console_scripts]
		loompy=loompy.commands:cli
	''',
	# metadata for upload to PyPI
	author="Linnarsson Lab",
	author_email="sten.linnarsson@ki.se",
	description="Work with Loom files for single-cell RNA-seq data",
	long_description=Path("README.md").read_text(encoding="utf-8"),
	long_description_content_type="text/markdown",
	license="BSD",
	keywords="loom omics transcriptomics bioinformatics",
	url="https://github.com/linnarsson-lab/loompy",
	download_url=f"https://github.com/linnarsson-lab/loompy/archive/{__version__}.tar.gz",
)
