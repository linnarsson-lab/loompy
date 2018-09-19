from setuptools import find_packages, setup

# First update the version in loompy/_version.py, then:

# cd ~/code  (the directory where loompy resides)
# cp -R loompy loompy-3.5
# cd loompy-3.5/loompy
# for f in *.py; do py-backwards -i $f -o $f -t 3.5; done
# cd ..
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
	python_requires='>=3.6',
	install_requires=['h5py', 'numpy', 'scipy', 'setuptools', 'pandas'],
	extras_require=dict(colors=['matplotlib']),
	# metadata for upload to PyPI
	author="Linnarsson Lab",
	author_email="sten.linnarsson@ki.se",
	description="Work with Loom files for single-cell RNA-seq data",
	license="BSD",
	keywords="loom omics transcriptomics bioinformatics",
	url="https://github.com/linnarsson-lab/loompy",
	download_url=f"https://github.com/linnarsson-lab/loompy/archive/{__version__}.tar.gz",
)
