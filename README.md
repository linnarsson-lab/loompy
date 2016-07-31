# loompy

★ This repository is under construction, and not yet ready for public use. Be patient.

`.loom` is an efficient file format for very large omics datasets, 
consisting of a main matrix and a variable number of row and column 
annotations. We use loom files to store single-cell gene expression 
data: the main matrix contains the actual expression values (one 
column per cell, one row per gene); row and column annotations 
contain metadata for genes and cells, such as `Name`, `Chromosome`, 
`Position` (for genes), and `Strain`, `Sex`, `Age` (for cells).

Loom files (`.loom`) are created in the [HDF5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format) file format, which 
supports an internal collection of numerical multidimensional datasets.
HDF5 is supported by many computer languages, including Java, MATLAB, 
Mathematica, Python, R, and Julia. `.loom` files are accessible from 
any language that supports HDF5.

## Installation

Use [pip](https://pip.pypa.io/en/stable/) from your terminal:

```bash
pip install loompy
```

**Note:** there are some prerequisites, which will be installed along with loompy. If you 
use the popular [Anaconda](https://www.continuum.io/why-anaconda) Python distribution, all prerequisites will have
already been installed. 

## Getting started
 
```python
import loom
ds = loom.connect("cortex.loom")
print ds.row_attrs.keys()
```

This will print the names of all the row attribute in the file. 

## Understanding the semantics of loom files

### Connecting, not loading and saving

Loom files are stored on disk and are never loaded entirely. They
are more like databases: you connect, retrieve some subset of the data,
maybe update some attributes.

### Reading and writing

Loom files are based on [HDF5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format), a file format suitable for large multidimensional
datasets. They are designed to be mostly created once, then used as read-only. 
They **do not** support writing and reading concurrently. They also
do no support journalling, so if something happens during a write, the 
**entire file can be lost**. Therefore, do not use loom files as 
your primary data storage. They are for working with data, not keeping 
it safe.

Loom files are great for distribution of large datasets, which are then
used as read-only for analytical purposes.

### Efficient indexing

The main matrix is stored in *chunked* format. That is, instead of being
stored by rows or by columns, it is stored as a sequence of little rectangles. 
As a consequence, both rows and columns (as well as submatrices) can be efficiently 
accessed. 

## Documentation

### Creating `.loom` files

Create from data:

```python
def create(filename, matrix, row_attrs, col_attrs):
	"""
	Create a new .loom file from the given data.

	Args:
		filename (str):			The filename (typically using a '.loom' file extension)
		
		matrix (numpy.ndarray):	Two-dimensional (N-by-M) numpy ndarray of float values
		
		row_attrs (dict):		Row attributes, where keys are attribute names and values are numpy arrays (float or string) of length N
		
		col_attrs (dict):		Column attributes, where keys are attribute names and values are numpy arrays (float or string) of length M

	Returns:
		Nothing. To work with the file, use loom.connect(filename).
	"""
```

Create from an existing CEF file:

```python
def create_from_cef(cef_file, loom_file):
	"""
	Create a .loom file from a legacy CEF file.

	Args:
		cef_file (str):		filename of the input CEF file
		
		loom_file (str):	filename of the output .loom file (will be created)
	
	Returns:
		Nothing.
	"""
```

Create from a Pandas DataFrame:

```python
def create_from_pandas(df, loom_file):
	"""
	Create a .loom file from a Pandas DataFrame.

	Args:
		df (pd.DataFrame):	Pandas DataFrame
		
		loom_file (str):	Name of the output .loom file (will be created)

	Returns:
		Nothing.

	The DataFrame can contain multi-indexes on both axes, which will be turned into row and column attributes
	of the .loom file. The main matrix of the DataFrame must contain only float values. The datatypes of the
	attributes will be inferred as either float or string. 
	"""
```

### Connecting

Establish a connection to an existing `.loom` file:

```python
def connect(filename):
	"""
	Establish a connection to a .loom file.

	Args:
		filename (str):		Name of the .loom file to open

	Returns:
		A LoomConnection instance.
	"""
```

Example:

```python
ds = loom.connect("filename.loom")
```

In the rest of the documentation below, `ds` is assumed to be an instance of `LoomConnection` obtained by connecting to a `.loom` file.

### Shape, indexing and slicing

The `shape` property returns the row and column count as a tuple:

```python
>>> ds.shape
(100,2345)
```

The data stored in the main matrix can be retrieved by indexing and slicing. The following are supported:

* Indices: anything that can be converted to a Python long
* Slices (i.e. `:` or `0:10`)
* Lists of the rows/columns you want (i.e. `[0, 34, 576]`)
* Mask arrays (i.e. numpy array of bool indicating the rows/columns you want)

Lists and mask arrays are supported along one dimension at a time only. Note that performance will
be poor if you select many rows (columns) out of a large matrix. It may be better to load the
entire matrix and then perform the sub-selection in memory (using numpy slicing).

Since the main matrix is two-dimensional, two arguments are always needed. Examples:

```python

ds[:, :]          # Return the entire matrix
ds[0:10, 0:10]    # Return the 10x10 submatrix starting at row and column zero 
ds[99, :]         # Return the 100th row 
ds[:, 99]         # Return the 100th column
ds[[0,3,5], :]    # Return rows with index 0, 3 and 5
ds[:, bool_array] # Return columns where bool_array elements are True
```

