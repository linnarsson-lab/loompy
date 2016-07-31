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

When you connect, all attributes are read into memory for quick access, but the 
main matrix remains on disk.

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

### Matrix and attributes

The data model that loom supports is a two-dimensional matrix with named row and column attributes. 
If the main matrix has N rows and M columns, then each row attribute has M elements (one per column) 
and vice versa.

By convention, the matrix represents measurements of some property (e.g. gene expression), genes are
stored in the rows and columns represent samples (e.g. cells).

The file can grow by the addition of attributes and columns (but not rows). For this reason, it's 
a good idea to put genes in the rows (since you will likely always work with the same  gene set).

### Data types

Loom supports a tiny subset of numpy data types:

* The main matrix is always a two-dimensional array of type `float32`
* Attributes are one-dimensional arrays of either `float64` or `string`

Note that there is no integer attribute type. However, float64s are large enough to 
represent all integers up to and including 9,007,199,254,740,992 without loss.

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

Note: there is usually no need to close the connection. The exception is if you need to write to the loom file from
two different processes (sequentially, not simultaneously). In that case, the first process needs to let go of the 
file by calling `close()` on the connection, before the second can start writing:

```python
def close(self):
    """
    Close the connection. After this, the connection object becomes invalid.
    """
```


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

### Attributes

Attributes are accessed as dictionaries on `row_attrs` and `col_attrs`, respectively. For example:

```python
ds.row_attrs.keys()       # Return list of row attribute names
ds.col_attrs.keys()       # Return list of column attribute names
ds.row_attrs["GeneName"]  # Return a numpy array of gene names (assuming the attirbute exists)
```

Note that these dictionaries are **read-only**. Do not attempt to change the value of an attribute (see below for
how to add and modify attributes).

### Adding attributes and columns

You can add attributes and columns to an existing loom file. It is not possible to add rows or to delete 
attributes or any part of the matrix.

To add an attribute, which also saves it to the loom file:

```python
def set_attr(self, name, values, axis = 0):
    """
    Create or modify an attribute.

    Args:
        name (str): 			Name of the attribute
        
        values (numpy.ndarray):	Array of values of length equal to the axis length
        
        axis (int):				Axis of the attribute (0 = rows, 1 = columns)

    Returns:
        Nothing.

    This will overwrite any existing attribute of the same name.
    """
```

**Note:** If you use an existing attribute name, the existing attribute will be overwritten.

You can also add an attribute by providing a lookup dictionary based on an existing attribute. 
This can be very useful to fill in missing values, fix typos etc:


```python
def set_attr_bydict(self, name, fromattr, dict, axis = 0, default = None):
    """
    Create or modify an attribute by mapping source values to target values.

    Args:
        name (str): 			Name of the destination attribute
        
        fromattr (str):			Name of the source attribute
        
        dict (dict):			Key-value mapping from source to target values
        
        axis (int):				Axis of the attribute (0 = rows, 1 = columns)
        
        default: (float or str):	Default target value to use if no mapping exists for a source value

    Returns:
        Nothing.

    This will overwrite any existing attribute of the same name. It is perfectly possible to map an
    attribute to itself (in-place).
    """
```

To add columns:

```python
def add_columns(self, submatrix, col_attrs):
    """
    Add columns of data and attribute values to the dataset.

    Args:
        submatrix (numpy.ndarray):	An N-by-M matrix of floats (N rows, M columns)
        
        col_attrs (dict):			Column attributes, where keys are attribute names and values are numpy arrays (float or string) of length M

    Returns:
        Nothing.

    Note that this will modify the underlying HDF5 file, which will interfere with any concurrent readers.
    """
```

You need to provide a submatrix corresponding to the columns, as well as a dictionary of column attributes
with values for all the new columns. 

**Note:** It is not possible to add rows. 


### Mathematical operations

#### Map

You can map a function across all rows (all columns), while avoiding loading the entire
dataset into memory:

```python
def map(self, f, axis = 0, chunksize = 100000000):
    """
    Apply a function along an axis without loading the entire dataset in memory.

    Args:
        f (func):		Function that takes a numpy ndarray as argument
        
        axis (int):		Axis along which to apply the function (0 = rows, 1 = columns)
        
        chunksize (int): Number of values to load per chunk

    Returns:
        numpy.ndarray result of function application
    """
```

The function will receive an array (of floats) as its only argument, and should return a single float value.


#### Correlation matrix

Numpy can compute correlation matrices but will internally cast float32 to float64. This leads to
unnecessarily large memory consumption. The following method computes the correlation matrix
while avoiding float64s:

```python
def corr_matrix(self, axis = 0, log=False):
    """
    Compute correlation matrix without casting to float64.

    Args:
        axis (int):			The axis along which to compute the correlation matrix.
        
        log (bool):			If true, compute correlation on log(x+1) values
    
    Returns:
        numpy.ndarray of float32 correlation coefficents

    This function avoids casting intermediate values to double (float64), to reduce memory footprint.
    If row attribute _Excluded exists, those rows will be excluded.
    """
```

The method will also respect the "_Excluded" row attribute, if it exists, and omit those rows from
the calculation.

#### Permutation

Permute the order of the rows (or columns):

```python
def permute(self, ordering, axis):
    """
    Permute the dataset along the indicated axis.

    Args:
        ordering (list of int): 	The desired order along the axis
        
        axis (int):					The axis along which to permute

    Returns:
        Nothing.
    """
```

#### Feature selection

Select genes (rows) based on a CV vs mean fit:

```python
def feature_selection(self, n_genes):
    """
    Fits a noise model (CV vs mean)
    
    Args:
        n_genes (int):	number of genes to include

    Returns:
        Nothing.
    
    This method creates new row attributes _LogMean, _LogCV, _Noise (CV relative to predicted CV), _Excluded (1/0) 
    and now column attribute _TotalRNA
    """
```

#### Projection (t-SNE and PCA)

To perform a projection of the columns onto the 2D plane:

```python
def project_to_2d(self, perplexity = 20):
    """
    Project to 2D and create new column attributes _tSNE1, _tSNE2 and _PC1, _PC2.

    Args:
        perplexity (int): 	Perplexity to use for tSNE

    Returns:
        Nothing.

    This method first computes a PCA using scikit-learn IncrementalPCA (which doesn't load the whole
    dataset in RAM), then uses the top principal components to compute a tSNE projection. If row 
    attribute '_Excluded' exists, the projection will be based only on non-excluded genes.
    """
```

The method will always perform both t-SNE and PCA and store the resulting coordinates as new column
attributes _tSNE1, _tSNE2, _PCA1 and _PCA2.

PCA will be performed incrementally (i.e. without loading the entire matrix), and t-SNE will be 
performed on the top principal components. Thus the full dataset is never loaded into memory.
