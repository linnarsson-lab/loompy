# loompy

`.loom` is an efficient file format for very large omics datasets, 
consisting of a main matrix, optional additional layers, and a variable number of row and column 
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

## Getting started

Run the following in a Jupyter notebook:

```python
>>> import loompy
>>> ds = loompy.connect("cortex.loom")
>>> ds
```

This shows the upper-left 10x10 corner of the matrix along with its attributes:

<p>(18539, 1715)</p>
<table><tbody><tr><td>&nbsp;</td><td><strong>Cell_type</strong></td><td>eNb1</td><td>eNb1</td><td>eNb1</td><td>eNb1</td><td>eNb1</td><td>eNb1</td><td>eNb1</td><td>eNb1</td><td>eNb1</td><td>eNb1</td><td>...</td></tr><tr><td>&nbsp;</td><td><strong>Cell_ID</strong></td><td>1772122_301_C02</td><td>1772122_180_E05</td><td>1772122_300_H02</td><td>1772122_180_B09</td><td>1772122_180_G04</td><td>1772122_182_E09</td><td>1772122_302_C04</td><td>1772122_302_D11</td><td>1772122_180_C11</td><td>1772122_298_A07</td><td>...</td></tr><tr><td>&nbsp;</td><td><strong>Timepoint</strong></td><td>day_35</td><td>day_35</td><td>day_35</td><td>day_35</td><td>day_35</td><td>day_35</td><td>day_35</td><td>day_35</td><td>day_35</td><td>day_35</td><td>...</td></tr><tr><td><strong>Gene</strong></td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>...</td></tr><tr><td>DDX11L1</td><td>&nbsp;</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>...</td></tr><tr><td>WASH7P_p1</td><td>&nbsp;</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>...</td></tr><tr><td>LINC01002_loc4</td><td>&nbsp;</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>...</td></tr><tr><td>LOC100133331_loc1</td><td>&nbsp;</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>...</td></tr><tr><td>LOC100132287_loc2</td><td>&nbsp;</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>...</td></tr><tr><td>LOC101928626</td><td>&nbsp;</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>...</td></tr><tr><td>MIR6723</td><td>&nbsp;</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>...</td></tr><tr><td>LOC100133331_loc2</td><td>&nbsp;</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>...</td></tr><tr><td>LOC100288069_p1</td><td>&nbsp;</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>...</td></tr><tr><td>FAM87B</td><td>&nbsp;</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>...</td></tr><tr><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr></tbody></table>

In this example, the total size of the dataset is 18539 rows (genes) by 1715 columns (cells). There are 
three column attributes (`Cell_type`, `Cell_ID` and `Timepoint`) and one row attribute (`Gene`).

Next, try this:

```python
ds[ds.Gene == "Actb", :]
```

This returns an array of the expression values for *Actb*. Note the use of `ds.Gene == ...` to pick out rows (columns) that match some criterion, and the use of ds[..., ...] to select subsets of the data. In this example, `ds.Gene == "Actb"` is used as the row selector to pick out the single row corresponding to *Actb*, and `:` is used to select all columns. Hence, the expression returns the expression values for *Actb* in every cell in the dataset.

Refer to the API Documentation below to learn more about creating and manipulating loom files.

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
do not support journalling, so if something happens during a write, the 
**entire file can be lost**. Therefore, do not use loom files as 
your primary data storage. They are for working with data, not keeping 
it safe.

Loom files do support concurrent reads, but only from separate processes (not threads), 
and (we think) only from a single compute node. On a compute cluster, you may encounter
bugs if you try to read the same Loom file from different compute nodes concurrently.

Loom files are great for distribution of large datasets, which are then
used as read-only for analytical purposes.

### Efficient indexing and compression

The main matrix is stored in *chunked* format. That is, instead of being
stored by rows or by columns, it is stored as a sequence of little rectangles. 
As a consequence, both rows and columns (as well as submatrices) can be efficiently 
accessed. 

By default, chunks are compressed and decompressed for you automatically. This makes Loom
a space-efficient format for storing large sparse datasets. We have found that Loom 
often uses less space than standard sparse matrix formats.

### Matrix and attributes

Loom files represent a two-dimensional matrix with named row and column attributes. 
If the main matrix has N rows and M columns, then each row attribute has M elements (one per column) 
and vice versa.

By convention, the matrix represents measurements of some property (e.g. gene expression), genes are
stored in the rows and columns represent samples (e.g. cells).

The file can grow by the addition of attributes and columns (but not rows). For this reason, it's 
a good idea to put genes in the rows (since you will likely always work with the same  gene set).

### Data types

Loom supports a subset of numpy data types:

* The main matrix, and any additional layers, are always a two-dimensional arrays of any valid numpy type
* Attributes are one-dimensional arrays of either `float64` or `string`

Note that there is no integer attribute type. However, float64s are large enough to 
represent all integers up to and including 9,007,199,254,740,992 without loss.

## API Documentation

### Creating `.loom` files

Create from data:

```python

def create(filename: str, matrix: np.ndarray, row_attrs: Dict[str, np.ndarray], col_attrs: Dict[str, np.ndarray], file_attrs: Dict[str, str] = None, chunks: Tuple[int, int] = (64, 64), chunk_cache: int = 512, dtype: str = "float32", compression_opts: int = 2) -> LoomConnection:
	"""
	Create a new .loom file from the given data.

	Args:
		filename (str):         The filename (typically using a `.loom` file extension)
		matrix (numpy.ndarray): Two-dimensional (N-by-M) numpy ndarray of float values
		row_attrs (dict):       Row attributes, where keys are attribute names and values
								are numpy arrays (float or string) of length N
		col_attrs (dict):       Column attributes, where keys are attribute names and
								values are numpy arrays (float or string) of length M
		chunks (tuple):         The chunking of the matrix. Small chunks are slow
								when loading a large batch of rows/columns in sequence,
								but fast for single column/row retrieval.
								Defaults to (64,64).
		chunk_cache (int):      Sets the chunk cache used by the HDF5 format inside
								the loom file, in MB. If the cache is too small to
								contain all chunks of a row/column in memory, then
								sequential row/column access will be a lot slower.
								Defaults to 512.
		dtype (str):           Dtype of the matrix. Default float32 (uint16, float16 could be used)
		compression_opts (int): Strenght of the gzip compression. Default None.
	Returns:
		LoomConnection to created loom file.
	"""
```

Create by combining existing .loom files:

```python

def combine(files: List[str], output_file: str, key: str = None, file_attrs: Dict[str, str] = None) -> None:
	"""
	Combine two or more loom files and save as a new loom file

	Args:
		files (list of str):	the list of input files (full paths)
		output_file (str):		full path of the output loom file
		key (string):			Row attribute to use to verify row ordering
		file_attrs (dict):		file attributes (title, description, url, etc.)

	Returns:
		Nothing, but creates a new loom file combining the input files.

	The input files must (1) have exactly the same number of rows, (2) have
	exactly the same sets of row and column attributes.
	"""
```


Create from a 10X Genomics [cellranger](http://support.10xgenomics.com/single-cell/software/pipelines/latest/what-is-cell-ranger) output folder:


```
def create_from_cellranger(folder: str, loom_file: str, cell_id_prefix: str = '', sample_annotation: Dict[str, np.ndarray] = None, genome: str = 'mm10') -> LoomConnection:
	"""
	Create a .loom file from 10X Genomics cellranger output

	Args:
		folder (str):				path to the cellranger output folder (usually called `outs`)
		loom_file (str):			full path of the resulting loom file
		cell_id_prefix (str):		prefix to add to cell IDs (e.g. the sample id for this sample)
		sample_annotation (dict): 	dict of additional sample attributes
		genome (str):				genome build to load (e.g. 'mm10')

	Returns:
		Nothing, but creates loom_file
	"""
```
	

You can use the *sample_annotation* dictionary to add column (cell) annotations to all cells in the dataset. For example, this is useful to add a
 sample ID to each of several datasets before combining them into a single .loom file. 
 
### Connecting

Establish a connection to an existing `.loom` file:

```python
def connect(filename: str, mode: str = 'r+') -> LoomConnection:
	"""
	Establish a connection to a .loom file.

	Args:
		filename (str):     Name of the .loom file to open
		mode (str):         read/write mode, accepts 'r+' (read/write) or
							'r' (read-only), defaults to 'r+'

	Returns:
		A LoomConnection instance.
	"""
```

Example:

```python
ds = loompy.connect("filename.loom")
```

In the rest of the documentation below, `ds` is assumed to be an instance of `LoomConnection` obtained by connecting to a `.loom` file.

Note: there is usually no need to close the connection. The exception is if you need to write to the loom file from
two different processes (sequentially, not simultaneously). In that case, the first process needs to let go of the 
file by calling `close()` on the connection, before the second can start writing:

```python
ds.close()
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

### Global attributes

Global attributes are available as

```python
>>> ds.attrs["title"]
"The title of the dataset"

>>> ds.attrs["title"] = "New title"
>>> ds.attrs["title"]
"New title"
```

The following global attributes are standard:

  * `title`, a short title for the dataset
  * `description`, a longer description of the dataset
  * `url`, a link to a web page for the dataset
  * `doi`, a DOI for the paper where the dataset was published

(They are standard in the sense that you are encouraged to use `title` rather than `Title` or `TITLE` for a title, but they are
not guaranteed to exist, or required)

The following global attributes are reserved:

  * `schema`, a type annotation schema (JSON-formatted string)

DO NOT attempt to set reserved global attributes to a different value.


### Row and column attributes

Row and column attributes are accessed as dictionaries on `row_attrs` and `col_attrs`, respectively. For example:

```python
ds.row_attrs.keys()       # Return list of row attribute names
ds.col_attrs.keys()       # Return list of column attribute names
ds.row_attrs["GeneName"]  # Return a numpy array of gene names (assuming the attribute exists)
```

Note that these dictionaries are **read-only**. Any modifications will not be saved in the .loom file and will cause internal inconsistencies in the `LoomConnection` object. Use *set_attr()* (below) to add or modify attributes.

For convenience, attributes are also available directly on the `LoomConnection` object:

```python
ds.GeneName		# Equivalent to ds.row_attrs["GeneName"]
```

Using attributes in this way results in a very compact and readable syntax for selecting subarrays:

```python
>>> ds[ds.Gene == "Actb",:]
array([[  2.,   9.,   9., ...,   0.,  14.,   0.]], dtype=float32)

>>> ds[np.logical_or(ds.Gene == "Actb", ds.Gene == "Gapdh"),:]
array([[  2.,   9.,   9., ...,   0.,  14.,   0.],
       [  0.,   1.,   4., ...,   0.,  14.,   3.]], dtype=float32)

>>> ds[:, ds.CellID == "AAACATACATTCTC-1"]
array([[ 0.],
       [ 0.],
       [ 0.],
       ..., 
       [ 0.],
       [ 0.],
       [ 0.]], dtype=float32)
```

There are some limitations: 

* Custom attributes do not override existing `LoomConnection` attributes, such as method names. For example, if your .loom file has a row attribute `shape`, then `ds.shape` will not return that attribute, but will still return the shape of the main matrix. 
* Column attributes take precedence. For example, if you have both `ds.row_attrs["Name"]` and `ds.col_attrs["Name"]`, then `ds.Name` returns the column attribute, not the row attribute.

Note again, that you should not assign to these attributes, because your assignment will not be saved in the .loom file and will cause internal inconsistencies in the `LoomConnection` object. Use *set_attr()* (below) to add or modify attributes.


### Adding attributes and columns

You can add attributes and columns to an existing loom file. It is not possible to add rows or to delete 
attributes or any part of the matrix.

To add an attribute, which also saves it to the loom file:

```python
	def set_attr(self, name, values, axis = 0, dtype=None):
		"""
		Create or modify an attribute.

		Args:
			name (str): 			Name of the attribute
			values (numpy.ndarray):	Array of values of length equal to the axis length		
			axis (int):				Axis of the attribute (0 = rows, 1 = columns)
			dtype (str):			Type ("float64", "int", or "string")

		Returns:
			Nothing.

		This will overwrite any existing attribute of the same name.
		"""
```

**Note:** If you use an existing attribute name, the existing attribute will be overwritten. This is pefectly fine, 
and is the only way to change an attribute or its type.



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

You can also add the contents of another .loom file:

```python
	def add_loom(self, other_file: str, key: str = None, fill_values: Dict[str, np.ndarray] = None) -> None:
		"""
		Add the content of another loom file

		Args:
			other_file (str):	filename of the loom file to append
			fill_values (dict): default values to use for missing attributes (or None to drop missing attrs, or 'auto' to fill with sensible defaults)

		Returns:
			Nothing, but adds the loom file. Note that the other loom file must have exactly the same
			number of rows, and must have exactly the same column attributes.
			The all the contents including layers but ignores layers in `other_file` that are not already persent in self
		"""
```

The content of the other file is added as columns on the right of the current dataset. The rows must match for this to work. That is, the two files must have exactly the same rows (genes). If `key` is given, the rows may be out of order, and will be aligned based on the key attribute. Furthermore, the two datasets must have the same column attributes (but of coure can have different *values* for those attributes at each column). Missing attributes can be given default values using `fill_values`
.

### Operations

#### Map

You can map a function across all rows (all columns), while avoiding loading the entire
dataset into memory:

```python
	def map(self, f_list: List[Callable[[np.ndarray], int]], axis: int = 0, chunksize: int = 1000, selection: np.ndarray = None) -> List[np.ndarray]:
		"""
		Apply a function along an axis without loading the entire dataset in memory.

		Args:
			f (list of func):		Function(s) that takes a numpy ndarray as argument

			axis (int):		Axis along which to apply the function (0 = rows, 1 = columns)

			chunksize (int): Number of rows (columns) to load per chunk

			selection (array of bool): Columns (rows) to include

		Returns:
			numpy.ndarray result of function application

			If you supply a list of functions, the result will be a list of numpy arrays. This is more
			efficient than repeatedly calling map() one function at a time.
		"""
```

The function will receive an array (of floats) as its only argument, and should return a single float value.

Example:

```python
>>> import numpy as np
>>> ds.map([np.mean])[0]
# Returns an array of row means
np.array([1.23, 0.32, ...])   
```


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

#### Batch scan

For very large loom files, it's very useful to scan across the file (along either rows or columns) in *batches*,
to avoid loading the entire file in memory. This can be achieved using the `batch_scan` method:

```
	def batch_scan(self, cells: np.ndarray = None, genes: np.ndarray = None, axis: int = 0, batch_size: int = 1000) -> Iterable[Tuple[int, np.ndarray, np.ndarray]]:
		"""Performs a batch scan of the loom file

		Args
		----
		cells: np.ndarray
			the indexes [1,2,3,..,1000] of the cells to select
		genes: np.ndarray
			the indexes [1,2,3,..,1000] of the genes to select
		axis: int
			0:rows or 1:cols
		batch_size: int
			the chuncks returned at every element of the iterator

		Returns
		------
		Iterable that yields triplets
		(ix, indexes, vals)

		ix: int
			first position / how many rows/cols have been yielded alredy
		indexes: np.ndarray[int]
			the indexes with the same numbering of the input args cells / genes (i.e. np.arange(len(ds.shape[axis])))
			this is ix + selection
		vals: np.ndarray
			the matrix corresponding to the chunk
		"""
```


### Working with layers

Loom supports multiple layers. There is always a single main matrix, but optionally one or more additional layers
having the same number of rows and columns. Layers are accessed using the `layer` property on the `LoomConnection`. 

#### Create a layer

```
def set_layer(self, name: str, matrix: np.ndarray, chunks: Tuple[int, int] = (64, 64), chunk_cache: int = 512, dtype: str = "float32", compression_opts: int = 2) -> None:
```

#### Access a layer

The `layer` property returns a Layer object, which can be sliced to get the data:

```
ds.layer["layer"][10, :]
```

The default layer can be accessed directly:

```
ds[10, :]
```

It can also be accessed using the empty string:

```
ds.layer[""]
```

Layers can be loaded in memory as sparse matrices, efficiently:

```
LoomLayer.as_coo() -> sparse.coo_matrix:
LoomLayer.as_csr() -> sparse.csr_matrix:
LoomLayer.as_csc() -> sparse.csc_matrix:
```
