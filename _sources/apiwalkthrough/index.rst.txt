.. _apiwalkthrough:

API Walkthrough
===============

.. _loomcreate:

Creating and connecting
-----------------------

Creating Loom files
~~~~~~~~~~~~~~~~~~~

To create a loom file from data, you need to supply a main matrix (numpy ndarray or scipy sparse matrix) and two dictionaries of row and column attributes (with attribute names as keys, and numpy ndarrays as values). If the main matrix is N×M, then the row attributes must have N elements, and the column attributes must have M elements. 

For example, the following creates a loom file with a 100x100 main matrix, one row attribute and one column attribute:

.. code:: python

  import numpy as np
  import loompy
  filename = "test.loom"
  matrix = np.arange(10000).reshape(100,100)
  row_attrs = { "SomeRowAttr": np.arange(100) }
  col_attrs = { "SomeColAttr": np.arange(100) }
  loompy.create(filename, matrix, row_attrs, col_attrs)

:func:`loompy.create` accepts numpy dense matrices (:class:`numpy.ndarray`) as well as scipy sparse matrices (:class:`scipy.sparse.coo_matrix`, 
:class:`scipy.sparse.csc_matrix`, or :class:`scipy.sparse.csr_matrix`). For example:

.. code:: python

  import numpy as np
  import loompy
  import scipy.sparse as sparse
  filename = "test.loom"
  matrix = sparse.coo_matrix((100, 100))
  row_attrs = { "SomeRowAttr": np.arange(100) }
  col_attrs = { "SomeColAttr": np.arange(100) }
  loompy.create(filename, matrix, row_attrs, col_attrs)

Note that :func:`loompy.create` does not return anything. To work with the newly created file, you must :func:`loompy.connect` to it.

You can also create an empty file using :func:`loompy.new`, which returns a connection to the newly created file. The file can then be populated with data. 
This is especially useful when you're building a dataset incrementally, e.g. by pooling subsets of other datasets:

.. code:: python

  with loompy.new("outfile.loom") as dsout:
      for sample in samples:
          with loompy.connect(sample) as dsin:
              logging.info(f"Appending {sample}.")
              dsout.add_columns(ds.layers, col_attrs=dsin.col_attrs, row_attrs=dsin.row_attrs)

You can also create a file by combining existing loom files (:func:`loompy.combine`). The files will be concatenated along the column
axis, and therefore must have the same number of rows. If the rows are potentially not in the same order, 
you can supply a ``key`` argument; the row attribute corresponding to the key will be used to sort the files. 
For example, the following code will combine files and use the "Accession" row attribute as the key: 

.. code:: python

  loompy.combine(files, output_filename, key="Accession")

You can import a 10X Genomics
`cellranger <http://support.10xgenomics.com/single-cell/software/pipelines/latest/what-is-cell-ranger>`__
output folder using :func:`loompy.create_from_cellranger`:

.. code:: python

  loompy.create_from_cellranger(folder, output_filename)


Connecting to Loom files
~~~~~~~~~~~~~~~~~~~~~~~~

In order to work with a loom file, you must first :func:`loompy.connect` to it. This does not load the data
or attributes, so is very quick regardless of the size of the file. It's more like connecting to a 
database than reading a file. Loom supports Python context management, so normally you should use 
a ``with`` statement to take care of the connection:

.. code:: python

  with loompy.connect("filename.loom") as ds:
    # do something with ds

The connection will be automatically closed at the end of the ``with`` block.

Sometimes, especially in interactive use in a Jupyter notebook, you may want
to just open the file and keep the connection around:

.. code:: python

  ds = loompy.connect("filename.loom")

In that case, you should close the file when you are done:

.. code:: python

    ds.close()

In most cases, forgetting to close the file will do no harm, but may (for example)
prevent concurrent processes from accessing the file, and will leak file handles.

In the rest of the documentation below, ``ds`` is assumed to be an
instance of :class:`.LoomConnection` obtained by connecting to a Loom
file.

.. _loommanipulate:

Manipulate data
---------------

Shape, indexing and slicing
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :attr:`.LoomConnection.shape` attribute returns the row and column count as a tuple:

.. code:: python

    >>> ds.shape
    (100, 2345)

The data stored in the main matrix can be retrieved by indexing and
slicing. The following are supported:

-  Indices: anything that can be converted to a Python long
-  Slices (i.e. ``:`` or ``0:10``)
-  Lists of the rows/columns you want (i.e. ``[0, 34, 576]``)
-  Mask arrays (i.e. numpy array of bool indicating the rows/columns you
   want)

Lists and mask arrays are supported along one dimension at a time only. Since 
the main matrix is two-dimensional, two arguments are always needed. Examples:

.. code:: python


    ds[:, :]          # Return the entire matrix
    ds[0:10, 0:10]    # Return the 10x10 submatrix starting at row and column zero 
    ds[99, :]         # Return the 100th row 
    ds[:, 99]         # Return the 100th column
    ds[[0,3,5], :]    # Return rows with index 0, 3 and 5
    ds[:, bool_array] # Return columns where bool_array elements are True

Note that performance will be poor if you select many individual rows (columns) out
of a large matrix. For example, in a dataset with shape (27998, 160796), loading ten 
randomly chosen individual full columns took 914 ms, 
whereas loading 1000 columns took 1 minute and 6 seconds, and loadingh 5000 columns took 13 minutes.
This slowdown is caused by a `performance bug <https://github.com/h5py/h5py/issues/293>`_ 
in h5py.

If the whole dataset fits in RAM, loading it in full and then selecting the row/columns you want
will be faster. If it doesn't, consider using the :func:`.LoomConnection.scan` method (see below), which in this example took
1 minute and 12 seconds regardless of how many columns were selected. As a rule of thumb,
:func:`.LoomConnection.scan` will be faster whenever you are loading more than about 1% of the rows
or columns (randomly selected). 

  
Sparse data
~~~~~~~~~~~

On disk, every layer is stored chunked and block-compressed, for efficient storage and access along both axes.

The main matrix and additional layers can be assigned from dense or sparse matrices.

You can load the main matrix or any layer as sparse:

.. code:: python

  ds.layers["exons"].sparse()  # Returns a scipy.sparse.coo_matrix
  ds.layers["unspliced"].sparse(rows, cols)  # Returns only the indicated rows and columns (ndarrays of integers or bools)

You can assign layers from sparse matrices:

.. code:: python

  ds.layers["exons"] = my_sparse_matrix


Modifying layers
~~~~~~~~~~~~~~~~

You can modify the data in any layer by assigning to a slice. For example:


.. code:: python

    ds[:, :] = newdata         # Assign a full matrix
    ds[3, 500] = 31            # Set the element at (3, 500) to the value 31
    ds[99, :] = rowdata        # Assign new values to row with index 99
    ds[:, 99] = coldata        # Assign new values to column with index 99


Global attributes
~~~~~~~~~~~~~~~~~

Global attributes are available at ``ds.attrs`` and can be accessed by name or
as a dictionary. You create new attributes by assignment, and delete them
using the ``del`` statement:

.. code:: python

    >>> ds.attrs.title
    "The title of the dataset"

    >>> ds.attrs.title = "New title"
    >>> ds.attrs["title"]
    "New title"

    >>> del ds.attrs.title

You can list the attributes and loop over them as you would with a dictionary:

.. code:: python

  >>> ds.attrs.keys()
  ["title", "description"]

  >>> for key, value in ds.attrs.items():
  >>>   print(f"{key} = {value}")
  title = New title
  description = Fancy dataset

Global attributes can be scalars, or multidimensional arrays of any shape, and 
the elements can be integers, floats or strings. See below for the exact types allowed.

Row and column attributes
~~~~~~~~~~~~~~~~~~~~~~~~~

Row and column attributes are accessed at ``ds.ra``
and ``ds.ca``, respectively, and support the same interface as global 
attributes. For example:

.. code:: python

  ds.ra.keys()       # Return list of row attribute names
  ds.ca.keys()       # Return list of column attribute names
  ds.ra.Gene = ...   # Create or replace the Gene attribute
  a = ds.ra.Gene     # Assign the array of gene names (assuming the attribute exists)
  del ds.ra.Gene     # Delete the Gene row attribute

Attributes can also be accessed by indexing:

.. code:: python

  a = ds.ra["Gene"]     # Assign the array of gene names (assuming the attribute exists)
  del ds.ra["Gene"]     # Delete the Gene row attribute

You can pick out multiple attributes into a single numpy array, as long as they have the same type:

.. code:: python

  a = ds.ra["Gene", "Attribute"]     # Returns a 2D array of shape (n_genes, 2)
  b = ds.ca["PCA1", "PCA2"]          # Returns a 2D array of shape (n_cells, 2)

Note that when you ask for multiple attributes, missing attributes are silently ignored. This can be
exploited to access attributes that may have different names:

.. code:: python

  a = ds.ra["Gene", "GeneName"]     # Return one or the other (if only one exists)
  b = ds.ca["TSNE", "PCA", "UMAP"]  # Return the one that exists (if only one exists)

(of course, if two or more attributes exists, they will be stacked as above)


Attributes can be any of the following:

* One-dimensional arrays of integers, floats or strings. The number of elements in the array must match the corresponding matrix dimension.

* Multidimensional arrays of any of the same element types. The length along the first dimension of a row attribute must equal the number of rows in the main matrix (and vice versa for column attributes). Remaining dimensions can be any size. 

For example, if the main matrix has M columns, the result of a dimensionality reduction
(for example, a PCA) to 20 dimensions could be stored as a column attribute with shape (M, 20).

You can assign attributes using almost any array or list-like type, but attributes will 
always return numpy array (``np.ndarray``). 

Using attributes as masks for indexing the main matrix results in a very compact and readable
syntax for selecting subarrays:

.. code:: python

    >>> ds[ds.ra.Gene == "Actb", :]
    array([[  2.,   9.,   9., ...,   0.,  14.,   0.]], dtype=float32)

    >>> ds[(ds.ra.Gene == "Actb") | (ds.ra.Gene == "Gapdh"), :]
    array([[  2.,   9.,   9., ...,   0.,  14.,   0.],
           [  0.,   1.,   4., ...,   0.,  14.,   3.]], dtype=float32)

    >>> ds[:, ds.ca.CellID == "AAACATACATTCTC-1"]
    array([[ 0.],
           [ 0.],
           [ 0.],
           ..., 
           [ 0.],
           [ 0.],
           [ 0.]], dtype=float32)

Note that numpy logical functions overload the bitwise, not the boolean operators. Use ``|`` 
for 'or', ``&`` for 'and' and ``~`` for 'not'. You also must place parentheses around the comparison 
expressions to ensure proper operator precedence. For example:

.. code:: python

  (a == b) & (a > c) | ~(c <= b)


Modifying attributes
~~~~~~~~~~~~~~~~~~~~

Unlike layers, attributes are always only read and written in their entirety. Thus, assigning to a slice 
does not modify the attribute on disk. To write new values for an attribute, you must assign a 
full list or ndarray to the attribute:

.. code:: python

  with loompy.connect("filename.loom") as ds:
    ds.ca.ClusterNames = values  # where values is a list or ndarray with one element per column
    # This does not change the attribute on disk:
    ds.ca.ClusterNames[10] = "banana"



Adding columns
~~~~~~~~~~~~~~

You can add columns to an existing loom file. It's not possible to add rows or to 
delete any part of the matrix.

.. code:: python

 ds.add_columns(submatrix, col_attrs)

You need to provide a submatrix corresponding to the new columns, as well as
a dictionary of column attributes with values for all the new columns. 

Note that if you are adding columns to an empty file, you must also provide row attributes:

.. code:: python

 ds.add_columns(submatrix, col_attrs, row_attrs={"Gene": genes})


You can also add the contents of another .loom file:

.. code:: python

  ds.add_loom(other_file, key="Gene")

The content of the other file is added as columns on the right of the
current dataset. The rows must match for this to work. That is, the two
files must have exactly the same number of rows. If ``key`` is given, the
rows will be ordered based on the key attribute. Furthermore, the two 
datasets must have the same column
attributes (but of course can have different *values* for those
attributes at each column). Missing attributes can be given default
values using the ``fill_values`` argument. If the files contain any global attribute
with conflicting values, you can automatically convert such attributes into column attributes
by passing ``convert_attrs=True`` to the method.



.. _loomlayers:

Layers
~~~~~~~~~~~~~~~~~~~

Loom supports multiple layers. There is always a single main matrix, but
optionally one or more additional layers having the same number of rows
and columns. Layers are accessed using the ``layers`` property on the
``LoomConnection`` object.

Layers support the same pythonic API as attributes:

.. code:: python

  ds.layers.keys()            # Return list of layers
  ds.layers["unspliced"]      # Return the layer named "unspliced"
  ds.layers["spliced"] = ...  # Create or replace the "spliced" layer
  a = ds.layers["spliced"][:, 10] # Assign the 10th column of layer "spliced" to the variable a
  del ds.layers["spliced"]     # Delete the "spliced" layer

The main matrix is availabe as a layer named "" (the empty string). It cannot be deleted but
otherwise supports the same operations as any other layer.

As a convenience, layers are also available directly on the connection object. The above
expressions are equivalent to the following:

.. code:: python

  ds["unspliced"]      # Return the layer named "unspliced"
  ds["spliced"] = ...  # Create or replace the "spliced" layer
  a = ds["spliced"][:, 10] # Assign the 10th column of layer "spliced" to the variable a
  del ds["spliced"]     # Delete the "spliced" layer

Sometimes you may need to create an empty layer (all zeros), to be filled later. Empty layers
are created by assigning a type to a layer name. For example:

.. code:: python

  ds["empty_floats"] = "float32"
  ds["empty_ints"] = "int64"


.. _loomoperations:

Graphs
~~~~~~

Loom supports sparse graphs with either the rows or the columns as nodes. For example,
a sparse graph of cells (stored in the columns) could represent a K nearest-neighbors 
graph of the cells. In that case, the cells are the nodes (so there are M nodes in the
graph if there are M columns in the main matrix), which are connected by an arbitrary
number of edges. The graph could be considered directed or undirected, and can have float-valued
weights on the edges. Loom even supports multigraphs (permitting multiple edges between pairs of nodes).
Graphs are stored as arrays of edges and the associated edge weights.

Row and column graphs are accessed at ``ds.row_graphs`` and ``ds.col_graphs``, respectively, 
and support the same interface as attributes. For example:

.. code:: python

  ds.row_graphs.keys()      # Return list of row graphs
  ds.col_graphs.KNN = ...   # Create or replace the column-oriented graph KNN
  a = ds.col_graphs.KNN     # Assign the KNN column graph to variable a
  del ds.col_graphs.KNN     # Delete the KNN graph

Graphs are returned as ``scipy.sparse.coo_matrix``, and can be created/assigned from any
scipy sparse format as well as from a numpy dense matrix or ndarray. In each case, the matrix
represents the adjacency matrix of the graph.

Views
~~~~~

Loompy views are in-memory views of a slice through the underlying loom file. Views can be created
explicitly by slicing:

.. code:: python

  ds.view[:, 10:20]

This will create a view, fully loaded in memory, containing all the rows of the underlying loom file,
but only columns 10 through 19 (zero-based). You can use fancy indexing including slices, arrays of integers
(to pick out specific rows/columns) and boolean arrays.

The power of the view is that it slices through *everything*: the main matrix, every layer, every attribute, 
and every graph. This hides a lot of messy and error-prone code,
and makes it easy to extract relevant subsets of a loom file.

The most common use of a ``view`` is in scanning through a file (see ``scan()`` below).

Operations
~~~~~~~~~~

Map
^^^

You can map one or more functions across all rows (all columns), while avoiding
loading the entire dataset into memory:

.. code:: python

  ds.map([np.mean, np.std], axis=1)

The functions will receive an array (of floats or integers) as their only argument, and
should return a single float or integer value. Internally, ``map()`` uses ``scan()`` to
loop across the file.

Note that you must always provide a list of functions, even if it has only one element, and
that the result is a list of vectors, one per function that was supplied. Hence the correct 
way to map a single function across the matrix is:

.. code:: python

  (means,) = ds.map([np.mean], axis=1)


Permutation
^^^^^^^^^^^

Permute the order of the rows or columns:

.. code:: python

  ordering = np.random.permutation(np.arange(ds.shape[1]))
  ds.permute(ordering, axis=1)

This permutes the order of rows or columns in the file, without loading
the entire file in RAM. The ``ordering`` argument should be a numpy array
of ds.shape[axis] elements, in the desired order.

Scan
^^^^

For very large loom files, it's very useful to scan across the file
(along either rows or columns) in *batches*, to avoid loading the entire
file in memory. This can be achieved using the ``scan()`` method:

.. code:: python

  for (ix, selection, view) in ds.scan(axis=1):
    # do something with each view

Inside the loop, you get access to the current ``view`` into the file. It has all the 
attributes, graphs and data of the original loom file, but only for the columns included 
in ``selection`` (or rows, if axis=0). 

In essence, you get a succession of slices through the loom file, corresponding to 
bands of columns (rows). The ``ix`` variable tells you the starting column of the band, whereas
the ``selection`` gives you the list of columns contained in the current view.

You can also scan across a selected subset of the columns or rows. For example:

.. code:: python

  cells = # List of columns you want to see
  for (ix, selection, view) in ds.scan(items=cells, axis=1):
    # do something with each view

This works exactly the same, except that each ``selection`` and ``view`` now include only 
the columns you asked for.


