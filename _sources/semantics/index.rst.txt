.. _semantics:

Understanding the semantics of loom files
-----------------------------------------

Connecting, not loading and saving
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Loom files are stored on disk and are never loaded entirely. They are
more like databases: you connect, retrieve some subset of the data,
maybe update some attributes.

When you connect, only some metadata is loaded,
but the main matrix, additional layers, graphs and
attribute values remain on disk. They are loaded only as 
needed when you access them through the connection object.

The canonical way to connect to a Loom file is:

.. code:: python

    with loompy.connect("filename.loom") as ds:
        # do something with the connection object ds

This ensures that the file is automatically closed when no longer needed.

For interactive usage, it may be more convenient to connect directly:

.. code:: python

    ds = loompy.connect("filename.loom")
    # do something with the connection object ds
    ds.close()

In that case, you need to close the connection explicitly.


Reading and writing
~~~~~~~~~~~~~~~~~~~

Loom files are based on
`HDF5 <https://en.wikipedia.org/wiki/Hierarchical_Data_Format>`__, a
file format suitable for large multidimensional datasets. Loom files 
are great for distribution of large datasets, and for use
in computational pipelines. They **do
not** support writing and reading concurrently. They also do not support
journalling, so if something happens during a write, the **entire file
can be lost**. Therefore, do not use loom files as your primary data
storage. They are for working with data, not keeping it safe.

Loom files do support concurrent reads, but only from separate processes
(not threads), and (we think) only from a single compute node. On a
compute cluster, you may encounter bugs if you try to read the same Loom
file from different compute nodes concurrently.

Efficient indexing and compression
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The main matrix is stored in *chunked* format. That is, instead of being
stored by rows or by columns, it is stored as a sequence of little
rectangles. As a consequence, both rows and columns (as well as
submatrices) can be efficiently accessed.

By default, chunks are compressed and decompressed for you
automatically. This makes Loom a space-efficient format for storing
large sparse datasets. We have found that Loom often uses the same
amount of space as standard sparse matrix formats.

Matrix and attributes
~~~~~~~~~~~~~~~~~~~~~

Loom files represent a two-dimensional matrix with named row and column
attributes. If the main matrix has N rows and M columns, then each column
attribute has M elements (one per column) and vice versa.

By convention, the matrix represents measurements of some property (e.g.
gene expression), genes are stored in the rows and columns represent
samples (e.g. cells).

The file can grow by the addition of attributes and columns (but not
rows). For this reason, it's a good idea to put genes in the rows (since
you will likely always work with the same gene set).

Data types
~~~~~~~~~~

Loom supports the following data types:

-  The main matrix, and any additional layers, are
   two-dimensional arrays of numbers, which can be any 8, 16, 32 or
   64-bit integer (signed or unsigned), or 16, 32 or 64-bit float.
-  Attributes are N-dimensional arrays, and the elements can be any 8, 16, 32 or
   64-bit integer (signed or unsigned), or 16, 32 or 64-bit float,
   or string. The size of the first dimension must match the
   corresponding matrix dimension. 

See the file format specification for detailed information on data types,
including a discussion on string handling, ASCII and Unicode support.
