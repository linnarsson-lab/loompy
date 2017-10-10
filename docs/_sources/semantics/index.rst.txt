.. _semantics:

Understanding the semantics of loom files
-----------------------------------------

Connecting, not loading and saving
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Loom files are stored on disk and are never loaded entirely. They are
more like databases: you connect, retrieve some subset of the data,
maybe update some attributes.

When you connect, all attributes are read into memory for quick access,
but the main matrix remains on disk.

Reading and writing
~~~~~~~~~~~~~~~~~~~

Loom files are based on
`HDF5 <https://en.wikipedia.org/wiki/Hierarchical_Data_Format>`__, a
file format suitable for large multidimensional datasets. They are
designed to be mostly created once, then used as read-only. They **do
not** support writing and reading concurrently. They also do not support
journalling, so if something happens during a write, the **entire file
can be lost**. Therefore, do not use loom files as your primary data
storage. They are for working with data, not keeping it safe.

Loom files do support concurrent reads, but only from separate processes
(not threads), and (we think) only from a single compute node. On a
compute cluster, you may encounter bugs if you try to read the same Loom
file from different compute nodes concurrently.

Loom files are great for distribution of large datasets, which are then
used as read-only for analytical purposes.

Efficient indexing and compression
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The main matrix is stored in *chunked* format. That is, instead of being
stored by rows or by columns, it is stored as a sequence of little
rectangles. As a consequence, both rows and columns (as well as
submatrices) can be efficiently accessed.

By default, chunks are compressed and decompressed for you
automatically. This makes Loom a space-efficient format for storing
large sparse datasets. We have found that Loom often uses less space
than standard sparse matrix formats.

Matrix and attributes
~~~~~~~~~~~~~~~~~~~~~

Loom files represent a two-dimensional matrix with named row and column
attributes. If the main matrix has N rows and M columns, then each row
attribute has M elements (one per column) and vice versa.

By convention, the matrix represents measurements of some property (e.g.
gene expression), genes are stored in the rows and columns represent
samples (e.g. cells).

The file can grow by the addition of attributes and columns (but not
rows). For this reason, it's a good idea to put genes in the rows (since
you will likely always work with the same gene set).

Data types
~~~~~~~~~~

Loom supports a subset of numpy data types:

-  The main matrix, and any additional layers, are always a
   two-dimensional arrays of any valid numpy type
-  Attributes are one-dimensional arrays of either ``float64`` or
   ``string``

Note that there is no integer attribute type. However, float64s are
large enough to represent all integers up to and including
9,007,199,254,740,992 without loss.