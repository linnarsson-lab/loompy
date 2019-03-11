.. _format:

Loom file format specs
======================

Versions
--------

This specification defines the Loom file format version ``2.0.1``.


.. _formatinfo:

Introduction
------------

The Loom file format is designed to efficiently hold large omics
datasets. Typically, such data takes the form of a large matrix of
numbers, along with metadata for the rows and columns. For example,
single-cell RNA-seq data consists of expression measurements for all
genes (rows) in a large number of cells (columns), along with metadata
for genes (e.g. ``Chromosome``, ``Strand``, ``Location``, ``Name``), and
for cells (e.g. ``Species``, ``Sex``, ``Strain``, ``GFP positive``).

We designed Loom files to represent such datasets in a way that
treats rows and columns the same. You may want to cluster both genes and
cells, you may want to perform PCA on both of them, and filter based on
quality controls. SQL databases and other data storage solutions almost
always treat data as a *table*, not a matrix, and makes it very hard to
add arbitrary metadata to rows and columns. In contrast, Loom makes
this very easy.

Furthermore, current and future datasets can have tens of thousands of
rows (genes) and hundreds of thousands of columns (cells). We designed
Loom for efficient access to arbitrary rows and columns.

The annotated matrix format lends itself to very natural representation
of common analysis tasks. For example, the result of a clustering
algorithm can be stored simply as another attribute that gives the
cluster ID for each cell. Dimensionality reduction such as PCA or t-SNE,
similarly, can be stored as two attributes giving the projection
coordinates of each cell.

Finally, we recognize the importance of graph-based analyses of such
datasets. Loom supports graphs of both the rows (e.g. genes) and the
columns (e.g. cells), and multiple graphs can be stored each file.

.. _hd5concepts:

HDF5 concepts
-------------

The Loom format is based on
`HDF5 <https://en.wikipedia.org/wiki/Hierarchical_Data_Format>`__, a
standard for storing large numerical datasets. Quoting from h5py.org:

    An HDF5 file is a container for two kinds of objects: datasets,
    which are array-like collections of data, and groups, which are
    folder-like containers that hold datasets and other groups. The most
    fundamental thing to remember when using h5py is: *Groups work like
    dictionaries, and datasets work like NumPy arrays*.

A valid Loom file is simply an HDF5 file that contains specific
*groups* containing the main matrix as well as row and column
attributes. Because of this, Loom files can be created and read by
any language that supports HDF5, including `Python <http://h5py.org>`__,
`R <http://bioconductor.org/packages/release/bioc/html/rhdf5.html>`__,
`MATLAB <http://se.mathworks.com/help/matlab/low-level-functions.html>`__,
`Mathematica <https://reference.wolfram.com/language/ref/format/HDF5.html>`__,
`C <https://www.hdfgroup.org/HDF5/doc/index.html>`__,
`C++ <https://www.hdfgroup.org/HDF5/doc/cpplus_RM/>`__,
`Java <https://www.hdfgroup.org/products/java/>`__, and
`Ruby <https://rubygems.org/gems/hdf5/versions/0.3.5>`__.

.. _specifications:

Specification
-------------

A valid Loom file conforms to the following:

Main matrix and layers
^^^^^^^^^^^^^^^^^^^^^^

-  There MUST be a single `HDF5 dataset <hdf5 dataset append>`_ at ``/matrix``, of dimensions (N, M)
-  There can OPTIONALLY be a `HDF5 group <https://support.hdfgroup.org/HDF5/doc/H5.intro.html#Intro-OGroups>`_ ``/layers`` containing additional
   matrices (called "layers")
-  Each additional layer MUST have the same (N, M) shape
-  Each layer can have a different data type, compression, chunking etc.

Global attributes
^^^^^^^^^^^^^^^^^

-  There can OPTIONALLY be at least one `HDF5
   attribute <https://www.hdfgroup.org/HDF5/Tutor/crtatt.html>`__ on the
   root ``/`` group, which can be any valid scalar or multidimensional datatype and should be
   interpreted as attributes of the whole Loom file. 
-  There can OPTIONALLY be an `HDF5
   attribute <https://www.hdfgroup.org/HDF5/Tutor/crtatt.html>`__ on the
   root ``/`` group named ``LOOM_SPEC_VERSION``, a string value giving the
   loom file spec version that was followed in creating the file. See top of this
   document for the current version of the spec.


Row and column attributes
^^^^^^^^^^^^^^^^^^^^^^^^^

-  There MUST be a group ``/row_attrs``
-  There can OPTIONALLY be one or more datasets at ``/row_attrs/{name}``
   whose first dimension has length N
-  There MUST be a group ``/col_attrs``
-  There can OPTIONALLY be one or more datasets at ``/col_attrs/{name}``
   whose first dimension has length M

 
The datasets under ``/row_attrs`` should be semantically interpreted as
row attributes, with one value per row of the main matrix, and in the
same order. Therefore, all datasets under this group must be
arrays with exactly N elements, where N is the number of
rows in the main matrix.

The datasets under ``/col_attrs`` should be semantically interpreted as
column attributes, with one value per column of the main matrix, and in
the same order. Therefore, all datasets under this group must be
arrays with exactly M elements, where M is the number of
columns in the main matrix.

Row and column sparse graphs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  There MUST be a group ``/col_graphs``
-  There can OPTIONALLY be one or more groups at ``/col_graphs/{name}``
-  Under each ``/col_graphs/{name}`` group, there MUST be three one-dimensional datasets
   called ``a`` (integer), ``b`` (integer) and ``w`` (float). These should
   be interpreted as a sparse graph in `coordinate list <https://en.wikipedia.org/wiki/Sparse_matrix>`_ 
   format. The lengths of the three datasets MUST be equal, which defines the number 
   of edges in the graph. Note that the number of columns in the dataset defines 
   the vertices, so an unconnected vertex is one that has no entry in ``a`` or ``b``.
-  There MUST be a group ``/row_graphs``
-  There can OPTIONALLY be one or more groups at ``/row_graphs/{name}``
-  Under each ``/row_graphs/{name}`` group, there MUST be three one-dimensional datasets
   called ``a`` (integer), ``b`` (integer) and ``w`` (float). These should
   be interpreted as a sparse graph in `coordinate list <https://en.wikipedia.org/wiki/Sparse_matrix>`_
   format. The lengths of the three datasets MUST be equal, which defines the number 
   of edges in the graph. Note that the number of rows in the dataset defines 
   the vertices, so an unconnected vertex is one that has no entry in ``a`` or ``b``.
-  Vertex indexing is zero-based. When an entry in ``a`` or ``b`` is zero, this denotes the first column 
   in the matrix. If there are N columns, then vertices are numbered from 0 to N - 1. 

Datatypes
---------

The main matrix and additional layers MUST be two-dimensional arrays of one of these numeric types: ``int8``, ``int16``, ``int32``, ``int64``, ``uint8``, ``uint16``, ``uint32``, ``uint64``, ``float16``, ``float32`` and ``float64``. Each layer can have its own datatype.

Row and column attributes are multidimensional arrays whose first dimension matches the corresponding main matrix dimension. The elements MUST be of one of the numeric datatypes ``int8``, ``int16``, ``int32``, ``int64``, ``uint8``, ``uint16``, ``uint32``, ``uint64``, ``float16``, ``float32`` and ``float64`` or fixed-length ASCII strings.

Global attributes are scalars or multidimensional arrays of any shape, whose elements are any of the numeric datatypes ``int8``, ``int16``, ``int32``, ``int64``, ``uint8``, ``uint16``, ``uint32``, ``uint64``, ``float16``, ``float32`` and ``float64`` or fixed-length ASCII strings.

All strings in Loom files are stored as fixed-length null-padded 7-bit ASCII. ``h5dump`` should report something like this:

.. code::

  DATATYPE  H5T_STRING {
    STRSIZE 24;
    STRPAD H5T_STR_NULLPAD;
    CSET H5T_CSET_ASCII;
    CTYPE H5T_C_S1;
  }


Unicode characters outside 7-bit ASCII are stored using `XML entity encoding <https://en.wikipedia.org/wiki/List_of_XML_and_HTML_character_entity_references>`_, to ensure maximum compatibility. Strings SHOULD be decoded when read and encoded when written. A compatible implementation may choose to encode/decode or not, but MUST decode on reading if it encodes on writing.

.. _loomexample:

Example
-------

Here's an example of the structure of a valid Loom file:

+----------------------+-------------------------------+---------------------------------------------+
| Group                | Type                          | Description                                 |
+======================+===============================+=============================================+
| /matrix              | float32[N,M] or uint16[N,M]   | Main matrix of N rows and M columns         |
+----------------------+-------------------------------+---------------------------------------------+
| /layers/             | (subgroup)                    | Subgroup of additional matrix layers        |
+----------------------+-------------------------------+---------------------------------------------+
| /row\_attrs/         | (subgroup)                    | Subgroup of all row attributes              |
+----------------------+-------------------------------+---------------------------------------------+
| /row\_attrs/Name     | string[N]                     | Row attribute "Name" of type string         |
+----------------------+-------------------------------+---------------------------------------------+
| /col\_attrs/         | (subgroup)                    | Subgroup of all column attributes           |
+----------------------+-------------------------------+---------------------------------------------+
| /col\_attrs/CellID   | float64[M]                    | Column attribute "CellID" of type float64   |
+----------------------+-------------------------------+---------------------------------------------+
| /col\_graphs/        | (subgroup)                    | Subgroup of all column graphs               |
+----------------------+-------------------------------+---------------------------------------------+
| /col\_graphs/KNN     | (subgroup)                    | A column graph "KNN"                        |
+----------------------+-------------------------------+---------------------------------------------+
| /col\_graphs/KNN/a   | int32[E]                      | Vector of edge 'from' vertices              |
+----------------------+-------------------------------+---------------------------------------------+
| /col\_graphs/KNN/b   | int32[E]                      | Vector of edge 'to' vertices                |
+----------------------+-------------------------------+---------------------------------------------+
| /col\_graphs/KNN/w   | float32[E]                    | Vector of edge weights                      |
+----------------------+-------------------------------+---------------------------------------------+
| /row\_graphs/        | (subgroup)                    | Subgroup of all row graphs                  |
+----------------------+-------------------------------+---------------------------------------------+



