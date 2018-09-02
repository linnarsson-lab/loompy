.. loompy documentation master file, created by
   sphinx-quickstart on Tue Oct  3 00:11:17 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Loompy documentation
====================

Loom is an efficient file format for very large omics datasets,
consisting of a main matrix, optional additional layers, a variable
number of row and column annotations, and sparse graph objects. We use loom files to store
single-cell gene expression data: the main matrix contains the actual
expression values (one column per cell, one row per gene); row and
column annotations contain metadata for genes and cells, such as
``Name``, ``Chromosome``, ``Position`` (for genes), and ``Strain``,
``Sex``, ``Age`` (for cells). Graph objects are used to store nearest-neighbor
graphs used for graph-based clustering.

The Loom logo illustrates how all the parts fit together:

|

.. image:: Loom_components.png
   :width: 976px
   :height: 345px
   :scale: 50 %
   :align: center

|

Loom files (``.loom``) are created in the
`HDF5 <https://en.wikipedia.org/wiki/Hierarchical_Data_Format>`__ file
format, which supports an internal collection of numerical
multidimensional datasets. HDF5 is supported by many computer languages,
including `Python <http://h5py.org>`__,
`R <http://bioconductor.org/packages/release/bioc/html/rhdf5.html>`__,
`MATLAB <http://se.mathworks.com/help/matlab/low-level-functions.html>`__,
`Mathematica <https://reference.wolfram.com/language/ref/format/HDF5.html>`__,
`C <https://www.hdfgroup.org/HDF5/doc/index.html>`__,
`C++ <https://www.hdfgroup.org/HDF5/doc/cpplus_RM/>`__,
`Java <https://www.hdfgroup.org/products/java/>`__, and
`Ruby <https://rubygems.org/gems/hdf5/versions/0.3.5>`__.


.. toctree::
   :hidden:

   self
   installation/index
   semantics/index
   apiwalkthrough/index
   cookbook/index
   conventions/index
   format/index
   fullapi/index
