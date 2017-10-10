.. loompy documentation master file, created by
   sphinx-quickstart on Tue Oct  3 00:11:17 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to loompy!
==================

``loompy`` is a python library to agilely work with the ``.loom`` file format.

``.loom`` is an efficient file format for very large omics datasets,
consisting of a main matrix, optional additional layers, and a variable
number of row and column annotations. We use loom files to store
single-cell gene expression data: the main matrix contains the actual
expression values (one column per cell, one row per gene); row and
column annotations contain metadata for genes and cells, such as
``Name``, ``Chromosome``, ``Position`` (for genes), and ``Strain``,
``Sex``, ``Age`` (for cells).

Loom files (``.loom``) are created in the
`HDF5 <https://en.wikipedia.org/wiki/Hierarchical_Data_Format>`__ file
format, which supports an internal collection of numerical
multidimensional datasets. HDF5 is supported by many computer languages,
including Java, MATLAB, Mathematica, Python, R, and Julia. ``.loom``
files are accessible from any language that supports HDF5.

Content
-------

* :ref:`genindex`

.. toctree::
   :hidden:

   self

.. toctree::

   installation/index

.. toctree::

   semantics/index

.. toctree::

   apiwalkthrough/index
   fullapi/index

.. toctree::
   format/index