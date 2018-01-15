

# loompy 2

⭐ Loompy v2.0 was released Dec. 24, 2017! ([what's new](https://github.com/linnarsson-lab/loompy/releases/tag/v2.0)?)

`.loom` is an efficient file format for very large omics datasets, 
consisting of a main matrix, optional additional layers, a variable number of row and column 
annotations. Loom also supports sparse graphs. We use loom files to store single-cell gene expression 
data: the main matrix contains the actual expression values (one 
column per cell, one row per gene); row and column annotations 
contain metadata for genes and cells, such as `Name`, `Chromosome`, 
`Position` (for genes), and `Strain`, `Sex`, `Age` (for cells).

![Illustration of Loom format structure](/doc/Loom-images.png)

Loom files (`.loom`) are created in the [HDF5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format) file format, which 
supports an internal collection of numerical multidimensional datasets.
HDF5 is supported by many computer languages, including Java, MATLAB, 
Mathematica, Python, R, and Julia. `.loom` files are accessible from 
any language that supports HDF5.

To get started, head over to [the documentation](http://linnarssonlab.org/loompy/)!

Loom, loompy, and the [loom-viewer](https://github.com/linnarsson-lab/loom-viewer) are being developed by members of the [Linnarsson Lab](http://linnarssonlab.org).

