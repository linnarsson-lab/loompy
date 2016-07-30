# loom.py

`.loom` is an efficient file format for very large omics datasets, consisting of a main matrix and a variable number of row and column annotations. We use loom files to store single-cell gene expression data: the main matrix contains the actual expression values (one column per cell, one row per gene); row and column annotations contain metadata for genes and cells, such as `Name`, `Chromosome`, `Position` (for genes), and `Strain`, `Sex`, `Age` (for cells).

## Getting started

```python
import loom
ds = loom.connect("cortex.loom")
print ds.row_attrs.keys()
```
