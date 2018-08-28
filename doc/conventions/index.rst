.. _conventions:

Conventions
===========

In order to maximize interoperability of Loom files between analysis pipelines, 
we suggest adhering to the following conventions. 

**Note:** This document is work in progress and subject to change! You can follow
the `discussion here <https://github.com/linnarsson-lab/loompy/issues/19>`_.


Single analysis per file
------------------------

Each Loom file stores a single analysis. If you want to try two different ways of clustering,
you store the results separately.

This convention simplifies things a lot, because we only need to keep track of one set of cluster
labels (for example). It also means that we don't need to store the relationship between different
attributes (e.g. that *this clustering* was done using *that PCA*). Such relationships must be 
stored external to the file.


Orientation
-----------

* Columns represent cells or aggregates of cells
* Rows represent genes

Loom files can grow along the column axis, but not the row axis, so this makes sense. 


Attribute naming conventions
----------------------------

We propose that algorithms always accept an argument specifying the name of each attribute it will use, with 
the defaults set as listed below. For example, an algorithm that needs a unique gene accession string would 
take an argument ``accession_attr="Accession"``, and a clustering algorithm that generates cluster labels would
take an argument ``cluster_id_attr="ClusterID"``. In this way, Loom files that conform to the conventions below
would work without fuss, but files that don't could still be made to work by supplying the non-standard attribute names.


Column attributes
^^^^^^^^^^^^^^^^^

``CellID`` a string label unique to each cell (preferrably globally unique across all datasets)

``Valid`` integers 1 or 0, indicating cells that are considered valid after some QC step

``ClusterID`` an integer label with values in [0, n_clusters].

``ClusterName`` a string label with arbitrary values representing cluster names. If both ``ClusterID`` and
``ClusterName`` are present, they should correspond 1:1.

``Outliers`` an integer label 1 or 0, indicating cells that are outliers relative to the clusters.

Any attribute that is an M-by-Y matrix (where M is the number of columns and Y is the dimensionality of the embedding) can be used
to store e.g. a PCA or t-SNE dimensionality reduction.

Row attributes
^^^^^^^^^^^^^^

``Gene`` a string human-readable gene name, not necessarily unique

``Accession`` a string, unique within the file (e.g. an ENSEMBL accession etc.)

``Selected`` integers 1 or 0, indicating genes that were selected by some previous step.

``Valid`` integers 1 or 0, indicating genes that are considered valid after some QC step




