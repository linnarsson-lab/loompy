.. _cookbook:

Cookbook
========

In this section, we will show by example how to complete common tasks with idiomatic use of loompy.

Working with a newly created file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Loompy 2 changes the behaviour of ``loompy.create``: it no longer returns a value. Thus in order to work with a newly created file you have to do this:

.. code:: python

  loompy.create("filename.loom", m, row_attrs, col_attrs)
  with loompy.connect("filename.loom") as ds:
      ....do something with ds
  # File closes automatically

The reason for the change is that we would often create files without closing the returned file handle, which causes issues especially in multi-process scenarios.

Note: if you simply want to create the file, and not access it, there is no need to use a ``with`` statement:

.. code:: python

  loompy.create("filename.loom", m, row_attrs, col_attrs)

This will leave the file closed.


Loading attributes from Pandas
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you have your metadata in a Pandas DataFrame, you can easily use it to create a Loom file. You will
need one DataFrame for the column metadata and one for the row metadata. Convert each DataFrame into a dictionary
of lists:

.. code:: python

  df_row_metadata = ... # A pandas DataFrame holding row metadata
  df_col_metadata = ... # A pandas DataFrame holding column metadata
  data = ... # A numpy ndarray holding the main dataset
  
  loompy.create(filename, data, df_row_metadata.to_dict("list"), df_col_metadata.to_dict("list"))


Combining data using scan() and new()
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We often want to scan through a number of input files (for example, raw
data files from multiple experiments), select a subset of the columns (e.g. cells passing QC)
and write them to a new file. This can be accomplished by creating an empty file using ``loompy.new()`` and
then filling it up using the ``scan()`` method.

For example, let's select cells that have more than 500 detected UMIs in each of several files:

.. code:: python

  with loompy.new(out_file) as dsout:  # Create a new, empty, loom file
    for f in input_files:
      with loompy.connect(f) as ds:
        totals = ds.map([np.sum], axis=1)[0]
        cells = np.where(totals > 500)[0] # Select the cells that passed QC (totals > 500)
        for (ix, selection, view) in ds.scan(items=cells, axis=1):
          dsout.add_columns(view.layers, col_attrs=view.ca, row_attrs=view.ra)

Note that by using ``new()`` we will first be creating the new file, then appending columns to it.

But what if the input files do not have their rows in the same order? ``scan()`` accepts a ``key`` argument 
to designate a primary key; each view is then returned sorted on the primary key on the *other axis*. 
For example, if you're scanning across columns, you should provide a row attribute as key, and each view will be sorted
on that attribute. 

Here's the same example, but this time we provide the ``key="Accession"`` argument to ensure that the input files
are sorted on the accession identifier along rows:

.. code:: python

  with loompy.new(out_file) as dsout:  # Create a new, empty, loom file
    for f in input_files:
      with loompy.connect(f) as ds:
        totals = ds.map([np.sum], axis=1)[0]
        cells = np.where(totals > 500)[0] # Select the cells that passed QC (totals > 500)
        for (ix, selection, view) in ds.scan(items=cells, axis=1, key="Accession"):
          dsout.add_columns(view.layers, col_attrs=view.ca, row_attrs=view.ra)


Fitting an incremental PCA
^^^^^^^^^^^^^^^^^^^^^^^^^^

Incremental algorithms are a powerful way of working with datasets that won't fit in RAM. For
example, we can use incremental PCA to learn a PCA transform by batch-wise partial fits:

.. code:: python

  from sklearn.decomposition import IncrementalPCA
  genes = (ds.ra.Selected == 1)
  pca = IncrementalPCA(n_components=50)
    for (ix, selection, view) in ds.scan(axis=1):
      self.pca.partial_fit(view[genes, :].transpose())
