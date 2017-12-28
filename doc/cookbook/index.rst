.. _cookbook:

Cookbook
========

In this section, we will show by example how to complete common tasks with idiomatic use of loompy.

Loading attributes from Pandas
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you have your metadata in a Pandas DataFrame, you can easily use it to create a Loom file. You will
need one DataFrame for the column metadata and one for the row metadata. Convert each DataFrame into a dictionary
of lists:

.. code:: python

  df_row_metadata = ... # A pandas DataFrame holding row metadata
  df_col_metadata = ... # A pandas DataFrame holding column metadata
  data = ... # A numpy ndarray holding the main dataset
  
  loompy.create(filename, data, df_row_metadata.todict("list"), df_col_metadata.todict("list"))


Combining data using scan()
^^^^^^^^^^^^^^^^^^^^^^^^^^^

We often want to scan through a number of input files (for example, raw
data files from multiple experiments), select a subset of the columns (e.g. cells passing QC)
and write them to a new file. This can be accomplished using the ``scan()`` method.

For example, let's select cells that have more than 500 detected UMIs in each of several files:

.. code:: python

  for f in input_files:
    with loompy.connect(f) as ds:
      totals = ds.map([np.sum], axis=1)[0]
      cells = np.where(totals > 500)[0] # Select the cells that passed QC (totals > 500)
      for (ix, selection, view) in ds.scan(items=cells, axis=1):
        loompy.create_append(out_file, view.layers, view.ra, view.ca)

Note that by using ``create_append`` we will first be creating the new file, then appending columns to it.

But what if the input files do not have their rows in the same order? ``scan()`` accepts a ``key`` argument 
to designate a primary key; each view is then returned sorted on the primary key on the *other axis*. 
For example, if you're scanning across columns, you should provide a row attribute as key, and each view will be sorted
on that attribute. 

Here's the same example, but this time we provide the ``key="Accession"`` argument to ensure that the input files
are sorted on the accession identifier along rows:

.. code:: python

  for f in input_files:
    with loompy.connect(f) as ds:
      totals = ds.map([np.sum], axis=1)[0]
      cells = np.where(totals > 500)[0] # Select the cells that passed QC (totals > 500)
      for (ix, selection, view) in ds.scan(items=cells, axis=1, key="Accession"):
        loompy.create_append(out_file, view.layers, view.ra, view.ca)

