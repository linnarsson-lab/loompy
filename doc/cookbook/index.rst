
Cookbook
--------

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


