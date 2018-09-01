new (function)
==============

Create a new empty Loom file and return it as a context manager. 

.. autofunction:: loompy.new


The ``new()`` function is especially useful when you want to combine data from multiple sources. You would 
first create an empty Loom file, then add data to it incrementally. For example, using ``new()`` and ``scan()``:

.. code:: python

  with loompy.new(out_file) as dsout:  # Create a new, empty, loom file
    for f in input_files:  # Loop over a list of input Loom files
      with loompy.connect(f) as ds:
        totals = ds.map([np.sum], axis=1)[0]  # Calculate the total molecule count for each cell
        cells = np.where(totals > 500)[0] # Select the cells that passed QC (totals > 500)
        for (ix, selection, view) in ds.scan(items=cells, axis=1, key="Accession"):
          dsout.add_columns(view.layers, col_attrs=view.ca, row_attrs=view.ra)

Note the use of ``key=="Accession"`` which makes ``scan()`` sort the rows of each file by the row 
attribute ``Accession``, ensuring that the resulting output is in a consistent order.
