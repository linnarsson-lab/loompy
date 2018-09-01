LoomLayer (class)
=================

A LoomLayer represents a layer of data and provides a numpy ndarray-like interface.

.. highlight:: python
.. code-block:: python

    with loompy.connect("mydataset.loom") as ds:
        # Names of all the layers ("" is the main matrix)
        ds.layers.keys()

        # Upper left corner of the main matrix
        ds.layers[""][0,0]  
        
        # Shorthand to slice the main matrix
        ds[0,0]

        # Load the entire layer named "spliced"
        ds.layers["spliced"][:, :]

        # Shorthand access to the layer named "spliced"
        ds["spliced"][:, :]

        # Assign a row of data to the named layer
        ds["spliced"][0, :] = new_data

        # Create a new empty layer of `ìnt32`` elements
        ds["empty_layer"] = "int32"

Layers can also be accessed using numpy fancy indexing, e.g. with a vector of bools or a list of element indices. 

.. tip::
    Note that 
    using fancy indexing to slice more than ~1% of the rows (or columns) is inefficient. If you want to extract a larger subset 
    of rows or columns, you're better of using :meth:`.LoomConnection.scan`.

.. autoclass:: loompy.LoomLayer
    :members: 
    :undoc-members:
