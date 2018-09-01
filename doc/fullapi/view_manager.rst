ViewManager (class)
===================

The ViewManager is used by the LoomConnection.view attribute to return an in-memory :class:`.LoomView` of a slice through the dataset. Examples:

.. highlight:: python
.. code-block:: python

    with loompy.connect("mydataset.loom") as ds:
        myview = ds.view[:100, :100]  # Top-left 100x100 corner of the dataset
        print(myview.ra.Gene)  # Will print the 100 genes that are in the view

Views can also be created manually:

.. highlight:: python
.. code-block:: python

    with loompy.connect("mydataset.loom") as ds:
        vm = loompy.ViewManager(ds)  # Create a view manager (this does not load any data)
        myview = vm[:100, :100]  # This returns a view of the top-left 100x100 corner

.. autoclass:: loompy.ViewManager
    :members:
    :undoc-members:
