LoomView (class)
================

A LoomView is in-memory Loom dataset, typically created by slicing the :meth:`.LoomConnection.view` attribute of 
a :class:`.LoomConnection`. It is designed to behave exactly like a LoomConnection, except that it's read-only. 
Examples:

.. highlight:: python
.. code-block:: python

    with loompy.connect("mydataset.loom") as ds:
        myview = ds.view[:100, :100]  # Top-left 100x100 corner of the dataset
        print(myview.ra.Gene)  # Will print the 100 genes that are in the view


.. autoclass:: loompy.LoomView
    :special-members: __init__ 
    :members:
    :undoc-members:
