MemoryLoomLayer (class)
=======================

A MemoryLoomLayer represents a layer of data residing in RAM only, and provides a numpy ndarray-like interface.
They are typically obtained by creating a :class:`.LoomView` on the LoomConnection.

.. highlight:: python
.. code-block:: python

    with loompy.connect("mydataset.loom") as ds:
        for (ix, selection, view) in ds.scan(axis=1):
            # Here, the matrix returned resides only in RAM and
            # each iteration gives a slab out of the full matrix
            print(view.layers["spliced"][0, :])


.. autoclass:: loompy.MemoryLoomLayer
    :special-members: __init__ 
    :members:
    :undoc-members:
