LayerManager (class)
====================

The LayerManager manages the ``layers`` attribute on a :class:`.LoomConnection` and provides a 
dict-like interface. Examples:

.. highlight:: python
.. code-block:: python

    with loompy.connect("mydataset.loom") as ds:
        print(ds.layers.keys())
        print(f"There are {len(ds.layers)} layers")
        for name, layer in ds.layers.items():
            print(name, layer.shape, layer.dtype)

.. autoclass:: loompy.LayerManager
    :special-members: __init__ 
    :members:
    :undoc-members:
