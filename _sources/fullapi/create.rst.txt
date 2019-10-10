create (function)
=================

Create a new Loom file. 

.. autofunction:: loompy.create


In order to work with a newly created file you have to connect to it:

.. code:: python

  loompy.create("filename.loom", m, row_attrs, col_attrs)
  with loompy.connect("filename.loom") as ds:
      ....do something with ds
  # File closes automatically

Note: if you simply want to create the file, and not access it, there is no need to use a ``with`` statement:

.. code:: python

  loompy.create("filename.loom", m, row_attrs, col_attrs)

This will leave the file closed.