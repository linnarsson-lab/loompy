LoomValidator (class)
=====================

Validate the content of a Loom file against the file format specification and/or attrubute conventions.

.. highlight:: python
.. code-block:: python

    >>> lv = loompy.LoomValidator()
    >>> if not lv.validate(d+f, strictness="conventions"):
    >>>     for warn in lv.warnings:
    >>>         print("WARNING:", warn)
    >>>     for err in lv.errors:
    >>>         print("ERROR:", err)

    WARNING: Optional global attribute 'Description' is missing
    WARNING: Optional global attribute 'Journal' is missing
    WARNING: Optional global attribute 'Authors' is missing
    WARNING: Optional global attribute 'Title' is missing
    WARNING: Optional global attribute 'Year' is missing
    WARNING: Optional column attribute 'ClusterName' is missing
    WARNING: Optional column attribute 'Valid' is missing
    WARNING: Optional row attribute 'Valid' is missing
    WARNING: Optional row attribute 'Selected' is missing
    ERROR: Column attribute 'ClusterID' is missing
    ERROR: Column attribute 'CellID' cannot contain duplicate values
    ERROR: For help, see http://linnarssonlab.org/loompy/conventions/


.. autoclass:: loompy.LoomValidator
    :special-members: __init__ 
    :members:
    :undoc-members:
