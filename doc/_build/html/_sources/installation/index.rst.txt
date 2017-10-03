.. _installation:

Installation
============

Using pip
---------

You can install the loompy package from PyPi with:

::

    pip install loompy

Tip: the package is updated often (don't worry, the format is stable;
even in the rare occasion that you need to update your code, your old
loom files won't break). To ensure that you have the latest version, do
this:

::

    pip install -U loompy

From source
-----------

Alternatively, you can install the latest version from source:

::

    git clone https://github.com/linnarsson-lab/loompy.git
    python setup.py install

If you just want to work with loom files within Python code, you should
be all set! We also made a web-app to make it easier to browse the data,
which you can install for local viewing, or set up for sharing loom
files from your own website. See the loom-viewer`<https://github.com/linnarsson-lab/loom-viewer/>`__repository for more information.


.. _gettingstarted:

Getting started
---------------

Run the following in a Jupyter notebook:

.. code:: python

    >>> import loompy
    >>> ds = loompy.connect("cortex.loom")
    >>> ds

This shows the upper-left 10x10 corner of the matrix along with its
attributes:

.. raw:: html

    <p>(18539, 1715)</p><table><tbody><tr><td>&nbsp;</td><td><strong>Cell_type</strong></td><td>eNb1</td><td>eNb1</td><td>eNb1</td><td>eNb1</td><td>eNb1</td><td>eNb1</td><td>eNb1</td><td>eNb1</td><td>eNb1</td><td>eNb1</td><td>...</td></tr><tr><td>&nbsp;</td><td><strong>Cell_ID</strong></td><td>1772122_301_C02</td><td>1772122_180_E05</td><td>1772122_300_H02</td><td>1772122_180_B09</td><td>1772122_180_G04</td><td>1772122_182_E09</td><td>1772122_302_C04</td><td>1772122_302_D11</td><td>1772122_180_C11</td><td>1772122_298_A07</td><td>...</td></tr><tr><td>&nbsp;</td><td><strong>Timepoint</strong></td><td>day_35</td><td>day_35</td><td>day_35</td><td>day_35</td><td>day_35</td><td>day_35</td><td>day_35</td><td>day_35</td><td>day_35</td><td>day_35</td><td>...</td></tr><tr><td><strong>Gene</strong></td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>&nbsp;</td><td>...</td></tr><tr><td>DDX11L1</td><td>&nbsp;</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>...</td></tr><tr><td>WASH7P_p1</td><td>&nbsp;</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>...</td></tr><tr><td>LINC01002_loc4</td><td>&nbsp;</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>...</td></tr><tr><td>LOC100133331_loc1</td><td>&nbsp;</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>...</td></tr><tr><td>LOC100132287_loc2</td><td>&nbsp;</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>...</td></tr><tr><td>LOC101928626</td><td>&nbsp;</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>...</td></tr><tr><td>MIR6723</td><td>&nbsp;</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>...</td></tr><tr><td>LOC100133331_loc2</td><td>&nbsp;</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>...</td></tr><tr><td>LOC100288069_p1</td><td>&nbsp;</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>...</td></tr><tr><td>FAM87B</td><td>&nbsp;</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>...</td></tr><tr><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr></tbody></table><br>


In this example, the total size of the dataset is 18539 rows (genes) by
1715 columns (cells). There are three column attributes (``Cell_type``,
``Cell_ID`` and ``Timepoint``) and one row attribute (``Gene``).

Next, try this:

.. code:: python

    ds[ds.Gene == "Actb", :]

This returns an array of the expression values for *Actb*. Note the use
of ``ds.Gene == ...`` to pick out rows (columns) that match some
criterion, and the use of ds[..., ...] to select subsets of the data. In
this example, ``ds.Gene == "Actb"`` is used as the row selector to pick
out the single row corresponding to *Actb*, and ``:`` is used to select
all columns. Hence, the expression returns the expression values for
*Actb* in every cell in the dataset.

Refer to the API Documentation below to learn more about creating and
manipulating loom files.