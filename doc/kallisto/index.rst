.. _kallisto:

Create from fastq
=================

This section introduces the ``loompy`` command-line tool, and shows how it is used to create .loom files directly from fastq files.

The ``loompy`` command-line tool
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After installing loompy 3, a new command-line tool ``loompy`` will be available. To verify, launch a Terminal window and type ``loompy``,
and you should see the following output:

.. code:: 

  Usage: loompy [OPTIONS] COMMAND [ARGS]...

  Options:
    --show-message / --hide-message
    --verbosity [error|warning|info|debug]
    --help                          Show this message and exit.

  Commands:
    fromfq


Using the ``loompy`` command
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Install `kallisto <https://pachterlab.github.io/kallisto/>`_

The excellent kallisto tool performs ultra-fast pseudoalignment, which loompy uses to count reads (UMIs) on genes.

2. Download an index (or build your own)

We provide a pre-built index of the `human genome <https://docs.python-guide.org/starting/install3/linux/>`_. 

.. warning::
  This index is only suitable for use with 10x Chromium data (it makes strong assumptions about the read distribution).

Unzip the index to a directory, which will have the following content:

.. code:: 

  10xv1_whitelist.txt
  10xv2_whitelist.txt
  10xv3_whitelist.txt
  fragments2genes.txt
  gencode.v31.fragments.idx
  gencode.v31.metadata.tab
  manifest.json
  spliced_fragments.txt
  unspliced_fragments.txt


3. Download metadata and fastq files

TBD

4. Run the ``loompy`` command

Run the following command (replace ``index_folder`` with the path to the index you just downloaded):

.. code::

  loompy fromfq 10X05_1.loom 10X05_1 index_folder metadata.tab fastq-files

After about half an hour, you will have a ``10X05_1.loom`` file with separate ``spliced`` and ``unspliced`` layers (the main matrix will be
the sum of the two), and rich metadata for both genes, cells and the sample itself stored as attributes.

