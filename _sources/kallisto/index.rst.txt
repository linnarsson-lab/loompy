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


Using the ``loompy fromfq`` command
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Install `kallisto <https://pachterlab.github.io/kallisto/>`_
---------------------------------------------------------------

The excellent kallisto tool performs ultra-fast pseudoalignment, which loompy uses to count reads (UMIs) on genes.

2. Download an index (or build your own)
----------------------------------------

We provide a pre-built index of the `human genome <https://storage.googleapis.com/linnarsson-lab-www-blobs/human_GRCh38_gencode.v31.tar.gz>`_. 

.. warning::
  This index is really only suitable for use with 10x Chromium data (it makes strong assumptions about the read distribution).

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

To build your own index, see the `build instructions <https://github.com/linnarsson-lab/loompy/blob/master/notebooks/build_index.ipynb>`_.

3. Download metadata and fastq files
-------------------------------------

You need to provide the input fastq files, and metadata for your samples. For this tutorial, please download the 1k PBMC v3 reference dataset 
from 10x Genomics: `pbmc_1k_v3_fastqs.tar <http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar>`_.

Metadata can be provided either in the form of a tab-delimited file with a single header row, or as a sqlite3 database with a table named ``sample``.

One metadata column must be named ``name`` and contain sample names (sample IDs). One column must be named ``technology`` and contain the technology name (currently one of ``10xv1``,
``10xv2`` or ``10xv3``). The technology name will be used to locate the barcode whitelist file (e.g. ``10xv3_whitelist.txt``). Finally, one
column must be named ``targetnumcells`` and contain the target (or expected) cell number in the samples. If absent, the default will be 5000.

The metadata file (or database) will typically contain metadata for all the samples you are working with, one row per sample. You can add any number of
additional columns of metadata, which will all be copied to the newly generated loom file as global attributes.

For this tutorial, create a tab-delimited file named ``metadata.tab`` with the following content:

.. code::

  name    technology  targetnumcells
  1kPBMC  10xv3       1000

(If you copy-paste the above, make sure the fields are tab delimited)


4. Run the ``loompy fromfq`` command
------------------------------------

Run the following command (replace ``human_GRCh38_gencode.v31`` with the path to the human genome index you just downloaded):

.. code::

  loompy fromfq 1kPBMC.loom 1kPBMC human_GRCh38_gencode.v31 metadata.tab pbmc_1k_v3_S1_L001_R1_001.fastq.gz pbmc_1k_v3_S1_L001_R2_001.fastq.gz pbmc_1k_v3_S1_L002_R1_001.fastq.gz pbmc_1k_v3_S1_L002_R2_001.fastq.gz

Note that the fastq files are listed in pairs of R1 (read 1) and R2 (read 2) files. The index files (I1) are not used.

After about half an hour, you will have a ``1kPBMC.loom`` file with separate ``spliced`` and ``unspliced`` layers (the main matrix will be
the sum of the two), and rich metadata for both genes, cells and the sample itself stored as attributes. The output should look something like this:

.. code::

  2019-09-29 16:28:08,186 - INFO - kallisto bus -i human_GRCh38_gencode.v31/gencode.v31.fragments.idx -o /tmp/tmp7yk3rf07 -x 10xv3 -t 56 pbmc_1k_v3_S1_L001_R1_001.fastq.gz pbmc_1k_v3_S1_L001_R2_001.fastq.gz pbmc_1k_v3_S1_L002_R1_001.fastq.gz pbmc_1k_v3_S1_L002_R2_001.fastq.gz
  2019-09-29 16:28:08,307 - INFO - [index] k-mer length: 31
  2019-09-29 16:28:08,307 - INFO - [index] number of targets: 845,338
  2019-09-29 16:28:08,307 - INFO - [index] number of k-mers: 178,605,364
  2019-09-29 16:28:29,537 - INFO - [index] number of equivalence classes: 4,191,221
  2019-09-29 16:28:40,951 - INFO - [quant] will process sample 1: pbmc_1k_v3_S1_L001_R1_001.fastq.gz
  2019-09-29 16:28:40,951 - INFO -                                pbmc_1k_v3_S1_L001_R2_001.fastq.gz
  2019-09-29 16:28:40,951 - INFO - [quant] will process sample 2: pbmc_1k_v3_S1_L002_R1_001.fastq.gz
  2019-09-29 16:28:40,951 - INFO -                                pbmc_1k_v3_S1_L002_R2_001.fastq.gz
  2019-09-29 16:31:44,144 - INFO - [quant] finding pseudoalignments for the reads ... done
  2019-09-29 16:31:44,145 - INFO - [quant] processed 66,601,887 reads, 46,119,840 reads pseudoaligned
  2019-09-29 16:31:52,543 - INFO - Loading gene metadata
  2019-09-29 16:31:52,818 - INFO - Loading fragments-to-gene mappings
  2019-09-29 16:31:53,426 - INFO - Indexing genes
  2019-09-29 16:31:53,846 - INFO - Loading equivalence classes
  2019-09-29 16:32:22,273 - INFO - Mapping equivalence classes to genes
  2019-09-29 16:32:32,817 - INFO - Loading fragment IDs
  2019-09-29 16:32:33,280 - INFO - Loading BUS records
  2019-09-29 16:33:46,692 - INFO - Sorting cell IDs
  2019-09-29 16:33:49,611 - INFO - Found 46,119,840 records for 60,662 genes and 551,892 uncorrected cell barcodes.
  2019-09-29 16:33:49,611 - INFO - Correcting cell barcodes
  2019-09-29 16:35:58,753 - INFO - Found 307,677 corrected cell barcodes.
  2019-09-29 16:35:58,754 - INFO - Removing redundant reads using UMIs
  2019-09-29 16:36:45,546 - INFO - 71% sequencing saturation.
  2019-09-29 16:36:45,546 - INFO - Counting pseudoalignments for main matrix
  2019-09-29 16:36:52,752 - INFO - Found 5,027,188 UMIs.
  2019-09-29 16:36:53,536 - INFO - Counting pseudoalignments for layer 'unspliced'
  2019-09-29 16:38:00,099 - INFO - Found 2,376,590 UMIs.
  2019-09-29 16:38:00,706 - INFO - Counting pseudoalignments for layer 'spliced'
  2019-09-29 16:39:09,718 - INFO - Found 3,231,999 UMIs.
  2019-09-29 16:39:09,718 - INFO - Calling cells
  2019-09-29 16:42:32,387 - INFO - Found 1189 valid cells and ~77 ambient UMIs.
  2019-09-29 16:42:32,388 - INFO - Creating loom file '1kPBMC.loom'
  2019-09-29 16:42:32,388 - INFO - Saving

As you can see, 46,119,840 of 66,601,887 reads pseudoaligned (~70%) which is typical. The sequencing saturation was 71%, and the cell
calling algorithm found 1189 valid cells (similar to the 1,222 cells reported by cellranger). Empty beads carried a median of 77
UMIs, presumably from cell-free ambient RNA.

Connect to the loom file and examine its global attributes:

.. code::

  import loompy
  with loompy.connect("1kPBMC.loom") as ds:
    print(ds.attrs.keys())

  ['AmbientPValue', 'AmbientUMIs', 'BarcodeTotalUMIs', 'CellBarcodes', 'CreationDate', 'KallistoCommand', 'KallistoVersion', 'LOOM_SPEC_VERSION', 'NumPseudoaligned', 'NumReadsProcessed', 'RedundantReadFraction', 'SampleID', 'Saturation', 'Species', 'name', 'targetnumcells', 'technology']


...column attributes...

.. code::

  import loompy
  with loompy.connect("1kPBMC.loom") as ds:
    print(ds.ca.keys())

  ['CellID', 'TotalUMIs']


...row attributes (see the `index build instructions <https://github.com/linnarsson-lab/loompy/blob/master/notebooks/build_index.ipynb>`_ for an explanation of these)...

.. code::

  import loompy
  with loompy.connect("1kPBMC.loom") as ds:
    print(ds.ra.keys())

  ['Accession', 'AccessionVersion', 'Aliases', 'CcdsID', 'Chromosome', 'ChromosomeEnd', 'ChromosomeStart', 'CosmicID', 'DnaBindingDomain', 'FullName', 'Gene', 'GeneType', 'HgncID', 'IsTF', 'Location', 'LocationSortable', 'LocusGroup', 'LocusType', 'MgdID', 'MirBaseID', 'OmimID', 'PubmedID', 'RefseqID', 'RgdID', 'UcscID', 'UniprotID', 'VegaID']


...and layers:

.. code::

  import loompy
  with loompy.connect("1kPBMC.loom") as ds:
    print(ds.layers.keys())

  ['', 'spliced', 'unspliced']

