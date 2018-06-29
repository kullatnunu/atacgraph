.. image:: https://img.shields.io/badge/python-2.7-blue.svg

ATAC-graph
=========

ATAC-seq pipeline for making analysis graph

ATAC-graph will produce 6 analysis results:

* Summary table
* Histogram
* Heatmap
* Profile map
* Fold enrichment graph
* Visualization of ATAC fragments IGV

System Requirement
==================

* Python 2.7

.. Note::
    AtacGraph needs `SAMtools <http://www.htslib.org/>`_ , `deepTools <https://deeptools.readthedocs.org>`_ and
    `BEDtools <http://bedtools.readthedocs.org/>`_ to run the script, we will need to install them on your server.

Installation
============

1. Download the source code and install the requirements.

  ::

  $ git clone https://github.com/kullatnunu/atacgraph.git
  $ pip install -r atacgraph/base.txt

  pip will install the following packages:

  * `NumPy <http://www.numpy.org/>`_
  * `matplotlib <http://matplotlib.org/>`_
  * `pandas <http://matplotlib.org/>`_
  
2. Add your ATAC-graph path to the PATH.


User's guide
============
  .. Note::
  * Before input bam file, please remove mitochondria
  * Input bam file is input.bam, mitochrondria name is mito_name
  
  ::
  
  $ samtools view -hq 10 input.bam| grep -v mito_name| samtools view -Sb - > output.bam

Arguments
---------
-p, --promoter <INT>
--------------------
  Size of promoter, default is 2,000

Input
-----
gtf, gff
--------
  Input gene annotation GTF or GFF file

bam
---
  Input atac-seq bam file


Tutorial
========
Demo file
---------

1. Into the atacgraph file and download the sample input file

  ::

  $ cd atacgraph
  $ wget -O data.tar.gz https://github.com/kullatnunu/atacgraph/blob/master/demo/data.tar.gz?raw=true
  $ tar xvfz data.tar.gz
  $ cd data

2. Run atacgraph script

  ::

  $ atac_graph.py genes_demo.gtf Ctrl_1_chr1.bam

