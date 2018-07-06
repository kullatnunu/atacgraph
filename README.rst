.. image:: https://img.shields.io/badge/python-2.7-blue.svg

ATAC-graph
=========

ATAC-seq pipeline for making analysis graph

ATAC-graph will produce 6 analysis results:

* ATAC-seq read length distribution
* Summary table of ATAC-seq peak abundance
* Fold enrichment analysis of chromatin accessibility
* Heatmap depicting accessibility for gene
* Profile map of accessibility for genes
* Visualization of ATAC fragments IGV

ATAC-graph Pipeline
=========

.. image:: https://github.com/RitataLU/atacgraph/blob/master/figure1.png

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

   (1) Enter your bash profile
   
   ::
  
   $ vi ~/.bash_profile
   
   (2) Add your ATAC-graph path to the PATH environment variable.
  
   ::
   
   $ PATH=$PATH:(put your ATAC-graph file path here)
   $ source ~/.bash_profile


Running ATAC-graph
==================
  .. Note::
  * Input bam file is input.bam, mitochrondria name is mito_name
  
  ::
  


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

