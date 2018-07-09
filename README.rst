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

* `SAMtools <http://www.htslib.org/>`_ 
* `deepTools <https://deeptools.readthedocs.org>`_
* `BEDtools <http://bedtools.readthedocs.org/>`_ 


* Python Modules 'Numpy', 'pandas' and 'Metplotlib'. To install the packages, use the following commands on an UNIX terminal:
  
  ::

  $ pip install numpy
  $ pip install pandas
  $ pip install matplolib
  
  
Installation
============

1. Download the source code and install the requirements.

  ::

  $ git clone https://github.com/kullatnunu/atacgraph.git
  $ pip install -r atacgraph/base.txt

  
2. Add your ATAC-graph path to the PATH.

   (1) Edit bash profile
  
   ::
  
   $ vi ~/.bash_profile
   
   (2) Add ATAC-graph path to the PATH environment variable.
 
   ::
  
   $ PATH=$PATH:(ATAC-graph file path)
   $ source ~/.bash_profile
   


Running ATAC-graph
==================

Running atacgraph.py
--------------------
Usage:
  
::

$ python atacgraph.py input_gtf_file input_bam_file
Usage: python atacgraph.py [input_gtf_file] [input_bam_file] [options]  

> Removing mitochondria? (y/n): 

  Enter "y" or "n" to remove or not remove the mitochondria. 

> Enter mitochondria name:

  Enter the mitochondria name that you want to remove.

  ex:"chrM"

Arguments
---------
-p, --promoter <INT>
--------------------
  Size of promoter, default is 2,000 bp before transcription start site

Input
-----
1. gene annotation
-----
  gene annotation in gtf or gff
  
2. bam
-----
  atac-seq bam file after mapping
  
Link
====

`Tutorial <https://github.com/kullatnunu/atacgraph/blob/master/Tutorial.rst/>`_ 
---------

