Atacgraph
=========

ATAC-seq pipeline



Installation
============

1. Install Python 2.7 and virturalenv.

  .. Note::
    AtacGraph needs `SAMtools <http://www.htslib.org/>`_ and
    `BEDtools <http://bedtools.readthedocs.org/>`_ to run the script, we will need to install them on your server.

2. Download the source code and install the requirements.

  ::

  $ git clone https://github.com/kullatnunu/atacgraph.git
  $ pip install -r atacgraph/base.txt

  pip will install the following packages:

  * `NumPy <http://www.numpy.org/>`_
  * `matplotlib <http://matplotlib.org/>`_
  * `pandas <http://matplotlib.org/>`_
  
3. Add your AtacGraph path to the PATH.



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

User's guide
============

Input
-----
gtf, gff
--------
Input gene annotation gtf or gff file

Arguments
---------
-p, --promoter <INT>
--------------------
Size of promoter, default is 2,000
