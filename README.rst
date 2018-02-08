# Atacgraph
ATAC-seq pipeline


Installation
============

1. Install Python 2.7 and virturalenv.

  .. Note::
    AtacGraph needs `SAMtools <http://www.htslib.org/>`_ and
    `BEDtools <http://bedtools.readthedocs.org/>`_ to run the script, we will need to install them on your server.

2. Create the virtual environment on your disk, and activate the environment.

  ::

  $ virtualenv --no-site-packages --python=python2.7 atacgraph_env
  $ source atacgraph_env/bin/activate


3. Download the source code and install the requirements.

  ::

  $ git clone https://github.com/kullatnunu/atacgraph.git
  $ pip install -r atacgraph/base.txt

  pip will install the following packages:

  * `NumPy <http://www.numpy.org/>`_
  * `matplotlib <http://matplotlib.org/>`_
  * `pandas <http://matplotlib.org/>`_
  
4. Add your MethGo path to the PATH environment variable.


Tutorial
========

1.Into the atacgraph file and download the sample input file

::

$ cd atacgraph
$ wget -O data.tar.gz https://github.com/kullatnunu/atacgraph/blob/master/demo/data.tar.gz?raw=true
$ tar xvfz data.tar.gz
$ cd data

2.Run atacgraph script

::

$ atac_graph.py genes_demo.gtf Ctrl_1_chr1.bam
