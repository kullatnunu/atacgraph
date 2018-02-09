Installation
============

1. Require Python 2.7 and virturalenv.

  .. Note::
    AtacGraph needs `SAMtools <http://www.htslib.org/>`_ , `deepTools <https://deeptools.readthedocs.org>`_ and
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
