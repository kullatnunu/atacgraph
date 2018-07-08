Tutorial
========
Run demo 
---------

1. Download the input sample in the atacgraph folder

  ::

  $ cd atacgraph
  $ wget -O data.tar.gz https://github.com/kullatnunu/atacgraph/blob/master/demo/data.tar.gz?raw=true
  $ tar xvfz data.tar.gz
  $ cd data

2. Run atacgraph script

  ::

  $ atac_graph.py genes_demo.gtf Ctrl_1_chr1.bam [-p 2000]
  
3. Output
  
   ATAC-seq read length distribution
.. image:: https://github.com/kullatnunu/atacgraph/blob/master/github/Ctrl_1_chr1.bam_hq.bam_readlen.png

   
   Summary table of ATAC-seq peak abundance
   
   
   Fold enrichment analysis of chromatin accessibility
.. image:: https://github.com/kullatnunu/atacgraph/blob/master/github/Ctrl_1_chr1.bam_hq.bam_integ_peak_peaks.broadPeak_Fold_Enrichment.png

   Heatmap depicting accessibility for gene & Profile map of accessibility for genes
.. image:: https://github.com/kullatnunu/atacgraph/blob/master/github/Ctrl_1_chr1.bam_hq.bam_coverage.bwgene_body_heatmap.png
   
   Visualization of ATAC fragments IGV
  

