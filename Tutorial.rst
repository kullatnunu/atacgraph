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
   
   Summary table of ATAC-seq peak abundance
   
   Fold enrichment analysis of chromatin accessibility
   
   Heatmap depicting accessibility for gene
   
   Profile map of accessibility for genes
   
   Visualization of ATAC fragments IGV
  

