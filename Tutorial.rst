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
  
   * **ATAC-seq read length distribution**
  
   .. image:: https://github.com/kullatnunu/atacgraph/blob/master/github/Ctrl_1_chr1.bam_hq.bam_readlen.jpg
   
   
   * **Summary table of ATAC-seq peak abundance** (The number of peaks in promoter and genebody)
    
     Ctrl_1_chr1.bam_hq.bam_integ_peak_peaks.broadPeak_gene_summary_table.xls
  
  
  
   * Fold enrichment analysis of chromatin accessibility
   
   .. image:: https://github.com/kullatnunu/atacgraph/blob/master/github/Ctrl_1_chr1.bam_hq.bam_integ_peak_peaks.broadPeak_Fold_Enrichment.jpg


   *  **Heatmap depicting accessibility for gene & Profile map of accessibility for genes**
   
   .. image:: https://github.com/kullatnunu/atacgraph/blob/master/github/Ctrl_1_chr1.bam_hq.bam_readlen.jpg
   
   
   *  **Visualization of ATAC fragments IGV**
     
     Ctrl_1_chr1.bam_hq.bam_junction.bed
     
     
   * other files:
     
  | Command | Description |
  | ---     | --- |
  | `git status` | List all *new or modified* files |
  | `git diff`   | Show file differences that **haven't been** staged |

