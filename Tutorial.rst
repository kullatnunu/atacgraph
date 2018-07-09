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
  
  
 
 
   * **Fold enrichment analysis of chromatin accessibility**
   
   .. image:: https://github.com/kullatnunu/atacgraph/blob/master/github/Ctrl_1_chr1.bam_hq.bam_integ_peak_peaks.broadPeak_Fold_Enrichment.jpg




   *  **Heatmap depicting accessibility for gene & Profile map of accessibility for genes**
   
   .. image:: https://github.com/kullatnunu/atacgraph/blob/master/github/Ctrl_1_chr1.bam_hq.bam_coverage.bwgene_body_heatmap.jpg
   
   
   
   
   *  **Visualization of ATAC fragments IGV**
     
      Ctrl_1_chr1.bam_hq.bam_junction.bed
     
     
   * **Other files**
     
=========================================================================  =====================================
Output file name	                                                         Descrition
=========================================================================  =====================================
Ctrl_1_chr1.bam_hq.bam	                                                   qulified bam file, (high quality score is at leat 10)
Ctrl_1_chr1.bam_hq.bam_coverage.bw	                                       bigwiggle format of input file 
Ctrl_1_chr1.bam_hq.bam_coverage.bwgene_body.matrix.gz	                     value of Heatmap depicting accessibility for gene
Ctrl_1_chr1.bam_hq.bam_coverage.bwgene_body.matrix.txt	                   value of Heatmap depicting accessibility for gene
Ctrl_1_chr1.bam_hq.bam_integ_peak_peaks.broadPeak	                         The location of ATAC-seq peaks
Ctrl_1_chr1.bam_hq.bam_integ_peak_peaks.broadPeak_3utr.txt	               The intersection site between 3UTR and peaks
Ctrl_1_chr1.bam_hq.bam_integ_peak_peaks.broadPeak_5utr.txt	               The intersection site between 5UTR and peaks
Ctrl_1_chr1.bam_hq.bam_integ_peak_peaks.broadPeak_cds.txt	                 The intersection site between CDS and peaks
Ctrl_1_chr1.bam_hq.bam_integ_peak_peaks.broadPeak_exons.txt	               The intersection site between exon and peaks
Ctrl_1_chr1.bam_hq.bam_integ_peak_peaks.broadPeak_gene_body.txt	           The intersection site between genebody and peaks
Ctrl_1_chr1.bam_hq.bam_integ_peak_peaks.broadPeak_gene_igr.txt	           The intersection site between IGR and peaks
Ctrl_1_chr1.bam_hq.bam_integ_peak_peaks.broadPeak_gene_promoter.txt	       The intersection site between promoter and peaks
Ctrl_1_chr1.bam_hq.bam_integ_peak_peaks.broadPeak_introns.txt	             The intersection site between intron and peaks
Ctrl_1_chr1.bam_hq.bam_integ_peak_peaks.broadPeak_gene_summary_table.xls	 The number of peaks in promoter and genebody
Ctrl_1_chr1.bam_hq.bam_readlen	                                           The number of ATAC-seq read length distribution
genes_demo.gtf.gtf_3utr_merge.bed	                                         The sites of 3UTR 
genes_demo.gtf.gtf_5utr_merge.bed	                                         The sites of 5UTR 
genes_demo.gtf.gtf_cds_merge.bed	                                         The sites of CDS 
genes_demo.gtf.gtf_exons_merge.bed	                                       The sites of exon
genes_demo.gtf.gtf_gene_body_merge.bed	                                   The sites of genebody 
genes_demo.gtf.gtf_gene_igr_merge.bed	                                     The sites of IGR 
genes_demo.gtf.gtf_gene_promoter_merge.bed	                               The sites of promoter
genes_demo.gtf.gtf_introns_merge.bed	                                     The sites of intron
=========================================================================  =====================================


