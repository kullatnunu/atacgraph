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
