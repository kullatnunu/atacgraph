#! /usr/bin/python
import pandas as pd
import numpy as np
import math
import os
import subprocess, sys
import glob
import time
import argparse
import matplotlib
import csv
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#time start
tstart = time.time()

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--promoter', type=int, default=2000)
parser.add_argument('-s', '--separate', type=int, default=1)
parser.add_argument('input_gene')
parser.add_argument('input_bam')
args = parser.parse_args()

origin_bam=args.input_bam

#ask for removing mitochondria
rm_mt = raw_input('Removing Mitochondria? (y/n): ')
rm_cp = raw_input('Removing Chlorophyll? (y/n): ')

if rm_mt == 'y':
	mt_name = raw_input('Enter Mitochondria name: ')

	subprocess.call('''samtools view %s | awk '{print$3}'>%s '''%(origin_bam,origin_bam+'_chr'),shell=True)
	a = pd.read_csv(origin_bam+'_chr',header=None)
	a.columns=['chr']
	a['chr']=a['chr'].astype(str)
	b=a.groupby(['chr']).size().reset_index(name='read')
	b_sum = b['read'].sum().astype(float)
	mt_num = b[b.chr == mt_name]
	if mt_num.empty :
		sys.exit("Mitochondria name not found.") 
	else:
		mt_num = mt_num.values[0][1]
		rmmt_ratio=(b_sum-mt_num)/b_sum
		rmmt_ratio = (rmmt_ratio*100)
		print "Original gene number: %.0f"%(b_sum)
		print "Remain gene number: %.0f ( %.2f%s gene remain)"%((b_sum-mt_num),rmmt_ratio,"%")
		print"----------------------------------------------------"

		subprocess.call('''samtools view -hq 10 %s| grep -v %s| samtools view -Sb - > %s'''%(origin_bam,mt_name,origin_bam+'_hq.bam'), shell=True)

		print "ATAC-seq_Pipeline_START"
elif rm_cp == 'y':
	cp_name = raw_input('Enter Chlorophyll name: ')

	subprocess.call('''samtools view %s | awk '{print$3}'>%s '''%(origin_bam,origin_bam+'_chr'),shell=True)
	a = pd.read_csv(origin_bam+'_chr',header=None)
	a.columns=['chr']
	a['chr']=a['chr'].astype(str)
	b=a.groupby(['chr']).size().reset_index(name='read')
	b_sum = b['read'].sum().astype(float)
	cp_name = b[b.chr == cp_name]
	if cp_name.empty :
		sys.exit("Chlorophyll name not found.")  
	else:
		cp_name = cp_name.values[0][1]
		rmcp_ratio=(b_sum-cp_name)/b_sum
		rmcp_ratio = (rmcp_ratio*100)
		print "Original gene number: %.0f"%(b_sum)
		print "Remain gene number: %.0f ( %.2f%s gene remain)"%((b_sum-cp_name),rmcp_ratio,"%")
		print"----------------------------------------------------"

		subprocess.call('''samtools view -hq 10 %s| grep -v %s| samtools view -Sb - > %s'''%(origin_bam,cp_name,origin_bam+'_hq.bam'), shell=True)

		print "ATAC-seq_Pipeline_START"

elif rm_cp == 'y' and rm_mt == 'y':
	mt_name = raw_input('Enter Mitochondria name: ')
	cp_name = raw_input('Enter Chlorophyll name: ')

	subprocess.call('''samtools view %s | awk '{print$3}'>%s '''%(origin_bam,origin_bam+'_chr'),shell=True)
	a = pd.read_csv(origin_bam+'_chr',header=None)
	a.columns=['chr']
	a['chr']=a['chr'].astype(str)
	b=a.groupby(['chr']).size().reset_index(name='read')
	b_sum = b['read'].sum().astype(float)
	mt_num = b[b.chr == mt_name]
	cp_name = b[b.chr == cp_name]
	if mt_num.empty :
		sys.exit("Mitochondria name not found.") 
	elif cp_name.empty :
		sys.exit("Chlorophyll name not found.")  
	else:
		mt_num = mt_num.values[0][1]
		rmmt_ratio=(b_sum-mt_num)/b_sum
		rmmt_ratio = (rmmt_ratio*100)
		print "Original gene number: %.0f"%(b_sum)
		print "Remain gene number: %.0f ( %.2f%s gene remain)"%((b_sum-mt_num),rmmt_ratio,"%")
		print"----------------------------------------------------"
		cp_name = cp_name.values[0][1]
		rmcp_ratio=(b_sum-cp_name)/b_sum
		rmcp_ratio = (rmcp_ratio*100)
		print "Original gene number: %.0f"%(b_sum)
		print "Remain gene number: %.0f ( %.2f%s gene remain)"%((b_sum-cp_name),rmcp_ratio,"%")
		print"----------------------------------------------------"

		subprocess.call('''samtools view -hq 10 %s| grep -vE %s| samtools view -Sb - > %s'''%(origin_bam,mt_name|cp_name,origin_bam+'_hq.bam'), shell=True)

		print "ATAC-seq_Pipeline_START"

else :
	subprocess.call('''samtools view -hq 10 %s| samtools view -Sb - > %s'''%(origin_bam,origin_bam+'_hq.bam'), shell=True)
	print "ATAC-seq_Pipeline_START"  

input_bam=origin_bam+'_hq.bam'
input_gene=args.input_gene

#clean bam files
subprocess.call('''samtools index %s'''%(input_bam),shell=True)

if glob.glob(input_gene+'.gtf'):
	print ""
else:
	subprocess.call('''gffread %s -T -o %s'''%(input_gene,input_gene+'.gtf'),shell=True)

integ_peak=input_bam+'_integ_peak_peaks.broadPeak'
long_peak=input_bam+'_long_peak_peaks.broadPeak'
short_peak=input_bam+'_short_peak_peaks.broadPeak'
peak_name = [integ_peak,long_peak,short_peak]

bam_coverage=input_bam+'_coverage.bw'
long_bam_coverage=input_bam+'_long_coverage.bw'
short_bam_coverage=input_bam+'_short_coverage.bw'
bam_coverage_name=[bam_coverage,long_bam_coverage,short_bam_coverage]

if (args.separate == 1):
	peak_name = [integ_peak]
	bam_coverage_name=[bam_coverage]
	print "*----------------------*"
	print "|Bam_read_length_filter|"
	print "*----------------------*"
	subprocess.call('samtools view -H %s > %s'%(input_bam,input_bam+'_header.sam'), shell=True)
	print "*----------------*"
	print "|Making peak file|"
	print "*----------------*"
	subprocess.call('''macs2 callpeak -t %s --nomodel --broad --shift -10 --extsize 20 -n %s'''%(input_bam,input_bam+'_integ_peak'), shell=True)
	#Making bam coverage
	subprocess.call('''bamCoverage -b %s -bs 10 --normalizeUsingRPKM --Offset 1 20 -o %s'''%(input_bam, input_bam+'_coverage.bw'),shell=True)
elif(args.separate == 2):
	peak_name = [integ_peak,long_peak,short_peak]
	bam_coverage_name=[bam_coverage,long_bam_coverage,short_bam_coverage]
	print "*----------------------*"
	print "|Bam_read_length_filter|"
	print "*----------------------*"
	subprocess.call('samtools view -H %s > %s'%(input_bam,input_bam+'_header.sam'), shell=True)
	subprocess.call('''samtools view %s -f 0x02| awk '{ if ( $9>200||$9<-200 ) print$0} '>%s'''%(input_bam,input_bam+'_long.sam'),shell=True)
	subprocess.call('''samtools view %s -f 0x02| awk '{ if ( $9<=200 && $9>= -200) print$0} '>%s'''%(input_bam,input_bam+'_short.sam'),shell=True)
	subprocess.call('''cat %s %s|samtools view -bS - > %s'''%(input_bam+'_header.sam',input_bam+'_long.sam',input_bam+'_long.bam'),shell=True)
	subprocess.call('''cat %s %s|samtools view -bS - > %s'''%(input_bam+'_header.sam',input_bam+'_short.sam',input_bam+'_short.bam'),shell=True)
	print "*--------------*"
	print "|Index_Bam_File|"
	print "*--------------*"
	subprocess.call('''samtools index %s'''%(input_bam+'_long.bam'),shell=True)
	subprocess.call('''samtools index %s'''%(input_bam+'_short.bam'),shell=True)
	print "*----------------*"
	print "|Making peak file|"
	print "*----------------*"
	subprocess.call('''macs2 callpeak -t %s --nomodel --broad --shift -10 --extsize 20 -n %s'''%(input_bam,input_bam+'_integ_peak'), shell=True)
	subprocess.call('''macs2 callpeak -t %s --format BAMPE --broad -n %s'''%(input_bam+'_long.bam',input_bam+'_long_peak'), shell=True)
	subprocess.call('''macs2 callpeak -t %s --format BAMPE --broad -n %s'''%(input_bam+'_short.bam',input_bam+'_short_peak'), shell=True)	
	#Making bam coverage
	subprocess.call('''bamCoverage -b %s -bs 10 --normalizeUsingRPKM --Offset 1 20 -o %s'''%(input_bam, input_bam+'_coverage.bw'),shell=True)
	subprocess.call('''bamCoverage -b %s -bs 10 --normalizeUsingRPKM -e -o %s'''%(input_bam+'_long.bam', input_bam+'_long_coverage.bw'),shell=True)
	subprocess.call('''bamCoverage -b %s -bs 10 --normalizeUsingRPKM -e -o %s'''%(input_bam+'_short.bam', input_bam+'_short_coverage.bw'),shell=True)
else:
	print "Please enter separation number"



#making histograms
subprocess.call('''samtools view %s | awk '{print$9}'>%s '''%(input_bam,input_bam+'_readlen'),shell=True)
readlen = pd.read_csv(input_bam+'_readlen',header=None)
readlen.columns=['readlen']
readlen=readlen[readlen['readlen']>0]
readlen=readlen.reset_index()

plt.style.use('ggplot')
fig=plt.figure()
ax1=fig.add_subplot(1,1,1)
n, bins, patches = ax1.hist(readlen['readlen'], bins=500, normed=False, color='orange',range=(0,1000))
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')
plt.xlabel('Read_length')
plt.ylabel('Number_of_Reads')
fig.suptitle(input_bam,fontsize=14, fontweight='bold')
ax1.set_title('Frequency_Distributions')
plt.savefig(input_bam+'_readlen.png',dpi=400,bbox_inches='tight')
plt.close(fig)


#annotation name
gene = "gene_body"
exon = "exons"
intron = "introns"
utr3 = "3utr"
utr5 = "5utr"
cds = "cds"
promoter = "gene_promoter"
igr = "gene_igr"
annotation_name=[promoter,gene,exon,intron,utr5,cds,utr3,igr]





#clean gtf
gene_gtf=pd.read_csv(input_gene+'.gtf', header=None, sep="\t")
gene_gtf=gene_gtf[gene_gtf[4]>gene_gtf[3]]
gene_gtf.to_csv(input_gene+'gene.gtf',sep="\t", index=False, header=None, quoting=csv.QUOTE_NONE)

print "*--------------------------------*"
print "|Extract UTR, exon, cds from gene|"
print "*--------------------------------*"

subprocess.call('''../extract_transcript_regions.py -i %s -o %s --gtf'''%(input_gene+'gene.gtf',input_gene+'.gtf'), shell=True)

print "*-------------------------------------*"
print "|Convert this blockbed (bed12) to bed6|"
print "*-------------------------------------*"

ano_filename = [exon,intron,utr5,cds,utr3]

for i in ano_filename:
	subprocess.call('''cat %s | bed12ToBed6 -i stdin -n > %s'''%(input_gene+'.gtf'+'_'+i+'.bed',input_gene+'.gtf'+'_'+i+'_bed6.bed'),shell=True)


#find gene_body.bed
genes = pd.read_csv(input_gene+'gene.gtf', header=None, sep="\t")
genes.columns=['chr','unknow', 'exon', 'g_str', 'g_end', 'g_score', 'g_dir','.', 'gene_name']
genes.chr=genes.chr.astype(str)
genes=genes[genes.exon=='exon']
gene_col=genes['gene_name'].str.split(';', expand=True)
gene_col.columns=gene_col.iloc[1,:]
gene_id = gene_col.filter(regex='gene_id')
gene_id = gene_id.ix[:,0].str.split(' ', expand=True)
gene_id[2] = gene_id[2].map(lambda x: x.lstrip('"').rstrip('"'))
gene_id.columns=['num','g_name','gene_id']
gene_bed = genes.ix[:,['chr', 'g_str', 'g_end', 'g_score','g_dir']].join(gene_id['gene_id'])
gene_bed = gene_bed.ix[:,['chr', 'g_str', 'g_end', 'gene_id','g_score','g_dir']]
gene_bed=gene_bed.drop_duplicates(subset=['g_str','g_end'],keep='first')
gene_bed=gene_bed.sort_values(['chr','g_str'],ascending=[True,True])
gene_bed=gene_bed.drop_duplicates(subset=['g_str'],keep='last')
gene_group=gene_bed.groupby(['chr','gene_id','g_score','g_dir']).agg({'g_str':'min', 'g_end':'max'}).reset_index()
gene_group = gene_group.drop_duplicates(subset=['g_str','g_end'],keep='first')
gene_body=gene_group.sort_values(['chr','g_str'],ascending=[True,True])
gene_body=gene_body.drop_duplicates(subset=['g_str'],keep='last')
gene_body = gene_body.ix[:,['chr', 'g_str', 'g_end', 'gene_id','g_score','g_dir']]
gene_body=gene_body[~gene_body.gene_id.str.contains('MI')]
gene_body.to_csv(input_gene+'.gtf'+'_gene_body_bed6.bed', sep='\t',index=False, header=None)
genome_bp=gene_body.groupby(['chr']).agg({'g_end':'max'}).reset_index()
genome_bp = genome_bp['g_end'].sum()
genome_bp=genome_bp*1.0


#find promoter.bed
gene_body['pro_str'] = np.where(gene_body.g_dir == '+', gene_body.g_str - args.promoter, gene_body.g_end - 0)
gene_body['pro_end'] = np.where(gene_body.g_dir == '+', gene_body.g_str + 0, gene_body.g_end + args.promoter)
num = gene_body._get_numeric_data()
num[num<0]=0
gene_promoter = gene_body.ix[:, ['chr','pro_str','pro_end','gene_id','g_score','g_dir']]
gene_promoter.columns=['chr','g_str','g_end','gene_id','g_score','g_dir']

gene_promoter.to_csv(input_gene+'.gtf'+"_gene_promoter_bed6.bed", sep='\t',index=False, header=None)

genome=gene_body.groupby(['chr']).agg({'g_end':'max'}).reset_index()
genome.to_csv(input_gene+'.gtf'+"_genome.bed", sep="\t", index=False, header=None)

gbname=[promoter,gene,exon,intron,utr5,cds,utr3]
for i in gbname:
	subprocess.call('''bedtools sort -i %s|bedtools merge -c 4,5,6 -o collapse,collapse,collapse >%s '''%(input_gene+'.gtf'+'_'+i+'_bed6.bed',input_gene+'.gtf'+'_'+i+'_merge.bed'),shell=True)


subprocess.call('''bedtools complement -i %s -g %s > %s '''%(input_gene+'.gtf'+'_'+'gene_body'+'_merge.bed',input_gene+'.gtf'+"_genome.bed", input_gene+'.gtf'+'_'+'gene_igr'+'.bed'),shell=True)
geneigr = pd.read_csv(input_gene+'.gtf'+'_'+'gene_igr'+'.bed', header=None, sep="\t")
geneigr['.'],geneigr['...'],geneigr['..']=['.','.','.']
geneigr.to_csv(input_gene+'.gtf'+'_'+'gene_igr'+'_bed6.bed', header=None,index=False, sep="\t")
subprocess.call('''bedtools sort -i %s|bedtools merge -c 4,5,6 -o collapse,collapse,collapse >%s '''%(input_gene+'.gtf'+'_'+'gene_igr'+'_bed6.bed',input_gene+'.gtf'+'_'+'gene_igr'+'_merge.bed'),shell=True)



#peak bp count
def peak_bp(peak):
	peak = pd.read_csv(peak, header=None, sep="\t").ix[:,0:5]
	peak.columns = ['chr', 'peak_str', 'peak_end', 'peak_id', 'peak_value', 'peak_dir']
	peak=peak.drop_duplicates(subset=['peak_str', 'peak_end'], keep='first')
	peak['bp'] = (peak.peak_end - peak.peak_str)
	peak_bp = peak.bp.sum()
	return peak_bp


#annotation bp count
def ano_bp(anno):
	ano = pd.read_csv(input_gene+'.gtf'+'_'+anno+'_merge.bed', header=None, sep="\t")
	ano.columns=['chr','g_str','g_end','gene_id','g_id','g_dir']
	ano=ano.drop_duplicates(subset=['g_str', 'g_end'], keep='first')
	ano['bp'] = (ano.g_end - ano.g_str)
	ano_bp = ano.bp.sum()
	return ano_bp


print "---------------------------------------"
print "|Making annotation_peak_associate file|"
print "---------------------------------------"

file_exist=[]
for i in annotation_name:
	for j in peak_name:
		subprocess.call('''bedtools intersect -nonamecheck -a %s -b %s -wo > %s'''%(input_gene+'.gtf'+'_'+i+'_merge.bed',j,j+'_'+i+'.txt'),shell=True)
	if os.stat( j+'_'+i+'.txt' ).st_size == 0 :
		file_exist.append(1)

#annotation_peak bp count
def annopeakbp(peak,anno):
	anno_peak_asso = pd.read_csv(peak+'_'+anno+'.txt', header=None, sep="\t")
	anno_peak_asso_bp = anno_peak_asso[len(anno_peak_asso.columns)-1].sum()
	anno_peak_asso_bp = float(anno_peak_asso_bp)
	return anno_peak_asso_bp


#making Enrichment Graph
def enrichment_num(peak,anno):
	try:
		enrichment=math.log((annopeakbp(peak,anno)/peak_bp(peak))/(ano_bp(anno)/genome_bp),2)
		return enrichment
	except pd.errors.EmptyDataError:
		return 0


#making associate table
def associate(inpeak):
	peak = pd.read_csv(inpeak, header=None, sep="\t").ix[:,0:5]
	gene_body = pd.read_csv(input_gene+'.gtf'+'_gene_body_merge.bed', header=None, sep="\t")
	peak.columns = ['chr', 'peak_str', 'peak_end', 'peak_name','peak_value','peak_dir']
	peak.chr=peak.chr.astype(str)
	gene_body.columns = ['chr', 'gbed_str', 'gbed_end', 'gene_id', 'gene_value', "gene_dir"]
	gene_body['pro_str'] = np.where(gene_body.gene_dir == '+', gene_body.gbed_str -2000, gene_body.gbed_end -0)
	gene_body['pro_end'] = np.where(gene_body.gene_dir == '+', gene_body.gbed_str +0, gene_body.gbed_end +2000)
	combined = pd.merge(gene_body, peak, on='chr')
	combined['Genebody'] = np.where((combined.gbed_str < combined.peak_end) & (combined.gbed_end > combined.peak_str), 1, 0)
	combined['Promoter'] = np.where((combined.pro_str < combined.peak_end) & (combined.pro_end > combined.peak_str), 1, 0)
	summary = combined[(combined.Promoter > 0) | (combined.Genebody > 0)]
	s1 = summary.drop(summary.columns[6:13], axis = 1).drop(summary.columns[4], axis = 1)
	s1_group=s1.groupby(['gene_id', 'chr', 'gbed_str', 'gbed_end', 'gene_dir']).agg({'Genebody':'sum', 'Promoter':'sum'}).reset_index().sort_values(["chr","gbed_str"])
	s1_group.to_csv(inpeak+"_gene"+"_summary_table.xls", sep='\t',index=False)


def coverage_heatmap(coverage):
	subprocess.call('''computeMatrix scale-regions -S %s -R %s --missingDataAsZero -bs 10 -a 1000 -b 1000 -out %s --outFileNameMatrix %s'''%(coverage,input_gene+'.gtf'+'_gene_body_merge.bed',coverage+'gene_body'+'.matrix.gz',coverage+'gene_body'+'.matrix.txt'),shell=True)
	subprocess.call('''plotHeatmap -m %s -out %s --legendLocation none'''%(coverage+'gene_body'+'.matrix.gz',coverage+'gene_body_heatmap.png'),shell=True)


#making table and peak_bed to bedgraph
print "-----------------------------------------------"
print "|Making Summary Table and Peak Density Heatmap|"
print "-----------------------------------------------"
for i in peak_name:
	summary_table = associate(i)

for i in bam_coverage_name:
	bamcoveragegraph=coverage_heatmap(i)

print "*--------------------------------------------*"
print "| Summary Table and Peak Density Heatmap Done|"
print "*--------------------------------------------*"


#making fold enrichment graph
print "*-----------*"
print "|Making plot|"
print "*-----------*"


plt.style.use('ggplot')
annotationname = ['Promoter','Genebody','Exon','Intron','5UTR','CDS','3UTR','IGR']
annotationname_index = range(len(annotationname))

for i in peak_name:
	enrichment_data = []
	for j in annotation_name:
		enrichment_data.append(enrichment_num(i,j))
	fold_enrich = enrichment_data
	fig = plt.figure()
	ax1 = fig.add_subplot(1,1,1)
	colors=['firebrick','darkblue','gold','darkviolet','darkgreen','darkcyan','c','purple']
	ax1.bar(annotationname_index, fold_enrich, align='center',color=colors)
	ax1.xaxis.set_ticks_position('bottom')
	ax1.yaxis.set_ticks_position('left')
	plt.ylim(ymin=-2,ymax=6)
	plt.xticks(annotationname_index, annotationname, fontsize='small')
	plt.ylabel("Fold Enrichment (log2)")
	plt.title(i+'_Fold_Enrichment')
	plt.savefig(i+'_Fold_Enrichment'+'.png',dpi=400,bbox_inches='tight')
	plt.close(fig)
	fe_table=pd.DataFrame([fold_enrich],columns=annotationname)
	fe_table.to_csv(i+'_Fole_Enrichment_Table', index=None, sep="\t")


#Making tophat junction for IGV
subprocess.call('''samtools view %s|awk '{if ($7 == "=" && $9>0){print $1"\t"$3"\t"$4"\t"$7"\t"$8"\t"$9}} ' >%s'''%(input_bam,input_bam+'paired'), shell="True")

with open(input_bam+'_junction.bed', 'w') as c:
	track_name = '''track name=junctions description="TopHat junctions"'''
	blank=""
	c.write("%s\n%s"%(track_name,blank))

c1_pe = pd.read_csv(input_bam+'paired', header=None, sep="\t")
c1_pe[7] = c1_pe[2]+c1_pe[4]
c1_pe = c1_pe.ix[:,(1,2,4,0)]
c1_pe.columns=['chr','str','end','name']
bin_size = int(10)
c1_pe.str=c1_pe.str//bin_size*bin_size
c1_pe.end=c1_pe.end//bin_size*bin_size
c1_score=c1_pe.groupby(['chr','str','end']).count().reset_index()
c1_score['chrname']=c1_score.chr.map(str)+c1_score.str.map(str)
c1_score['dir']='+'
c1_score['thickstart']=c1_score.str
c1_score['thickend']=c1_score.end
c1_score['block_count']=2
c1_score['block_size']="1,1"
c1_score['location']=c1_score.end-c1_score.str-1

#rgb_200_up_below
c1_score['rgb']=np.where(c1_score['location']>200,'255,0,0','0,0,255')
c1_score['zero']=0
c1_score['block_location']=c1_score['zero'].map(str)+","+c1_score['location'].map(str)
c1_junction_pd = c1_score.ix[:,('chr','str','end','chrname','name','dir','thickstart','thickend','rgb','block_count','block_size','block_location')]
c1_junction_pd.to_csv(input_bam+'_junction.bed',mode='a', header=None, index=None, sep="\t")

subprocess.call('''rm *bed6*|rm *3utr.bed|rm *5utr.bed|rm *_start.bed|rm *_cds.bed|rm *_codingexons.bed|rm *_codingintrons.bed|rm *_exons.bed|rm *_igr.bed|rm *_genome.bed|rm *_introns.bed|rm *_noncodingexons.bed|rm *_noncodingintrons.bed|rm *bam.bai|rm *sam|rm *bam_chr|rm *gappedPeak|rm *peaks.xls|rm *.bampaired ''',shell=True)

tend = time.time()#time stop

if rm_mt == 'y' and mt_num:
    print "Original gene number: %.0f"%(b_sum)
	print "Remain gene number: %.0f ( %.2f%s gene remain)"%((b_sum-mt_num),rmmt_ratio,"%")	
	print " "
elif rm_cp == 'y' and cp_name:
    print "Original gene number: %.0f"%(b_sum)
    print "Remain gene number: %.0f ( %.2f%s gene remain)"%((b_sum-cp_name),rmcp_ratio,"%")
    print " "


print "***----------Processing Time: %s seconds ----------***" %(tend-tstart)
