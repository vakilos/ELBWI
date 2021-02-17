print("Import modules and set variables!")
from functions import *
import sys, os, re, argparse, subprocess, itertools
import numpy as np
import pandas as pd
import skbio.diversity as div
from Bio import SeqIO
class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

parser=argparse.ArgumentParser(description="This script will process raw ASV data..")
parser.add_argument("-i", dest="rootDir", help="Project directory should contain qiime2_dada2 subdirectory)", required=True)
args=parser.parse_args()
### Global variables
rootDir = args.rootDir  
if rootDir.endswith('/') is False:
   rootDir += "/"
tmpDir = rootDir+'tmp/'
if os.path.isdir(tmpDir) is False:
   os.mkdir(tmpDir)
figDir = rootDir+"figures/"
if os.path.isdir(figDir) is False:
   os.mkdir(figDir)
inputDir = rootDir+'qiime2_dada2/'
os.chdir(inputDir)
remove_asvs = list() ### give here asv names (i.e. ASV_19) to remove from the dataset if you like
minimum_reads = 600

###########################################################################################################################################################
print('\nImport and prepare data....')
asv_df = pd.read_csv("final_asv_table.tsv", sep='\t', index_col='SampleID')
tax_df = pd.read_csv("taxonomy_table.tsv", sep='\t', index_col='asvid')
fasta = list(SeqIO.parse('final_asv_seqs.fasta','fasta'))
stats_df = pd.read_csv('ReadStats.csv')

meta_df = pd.read_csv(rootDir+'ELBWI_metatable.tsv', sep='\t', index_col='SampleID')
nut_df  = pd.read_csv(rootDir+'nutrition.tsv', sep='\t')
nut_df = nut_df.merge(meta_df.drop_duplicates('Patient'), on="Patient", how='left')


"""Change Escherichia taxonomy to Escherichia/Shigella"""
esc = re.compile(r'(?<=.)*Escherichia')
tax_df = tax_df.reset_index()
tax_df['ASVg'] =   tax_df.ASVg.apply(lambda x:re.sub(esc, 'Escherichia/Shigella',x))
tax_df['genus'] = tax_df.genus.apply(lambda x:re.sub(esc, 'Escherichia/Shigella',x))
tax_df.set_index('asvid', inplace=True)


print("Minimum confidence for taxonomic clasifficatiohn is 70...")


""" Remove ASVs that are shorter than expected  """
fasta_lengths = [len(rec.seq) for rec in fasta]
length_cutoff = 0.8*466
Newfasta = [rec for rec in fasta if len(rec.seq) >= length_cutoff]

SeqIO.write(Newfasta, 'final_asv_seqs_clean.fasta', 'fasta')
print('\n'+color.BOLD+color.RED+'Remove sequences that their length is < 80% of the expected amplicon length(466bp))'+color.END)
print(color.BOLD+'{} short ASV sequences were removed from the dataset\n'.format(len(fasta)-len(Newfasta))+color.END)
asv_df = asv_df[[rec.id for rec in Newfasta]]

## update stats table
tmp = asv_df.sum(axis=1).to_frame().rename(columns={0:"After short amplicons removal (<80% exp. length)"}).reset_index().rename(columns={'SampleID':'sample-id'})
stats_df = stats_df.merge(tmp, on='sample-id', how='left')

"""     Remove ASVs that are present in the negative controls         """
neg_genera  =  ['Caldalkalibacillus', 'Nesterenkonia']
print(asv_df.columns)
neg_asvs = [c for c in asv_df.columns if tax_df.loc[c,'genus'] in neg_genera]
neg_asvs_read_per = asv_df[neg_asvs].sum().sum()*100 / asv_df.sum().sum()
neg_asvs_pre = asv_df[neg_asvs].astype(bool).any(axis=1).sum()  *100/ asv_df.shape[0]
print(color.BOLD+color.RED+"Remove ASVs which were detected in high abundance in the negative controls after blastn against the 16S Ribosomal RNA NCBI database."+color.END)
print(color.BOLD+"{} ASVs were removed from the dataset ({:.2f}% of total reads, detected in {:.2f}% of the samples)".format(len(neg_asvs),neg_asvs_read_per,neg_asvs_pre)+color.END)
print(tax_df.loc[neg_asvs,'ASVg'].values,'\n')
asv_df.drop(columns= neg_asvs, inplace=True)
tmp = asv_df.sum(axis=1).to_frame().rename(columns={0:"negative control ASV removal"}).reset_index().rename(columns={'SampleID':'sample-id'})
stats_df = stats_df.merge(tmp, on='sample-id', how='left') ### update stats table

"""      Remove samples that are not assigned to individuals in the metadata """
non_subject_samples = [i for i in asv_df.index if i not in meta_df.index]
asv_df.drop(index=non_subject_samples, inplace=True)
print(color.BOLD+color.RED+"Remove non-subject samples from the dataset"+color.END)
print(color.BOLD+"{} samples were removed...".format(len(non_subject_samples))+color.END)
print(non_subject_samples,'\n')



print(color.BOLD+color.RED+"Remove ASVs based on read count for a given genus between my dataset and other datasets in the same run "
      "ASVs with much higher count that are not detected in the first (check_cross_contamination.py)"+color.END)
remove_cross_genus = ["Faecalibacterium",'Cetobacterium','Clostridioides','Planoglabratella'] ### selected from check_cross_contamination.py
remove_cross_unassigned = ['ASV_127', 'ASV_176', 'ASV_60', 'ASV_82', 'ASV_92', 'ASV_93', 'ASV_94'] ### those are only detected in the second run
remove_cross_phylum = ['Cyanobacteria']
remove_cross_asvs = [t for t in tax_df.loc[tax_df['phylum'].isin(remove_cross_phylum)].index if t in asv_df.columns] + [t for t in tax_df.loc[tax_df['asv'].isin(remove_cross_unassigned)].index if t in asv_df.columns]
remove_cross_asvs += [t for t in tax_df.loc[tax_df['genus'].isin(remove_cross_genus)].index if t in asv_df.columns]
cross_read_per = asv_df[remove_cross_asvs].sum().sum() *100 / asv_df.sum().sum()
cross_read_pre = asv_df[remove_cross_asvs].astype(bool).any(axis=1).sum() *100 / asv_df.shape[0]
print(color.BOLD+"{} ASVs removed as cross-contaminants ({:.2f}% of total reads, detected in {:.2f}% of the samples)\n".format(len(remove_cross_asvs), cross_read_per, cross_read_pre)+color.END,tax_df.loc[remove_cross_asvs,'ASVg'].values)
asv_df.drop(columns=remove_cross_asvs, inplace=True)
## update stats table
tmp = asv_df.sum(axis=1).to_frame().rename(columns={0:"After sequencing run contaminant removal"}).reset_index().rename(columns={'SampleID':'sample-id'})
stats_df = stats_df.merge(tmp, on='sample-id', how='left')


""" Taq contamination -  remove after check_cross_contamination.py analysis """
contaminant_genus = ['Aeribacillus', 'Caldalkalibacillus', 'Halomonas', 'Oceanobacillus']
contaminant_unassigned = ['ASV_13']
reagent_contam_asvs = [t for t in tax_df.loc[tax_df['genus'].isin(contaminant_genus)].index if t in asv_df.columns] + [t for t in tax_df.loc[tax_df['asv'].isin(contaminant_unassigned)].index if t in asv_df.columns]
reagent_contam_read_per = asv_df[reagent_contam_asvs].sum().sum() *100 / asv_df.sum().sum()
reagent_contam_pre = asv_df[reagent_contam_asvs].astype(bool).any(axis=1).sum()*100/ asv_df.shape[0]
asv_df.drop(columns=reagent_contam_asvs, inplace=True)
#asv_df.drop(columns=[t for t in asv_df.columns if t in tax_df.loc[tax_df['family'].isin(contaminant_family)].index], inplace=True)
print(color.BOLD+color.RED+"\nRemove ASVs their abundance correlates (r > 0.90) with the most abundant genus in the negative control - reagents contamintants (Taq polymerase) - "+color.END)
print(color.BOLD+'{} ASVs were removed as reagent contaminants ({:.2f}% of total reads, detected in {:.2f}% of the samples)'.format(len(reagent_contam_asvs), reagent_contam_read_per,reagent_contam_pre)+color.END)
print([tax_df.loc[x, "ASVg"] for x in reagent_contam_asvs])


find_best_min_read_thres(asv_df, minthres=100, maxthres=1500, ExportFolder=figDir)



### export this asv table to work with the threshold for clustering...
asv_df.to_csv('/media/christos/ssd/work/Infants/tmp/non-rarefied-table.tsv', sep='\t')
tax_df.to_csv(tmpDir+'non-rarefied-taxonomy.tsv',sep='\t')


metagenomeSeq_count = convert_df_to_taxlevel_from_myClassifier(asv_df.loc[asv_df.astype(bool).sum(axis=1) > 0],tax_df,taxlevel='asv').T ## subset for samples that have at least one read and transpose table
metagenomeSeq_count = metagenomeSeq_count.loc[:,(metagenomeSeq_count!=0).any(axis=0)] ## drops columns with all zeros (OTUS)
metagenomeSeq_count.to_csv(tmpDir+'MetagenomeSeq_count.tsv', sep='\t')
metagenomeSeq_tax = tax_df.set_index('asv').reindex(metagenomeSeq_count.index)
metagenomeSeq_meta = meta_df.reindex(metagenomeSeq_count.columns)
metagenomeSeq_tax.to_csv(tmpDir+'MetagenomeSeq_tax.tsv', sep='\t')
metagenomeSeq_meta.to_csv(tmpDir+'MetagenomeSeq_meta.tsv', sep='\t')






print(color.BOLD+color.GREEN+"\nRarefying dataset to {} reads....".format(minimum_reads)+color.END)
""" Rarefraction - The function will automatically rarefy to the depth of the smallest library. """
asv_df = rarefy_dataframe(asv_df, seed=0, thres=minimum_reads)
## update stats table
tmp = asv_df.sum(axis=1).to_frame().rename(columns={0:"After rarefaction"}).reset_index().rename(columns={'SampleID':'sample-id'})
stats_df.merge(tmp, on='sample-id', how='left').to_csv('ReadsStats.csv', index=False)

### Remove ASVs that have zero counts after rarefying
print("{} ASVs with 0 reads after rarefying are removed.".format(asv_df.loc[: , asv_df.sum(axis=0) == 0].shape[1]))
asv_df = asv_df.loc[:, asv_df.sum(axis=0) > 0]

### Export data for R or qiime2 (to /tmp )
convert_df_to_taxlevel_from_myClassifier(asv_df, tax_df, taxlevel='genus').to_csv(tmpDir+'GenusCountTable.tsv',sep='\t',index_label="SampleID")
convert_df_to_taxlevel_from_myClassifier(asv_df, tax_df, taxlevel='phylum').to_csv(tmpDir+'PhylumCountTable.tsv',sep='\t',index_label="SampleID")

### Export final files - CountTable, MetaTable, TaxonomyTable, fasta with ASV sequences

Newfasta = [rec for rec in fasta if rec.id in asv_df.columns]
SeqIO.write(Newfasta, 'final_asv_seqs_clean.fasta', 'fasta') ### include only ASVs used in the downstream analysis
prepare_fasta_for_iqtree('final_asv_seqs_clean.fasta',tax_df)
asv_df = convert_df_to_taxlevel_from_myClassifier(asv_df,tax_df, taxlevel='ASVg') #### change ASVs hash to ASVg label for readability in downstream analysis !!!!!!
tax_df = tax_df.reset_index().set_index('ASVg') ## same for taxonomy
tax_df = tax_df.loc[asv_df.columns] ## include only ASVs in the CountTable
asv_df.to_csv(tmpDir+'CountTable.tsv',sep='\t', index_label='SampleID') ## use this table in other notebooks for this analysis
tax_df.to_csv(tmpDir+'TaxonomyTable.tsv', sep='\t', index_label="ASVg") ## export taxonomy table
prepare_counts_for_ecol_inference(asv_df).to_csv(tmpDir+'comm.tsv', sep='\t') ### change the ASV names for R


meta_df = meta_df.reindex(asv_df.index) ## include only samples in the CountTable
""" Add alpha-diversity indexes to the metatable  """
### Update the meta table with alpha diversities for each sample.
meta_df = meta_df.assign(AlphaShannon = div.alpha_diversity("shannon", asv_df.astype(int).values, asv_df.index, base=np.e))
meta_df = meta_df.assign(AlphaSimpson = div.alpha_diversity("simpson", asv_df.astype(int).values, asv_df.index))
meta_df = meta_df.assign(AlphaChao = div.alpha_diversity("chao1", asv_df.astype(int).values, asv_df.index))
#meta_df = meta_df.assign(Richness = asv_df.astype(bool).sum(axis=1))
meta_df = meta_df.assign(Richness = div.alpha_diversity('observed_otus', asv_df.astype(int).values,asv_df.index))
meta_df['Eveness'] = meta_df.AlphaShannon / (np.log(meta_df.Richness))
### Update the meta table that was exported earlier
meta_df.to_csv(tmpDir+'MetaTable.tsv',sep='\t', index_label='SampleID')


print(color.BOLD+color.BLUE+"\nClustering analysis...."+color.END)
find_best_number_of_clusters(asv_df, meta_df, exportFolder=figDir) #6

#### Create clusters
"""K-means clustering with 4 clusters """
master_df = clustering(asv_df, meta_df,tax_df, nclusters=4, eps=0.1, outfolder=figDir)

### Find cluster with more richness and rename it. If you divide Shannon index with ln(Richness) you get eveness.
InterCluster = master_df.groupby('Cluster').agg({'Eveness': 'mean'}).idxmax().iloc[0]
print('Define Intermediate cluster as the cluster with highest mean Eveness')
master_df.loc[master_df['Cluster'] == InterCluster, 'ClusterASV'] = 'Intermediate' ### This cluster is renamed manually !!!!
print(color.BOLD+color.BLUE+"These are the cluster labels:"+color.END)
print(master_df.drop_duplicates('Cluster')[['Cluster', 'ClusterASV']].values) ### ClusterASV is a column that holds labels for all clusters

master_df.to_csv(tmpDir+'masterTable.tsv',sep='\t')
print("{} samples and {} ASVs in the count table.".format(*asv_df.shape))


metagenomeSeq_meta = pd.read_csv(tmpDir+'MetagenomeSeq_meta.tsv',sep='\t',index_col="SampleID")
metagenomeSeq_meta = metagenomeSeq_meta.merge(master_df['ClusterASV'], on='SampleID',how='left')
metagenomeSeq_meta.to_csv(tmpDir+'MetagenomeSeq_meta.tsv', sep='\t')







