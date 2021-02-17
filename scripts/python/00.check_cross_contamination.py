import pandas as pd
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
pd.set_option('display.max_rows',None)
import re,os,sys
from functions import *

parser=argparse.ArgumentParser(description="This script will detect ASVs which originate from other datasets in the same seqeuncing run (cross-contamination) ")
parser.add_argument("-i", dest="rootDir", help="Run in the project directory - should contain the qiime2_dada2 subdirectory and the crosstTabulated.txt file.", required=True)
args=parser.parse_args()
rootDir = args.rootDir
if rootDir.endswith("/") is False:
    rootDir+= "/"

os.chdir(rootDir)
cr = pd.read_csv('crossTabulated.txt', sep='\t', index_col='OTUnumber')
tax = pd.read_csv('qiime2_dada2/taxonomy_table.tsv',sep='\t', index_col='asvid')
"""Change Escherichia taxonomy to Escherichia/Shigella"""
esc = re.compile(r'(?<=.)*Escherichia')
tax = tax.reset_index()
tax['ASVg'] =   tax.ASVg.apply(lambda x:re.sub(esc, 'Escherichia-Shigella',x))
tax['genus'] = tax.genus.apply(lambda x:re.sub(esc, 'Escherichia-Shigella',x))
tax.set_index('asvid', inplace=True)
asv = pd.read_csv('qiime2_dada2/final_asv_table.tsv',sep='\t')
renameFile = 'reseq_ids.tsv'
ren = pd.read_csv(renameFile, sep='\t')
ren.set_index('ID',inplace=True)
#ren = ren.to_dict()
#match = re.compile(r'SAMPLE\.\d+\.(\w+)')
#sampleIDs = [ren.loc[re.match(match,col).group(1),"SampleID"] for col in cr.columns if re.match(match,col)]
ms = cr.filter(regex="SAMPLE")
cr['MyMax'] = ms.max(axis=1).values
cr  = cr[["crossedDataset","highestSample", "highestSamplecounts", "MyMax", "genus_id"]]
crossGenus = cr.loc[cr['MyMax'] < cr["highestSamplecounts"],'genus_id'].unique()
#print(set([x for x in cr['genus_id'].values if  x  in cr.loc[cr['MyMax']*2 < cr["highestSamplecounts"],'genus_id'].values]))
myGenus = tax.genus.unique()
removeGenus = []
for x in myGenus:
    for y in crossGenus:
        if re.match(re.compile(r"("+x+")"+".*"), y):
            if x not in removeGenus:
                removeGenus.append(x)

print("Genera that have higher read count in samples from other datasets in the run")
print(len(pd.Series(removeGenus).unique()))
print(pd.Series(removeGenus).unique())

### check if any of those genera is detected in the first run as well - if so remove it from list

ren = pd.read_csv(renameFile, sep='\t')
replace_IDs = ['S113', 'S129', 'S131', 'S137', 'S145', 'S153', 'S155', 'S160', 'S162', 'S163', 'S164', 'S165', 'S171', 'S173', 'S174', 'S175', 'S177', 'S178', 'S179', 'S27', 'S42', 'S9', 'S91', 'S99',"S59"]
ren.set_index('ID',inplace=True)
ren = ren.loc[replace_IDs,] ### these are the samples we included from resequencing
asv =asv.loc[~asv['SampleID'].isin(ren.SampleID.values)] ## get samples from first run
asv.set_index("SampleID",inplace=True)
asv.drop(index=['neg2','neg3'],inplace=True)
asv = convert_df_to_taxlevel_from_myClassifier(asv,tax,taxlevel='genus')
asvmax = asv.max(axis=0).to_frame().reset_index().rename(columns={'index':'genus', 0:'FirstRunMax'})
asvpre = asv.astype(bool).sum(axis=0).to_frame().reset_index().rename(columns={'index':'genus',0:'FirstRunPrevalence'})
asvpre['FirstRunPrevalence'] = (asvpre['FirstRunPrevalence'] / asv.shape[0]).apply(lambda x:round(x,2))


#print(asvmax.loc[asvmax['genus'].isin(removeGenus)].sort_values('FirstRunMax', ascending=False))
cr.rename(columns={'genus_id':'genus'}, inplace=True)
cr = cr.merge(asvmax,on='genus',how='outer')
cr = cr.merge(asvpre,on='genus',how='outer')

#print(cr.sort_values(["MyMax","highestSamplecounts"],ascending=False))
cr.sort_values(["MyMax","highestSamplecounts"],ascending=False).to_excel("/home/christos/Desktop/cross.xlsx")

asv = pd.read_csv('/media/christos/ssd/work/Infants/publication/qiime2_dada2/final_asv_table.tsv',sep='\t',index_col="SampleID")
asv = convert_df_to_taxlevel_from_myClassifier(asv,tax,taxlevel='genus')
"""Select potential contaminants - they appear only in the reseq run and for those that have a genus taxonomy there are samples 
from other datasets with much higher read count. (After manually inscpecting cross.xlsx file..)"""
#pc = asv[["Faecalibacterium",'Cetobacterium','ASV_124',"ASV_59","ASV_81",'ASV_92','ASV_93','ASV_94','Planoglabratella']]
potential_contaminants = cr.loc[(cr["highestSamplecounts"] > cr["MyMax"]) & (cr.FirstRunMax == 0)]['genus'].unique().tolist()
print(color.BOLD+color.CYAN+"Use this list in the 2.process_asv_data.py for remove_cross_genus variable!")
print(potential_contaminants,color.END) ###
pc = asv[potential_contaminants]
pd.set_option('display.max_columns',None)
print(pc.agg({'max','sum','count_nonzero'},axis=0))

potential_contaminants_unassigned = cr.loc[(cr["MyMax"].isna()) & (cr['FirstRunMax'] == 0)]['genus'].unique().tolist() ### those are asvs that are not part of the crosstab table and also not detected in the first run
pc_un = asv[potential_contaminants_unassigned]
print(color.BOLD+color.CYAN+"Use this list in the 2.process_asv_data.py for remove_cross_unassigned variable!")
print(potential_contaminants_unassigned, color.END) ### use this list in the 2.process_asv_data.py for remove_cross_unassigned variable
print(pc_un.agg({'max','sum','count_nonzero'},axis=0))



##### check how Taq contaminants corelate with each other and use this as an indicator for other potential contaminants.

from scipy.stats import pearsonr
import itertools
print("\nCheck for Taq polymerase known contaminants")
contaminant_genus = ['Halomonas','Chromohalobacter', 'Caldalkalibacillus', 'Geobacillus', 'Aeribacillus']
cor = pd.DataFrame()
"""Find if reagent contaminants correlate to each other..."""

### All but Geobacillus not correlated - maybe due to low read numbers
for x,y in itertools.combinations(contaminant_genus, 2):
    r,p = pearsonr(asv[x],asv[y])
    tf = pd.DataFrame({"x":[x],'y':[y],'r':[r], 'p':[p]})
    cor = cor.append(tf)
print(cor.sort_values('r',ascending=False))
#print(asv['Geobacillus'].sum(), round(asv['Geobacillus'].astype(bool).sum()/asv.shape[0], 2))


"""Check if Caldalkalibacillus that is the most abundant genus in the negative controls correlates with any of the other genera in the dataset."""
cor_halo = pd.DataFrame()
import seaborn as sns
for c in asv.columns:
    r,p = pearsonr(asv["Caldalkalibacillus"],asv[c])
    pre = asv[c].astype(bool).sum() *100 / asv.shape[0]
    tf = pd.DataFrame({"x":["Caldalkalibacillus"],'y':[c],'r':[r], 'p':[p]})
    cor_halo = cor_halo.append(tf)
cor_halo = cor_halo.melt(id_vars=['x','y'], value_vars='r').rename(columns={'y':'genus', 'value':'r'})
cor_halo['contaminant'] = ["yes" if x > 0.90 else 'no' for x in cor_halo.r]
cor_halo.to_csv('cor_calda.csv')
#print(cor_halo.loc[cor_halo.r > 0.7].sort_values('r',ascending=False)['genus'].values)
print(cor_halo.sort_values('r',ascending=False))
print("\nAssume genera for with pearson R > .9 for Caldalkalibacillus (most abundnant in negative and also reported as Taq polymerase contaminant) are coming from the same source.")
print(color.BOLD+color.CYAN+"Use this list as reagent contaminant genera to remove asvs (contaminant_genus and contaminant_unassigned variables)")
print(cor_halo.loc[cor_halo['r'].ge(0.9),'genus'].unique().tolist(), color.END)
