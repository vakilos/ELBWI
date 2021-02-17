from functions import *
import scipy
import sys, os, re, argparse, subprocess, itertools
import numpy as np
import pandas as pd

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)


parser=argparse.ArgumentParser(description="This script run most of the analysis and export datasets for plotting Figures")
parser.add_argument("-i", dest="rootDir", help="Project directory - (should contain tmp", required=True)
args=parser.parse_args()

"""Import files that have been processed with asv_analysis.ipby)"""
rootDir = args.rootDir
if rootDir.endswith("/") is False:
    rootDir += '/'

tmpDir = rootDir+'tmp/'

os.chdir(tmpDir)


asv = pd.read_csv('CountTable.tsv',sep='\t', index_col='SampleID')
tax = pd.read_csv('TaxonomyTable.tsv', sep='\t', index_col='ASVg')
meta = pd.read_csv('MetaTable.tsv',sep='\t', index_col='SampleID')
master = pd.read_csv('masterTable.tsv', sep='\t',index_col='SampleID')
clinc = pd.read_csv('meta_with_clinical.tsv', sep='\t', index_col="SampleID")


### Use this table for ecological processes inference...
clinc = clinc[meta.columns.tolist() + ['oxygenCumIndex','AntiCumIndex', 'il6WinIndex']].rename(columns={'AntiCumIndex':'Antibiotics','oxygenCumIndex':'Oxygenation', "il6WinIndex":'IL-6'}).to_csv('/media/christos/ssd/work/Infants/tmp/metaTable_clinical.tsv', sep='\t')


#### FIGURE 1a - genera abundances and prevalence
print('\nFig. 1A')
"""Genera with > 0.2 relative abundance in at least one sample are shown per bodysite """

most_abundant_per_bodysite(asv,tax,meta, Taxlevel='genus',filter=0.2, exportDir=tmpDir)

""" Species prevelance """
species_prevalence(asv,meta,tax, Taxlevel='genus', exportDir=tmpDir)
"""Outputs  /tmp/mostabuperbodysitegenus.tsv and /tmp/prev_toR.tsv are parsed in R ( Fig1A_abundance_prevelance.R )"""

### make one function for abundance and prevalence.

#### FIGURE 1b - phyla abundances
print('\nFig. 1B')
most_abundant_per_bodysite(asv,tax,meta, Taxlevel='phylum',filter=0.2, exportDir=tmpDir)

"""Outputs  /tmp/mostabuperbodysitephylum.tsv to be parsed in R (Fig1BC_plot_most_abundant_taxa_per_bodysite.R ) """


### FIGURE 1c - envfit
## command in R -- > envfit(ord=mds.values , env=M, data =M , perm=10000) (envfit.R)
### Save envfit results in a a file --> envFit_to_R.tsv 

"""'envFit_to_R.tsv' to be used with envfit_final_plot.R"""





### FIGURE 2a
""" This is plotted by plot_alpha_diversities.R """


### FIGURE 2b - inter-intra body site distance per patient!
print('\nFig. 2B')
from scipy.spatial.distance import pdist, squareform
ddf = pd.DataFrame(squareform(pdist(asv, metric='braycurtis')), index=asv.index, columns=asv.index)
bd = pd.DataFrame()
for name, group in meta.reindex(asv.index).groupby(['Patient']):
    for sample1, sample2 in itertools.combinations(group.index, 2):
        if group.loc[sample1,"SampleType"] == group.loc[sample2, 'SampleType']:
            #comp = group.loc[sample1, 'SampleType']
            comp = 'intra'
        else:
            comp = 'inter'
        bd = bd.append(pd.DataFrame({"Comparison":[comp], 'Patient':[name], 'Distance':[ddf.loc[sample1, sample2]]}), ignore_index=True)

bd.groupby(['Patient','Comparison']).agg({'Distance':'mean'})
pvalue_dict = {}
for name, group in bd.groupby(['Patient']):
    #print(group.Comparison.shape[0], group.Comparison ==)
    stat, pvalue = scipy.stats.mannwhitneyu(group.loc[group.Comparison == 'inter', 'Distance'], group.loc[group.Comparison == 'intra','Distance'], alternative='two-sided')
    pvalue_dict[name] = pvalue
bd['pvalue'] = bd.Patient.map(pvalue_dict).round(5)
bd['asterisk'] = bd.pvalue.apply(lambda x: '*' if (x > 0.01 and x <= 0.05) else("**" if (x <= 0.01 and x > 0.001) else ("***" if x <= 0.001 else "")) )
bd.to_csv('DistanceBodysites_toR.tsv', index=False, sep='\t')
print(bd.drop_duplicates('Patient').pvalue)
print(bd.loc[bd.asterisk != ''].drop_duplicates("Patient"))

""" Output is used with plot_BodysiteDistanceComparison.R """



### FIGURE 3a + 3b
print("Fig 3A, 3B")
""" masterTable.tsv nd DominantASVperCluster.tsv are used with plot_DominantASVPerCluster_fig3.R"""
print("\n Most Abundant ASVs per cluster.")
print(master.ClusterASV.unique())

### Relative abundance of the most dominants taxa per cluster
r = relativeAbundance_from_dataframe(asv)
r = r[[x for x in master.ClusterASV.unique() if x != 'Intermediate']]
r = r.merge(master[['Cluster','ClusterASV']], on='SampleID', how='inner')
m = r.melt(id_vars=['Cluster','ClusterASV'])
m.to_csv('DominantASVperCluster.tsv',sep='\t', index=False)

#print( r.groupby('ClusterASV').mean())
import numpy as np, scipy.stats as st
for n, g in r.groupby('ClusterASV'):
    if n != 'Intermediate':
        print("{0} -  mean: {1}, CI: ({2},{3}) ".format(n, g[n].mean(), *st.t.interval(0.95, len(g[n])-1, loc=np.mean(g[n]), scale=st.sem(g[n]))))
    else:
        print("Intermediate:")
        for x in [x for x in master.ClusterASV.unique() if x != 'Intermediate']: ## print abundance of all pre dominant ASVs for Intermediate cluster
            print("\t{0} -  mean: {1}, CI: ({2},{3}) ".format(x, g[x].mean(), *st.t.interval(0.95, len(g[x]) - 1, loc=np.mean(g[x]), scale=st.sem(g[x]))))

### FIGURE 3c - Richness across clusters

"""masterTable.tsv is used with Richness_for_Clusters_fig3.R """


### FIGURE 3d - Heatmap for clusters and body sites
print("\nFig. 3D")
CBS = master.groupby(['SampleType','ClusterASV']).size().reset_index().rename(columns={0:'value'})
CDict = master.groupby('ClusterASV').size().to_dict()
CBS['ClusterSamples'] = CBS['ClusterASV'].map(CDict)
CBS.value = CBS.value / CBS.ClusterSamples
CBS.ClusterASV.replace({'ASV_208:Escherichia/Shigella':'EC',
                        'Intermediate':"IC",
                        "ASV_184:Lactobacillus":"LC",
                        "ASV_66:Staphylococcus":"SC"}, inplace=True)
CBS.SampleType = CBS.SampleType.str.replace('_'," ")
CBS.to_csv('ClustersInBodySites.csv', index=False)

"""ClustersInBodySites.csv is used with ClusterInBodysiteHeatmap_fig3.R  """


### FIGURE ABC- ecological processes inference



### FIGURE 4d - Posterior probability for detection in multiple body sites
"""Here i estimate prevalence as a percentage of samples in a particular body site. Also create files for plotting heatmaps in R."""
print('\nFig. 4D')
###Prevalence, this function estimates global prevalence and not body site specific
sp = asv.astype(bool).sum(axis=0).to_frame().rename(columns={0:'per'}).div(asv.shape[0]).reset_index().rename(columns={'index':'ASV'})
bayes = bayes_per_bodysite(asv,meta, tax, Taxlevel=None, exportDir=tmpDir)
print(bayes.groupby('Posterior').agg({'value':'mean'}))

#shared = prob_of_sharing(asv, meta)
odds = odds_ratio(asv,meta)
print(odds)
#print(shared.groupby('Sites').agg({'value':'mean'}))

ssp = master.groupby('SampleType')[asv.columns].agg(lambda x: x.astype(bool).sum(axis=0)).T
ssp.Gut = ssp.Gut / master.loc[master.SampleType == 'Gut'].shape[0]
ssp.Oral_cavity = ssp["Oral_cavity"] / master.loc[master.SampleType == 'Oral_cavity'].shape[0]
ssp.Skin = ssp.Skin / master.loc[master.SampleType == 'Skin'].shape[0]
ssp.index.rename("ASV", inplace=True)
pp = ssp.merge(sp, on='ASV', how='right').merge(bayes, on ='ASV', how ='right')#.to_csv('../tmp/PosteriorPrevalence.csv', index=False)
#pp = ssp.merge(sp, on='ASV', how='right').merge(shared, on ='ASV', how ='right')#.to_csv('../tmp/PosteriorPrevalence.csv', index=False)
pp = pp.pivot(columns='Posterior', values='value', index='ASV').fillna(0)
#pp = pp.pivot(columns='Sites',values='value',index="ASV").fillna(0)
pp.to_csv('PosteriorHeatmap.csv')
ssp.rename(columns={'Oral_cavity': "Oral cavity"}).reindex(pp.index).to_csv('BodysitePrevalenceAsPercentage.csv')





### Fig 5a - Schematic represantation of transition probabilities
print('\nFig. 5A')

transitions = transition_prob_from_master(master) 

transitions.to_csv('Transitions.csv',index=False)



### Figure 5b - Steady state frequencies of clusters.
"""We manually imported transition probabilities in the markov.R script to estimate the steady state frequencies.
This script also plots the respective plot."""

### Figure 5c - Mono-dominated cluster frequency in IVH groups
print('\nFig. 5C')
af = master[['Patient','IVH', "BPD", "SampleType", 'ClusterASV']]

labels = {y:"Monodominated" for y in [x for x in af.ClusterASV.unique() if x != 'Intermediate']}
af.ClusterASV.replace(labels, inplace=True)
af.IVH.replace({1:'yes',2:'yes',0:'no'},inplace=True)
af.BPD.replace({1:'no',2:'yes',0:'no'},inplace=True)
af.SampleType = af.SampleType.str.replace("_", " ")
bf = af
af = af.groupby(['IVH', 'SampleType','Patient']).agg({'ClusterASV':'value_counts'}).unstack().fillna(0).reset_index()
bf = bf.groupby(['BPD','SampleType',"Patient"]).agg({"ClusterASV":'value_counts'}).unstack().fillna(0).reset_index()
print(af)
af['MonoFreq'] = af[('ClusterASV','Monodominated')] / (af[('ClusterASV','Monodominated')] + af[('ClusterASV','Intermediate')])
bf["MonoFreq"] = bf[("ClusterASV","Monodominated")]  / bf[("ClusterASV","Monodominated")] + bf[("ClusterASV", "Intermediate")]

af.columns = [x[0]+x[1] for x in af.columns]
bf.columns = [x[0]+x[1] for x in bf.columns]
af.rename(columns={'ClusterASVIntermediate':'Intermediate', 'ClusterASVMonodominated':'Monodominated'}, inplace=True)
bf.rename(columns={'ClusterASVIntermediate':'Intermediate', 'ClusterASVMonodominated':'Monodominated'}, inplace=True)
bf.to_csv("MonoFreq_bpd_toR.csv",index=False)
af.to_csv('MonoFreq_toR.csv', index=False)
"""MonoFreq_toR.csv is then imported in cluster_risk.R for plotting."""










