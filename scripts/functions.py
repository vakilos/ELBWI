import sys, os, re, argparse, subprocess, itertools
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import skbio.diversity as div
from skbio.stats.ordination import pcoa
import datetime as dt
from Bio import SeqIO
import scipy
#import ecopy as ep
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
"""Check the Read numbers for different minimum thresholds.."""
def find_best_min_read_thres(df, minthres=100, maxthres=10000, ExportFolder=None):
    samples_left = [df.loc[df.sum(axis=1) > threshold].shape[0] for threshold in range(minthres, maxthres,100)]
    #dataset_otus = [asv_df.loc[asv_df.sum(axis=1) > threshold , asv_df.sum(axis=0) > 0].shape[1] for threshold in range(100,10000,100)]   ##(df != 0).any(axis=0) returns a boolean Series indicating which columns have nonzero entries. (The any operation aggregates values along the 0-axis -- i.e. along the rows -- into a single boolean value. Hence the result is one boolean value for each column.)
    dataset_otus = []
    for threshold in range(minthres, maxthres, 100):
        tdf = df.loc[df.sum(axis=1) > threshold]
        tdf = tdf.loc[:, tdf.sum(axis=0) > 0]
        dataset_otus.append(tdf.shape[1])
    fig,ax1 = plt.subplots()
    ax1.plot(range(minthres, maxthres, 100), dataset_otus, color='blue')
    ax1.set(xlabel='Minimum Read Number')
    ax1.set_ylabel('Number of ASVs', fontsize=14, color='blue')
    plt.title("ASVs and samples left in the dataset as function of number of reads")
    #ax1.set_ylabel("Drug administered volume (ml)", fontsize=16, color='blue')
    #plt.xticks([0,30,60,90,final_time], [str(x) for x in [0,30,60,90,final_time]])
    ax1.tick_params(axis='y', labelcolor='blue')
    ax2 = ax1.twinx()
    ax2.plot(range(minthres, maxthres, 100), samples_left, color='green')
    ax2.set_ylabel("Number of samples", fontsize=14, color='green')
    ax2.tick_params(axis='y', labelcolor='green')
    if ExportFolder is None:
        fig.savefig('/media/christos/ssd/work/Infants/figures/1.Number-of-samples-ASVs-for-different-read-count-thresholds.png')
    else:
        fig.savefig(ExportFolder+"/Number_of_samples_and_ASVS_for_min_read_threshold.png")
    #plt.show()

def convert_df_to_taxlevel_from_myClassifier(CountsDf, TaxDf, taxlevel="phylum"):
    """"CountsDf should have OTU names in columns"""
    #print("\nConverting OTU counts for {} level".format(taxlevel))
    Df = CountsDf.copy()
    tax_levels = TaxDf.columns
    Df.rename(columns={k:v for k,v  in zip(Df.columns, [TaxDf.loc[x, taxlevel] if TaxDf.loc[x, taxlevel] !='Unassigned' else TaxDf.loc[x, 'asv'] for x in Df.columns])}, inplace=True)
    Df = Df.groupby(Df.columns, axis=1).sum()
    return Df
def relativeAbundance_from_dataframe(countDf):
    "Count dataframe should have samples in rows and OTUs in columns"
    return countDf.fillna(0).apply(lambda x: x/float(np.sum(x)), axis=1)

def rarefy_dataframe(Df, thres=200, seed=0):
    from numpy.random import RandomState
    prng =RandomState(seed)
    if thres == "min":
        M = Df.values
        noccur = np.sum(M, axis=1)
        nvar = M.shape[1]
        depth = int(np.min(noccur)) ### find sample with minimum reads to rarefy
        old_samples = Df.shape[0]
        Df = Df[Df.sum(axis=1).ge(depth)]  ## remove samples below threshold
        print('\n{} samples were removed after rarefying'.format(old_samples - Df.shape[0]))
    else:
        depth = thres
        old_samples = Df.shape[0]
        Df = Df[Df.sum(axis=1).ge(depth)]  ## remove samples below threshold
        print('\n{} samples were removed after rarefying'.format(old_samples - Df.shape[0]))
        M = Df.values
        noccur = np.sum(M, axis=1)
        nvar = M.shape[1]

    Mrarefied = np.empty_like(M)
    for i in range(M.shape[0]): # for each sample
        p = M[i] / float(noccur[i]) # relative frequency / probability
        choice = prng.choice(nvar, depth, p=p)
        Mrarefied[i] = np.bincount(choice, minlength=nvar)
    return pd.DataFrame(Mrarefied, columns=Df.columns, index=Df.index)

def most_abundant_per_bodysite(CountsDf,TaxDf, MetaDf, Taxlevel='phylum', filter=0.2, exportDir=''):
    Df = relativeAbundance_from_dataframe(convert_df_to_taxlevel_from_myClassifier(CountsDf, TaxDf, taxlevel=Taxlevel))
    MasterDf = Df.merge(MetaDf, on='SampleID', how='inner')
    melted = MasterDf.melt(id_vars=['SampleType'], var_name='Taxonomy', value_name='Relative_abundance', value_vars=Df.columns)
    most_abundant_taxa = melted.loc[melted.Relative_abundance >= filter].Taxonomy.unique()
    melted = melted.loc[melted.Taxonomy.isin(most_abundant_taxa)]
    #prev = MasterDf.groupby("SampleType").agg({Df.columns:'count_nonzero'})#.T.reset_index().rename(columns={'index':"Taxonomy","count_nonzero":"Prevalence"})
    #melted = melted.merge(prev, on='Taxonomy', how='left')
    melted.SampleType.replace({'Oral_cavity':'Oral cavity'}, inplace=True)
    categories = melted.groupby(['Taxonomy','SampleType']).agg('mean').reset_index().groupby(['Taxonomy']).agg({"Relative_abundance":'max'}).reset_index().sort_values("Relative_abundance",ascending=False).Taxonomy.values.tolist()
    print(categories)
    melted.Taxonomy = pd.Categorical(melted.Taxonomy, categories=categories, ordered=True)
    melted.sort_values('Taxonomy',ascending=False,inplace=True)
    os.chdir(exportDir)
    melted.to_csv('mostabuperbodysite'+Taxlevel+'.tsv',sep='\t', index=False)
    return melted

def group_non_abundant_taxa(CountsDf, taxdf, threshold=0.05):
    """Groups together OTUs or (whatever the columns are in the df) that are detected
    below a certain relative abundance"""
    CountsDf["Other ( < "+str(threshold)+")"] = [0 for x in CountsDf.index]
    #print(CountsDf.sum(axis=1))
    #CountsDf.fillna(0, inplace=True)
    for index, row in CountsDf.iterrows():
        CountsDf.loc[index,"Other ( < "+str(threshold)+")"] = row.where(row.le(threshold, fill_value=0)).sum(axis=0)

    mask = (CountsDf["Other ( < "+str(threshold)+")"] > 0)
    maskDf = CountsDf[mask]
    CountsDf[maskDf < threshold] = 0
    #print(CountsDf["Other ( < "+str(threshold)+")"].to_list() )
    CountsDf = CountsDf.apply(lambda x: x/float(np.sum(x)), axis=1)
    print("Grouping taxa with relative abundance less than {}!".format(threshold))
    return CountsDf

def barplot_rel_abundance_simple(CountsDf, outfilename="relAbuggplot.svg"):
    plt.close()
    Df = CountsDf.copy()
    Df = relativeAbundance_from_dataframe(Df)
    Df.reset_index(inplace=True)
    Dfmelt = pd.melt(Df, id_vars=["SampleID"], value_vars=list(Df.columns).remove('SampleID'),value_name="Relative_abundance", var_name="Taxonomy")
    Dfmelt = Dfmelt.loc[Dfmelt["Relative_abundance"] > 0]
    Dfmelt.to_csv("/media/christos/ssd/work/Infants/tmp/RelAbuNozeros.tsv", sep='\t', index=False)
    subprocess.run('Rscript /media/christos/ssd/work/Infants/scripts/plot_relative_abundance_simple.R /media/christos/ssd/work/Infants/tmp/RelAbuNozeros.tsv '+outfilename, shell=True)
    os.remove("/media/christos/ssd/work/Infants/tmp/RelAbuNozeros.tsv")
def barplot_rel_abundance_patient(CountsDf, TaxDf, MetaDf, taxlevel='genus',  outfilename='relAbuPatient.png'):
    plt.close()
    Df = CountsDf.copy()
    Df = convert_df_to_taxlevel_from_myClassifier(Df, TaxDf, taxlevel='genus')

    Df = relativeAbundance_from_dataframe(Df)
    Df = group_non_abundant_taxa(Df, TaxDf, threshold=0.05)
    MasterDf = Df.merge(MetaDf, on='SampleID', how='inner')
    MasterDf.reset_index(inplace=True)
    Dfmelt = pd.melt(MasterDf, id_vars=list(set(MasterDf.columns)-set(Df.columns)), value_vars=Df.columns, value_name="Relative_abundance", var_name="Taxonomy")
    Dfmelt = Dfmelt.loc[Dfmelt["Relative_abundance"] > 0]
    Dfmelt.to_csv("/media/christos/ssd/work/Infants/tmp/RelAbuNozeros.tsv", sep='\t', index=False)
    subprocess.run('Rscript plot_relative_abundance_patient.R /media/christos/ssd/work/Infants/tmp/RelAbuNozeros.tsv '+outfilename, shell=True)
    #os.remove("/media/christos/ssd/work/Infants/tmp/RelAbuNozeros.tsv")

def make_colors_for_otus(N, seed=0):
    import random
    import pandas as pd
    random.seed(seed)
    colors = ["#%06x" % random.randint(i, 0xFFFFFF) for i in range(N)]
    C = pd.DataFrame({'otu_color':colors})
    return C

def plot_neotype_microbiome(CountsDf, MasterDf, TaxDf, clusterCol='Cluster', Taxlevel=None, outfilename='NeotypeMicrobiome.png'):
    Df = CountsDf.copy()
    if Taxlevel is not None:
        Df = convert_df_to_taxlevel_from_myClassifier(Df, TaxDf, taxlevel=Taxlevel)
    taxaCol = Df.columns
    make_colors_for_otus(Df.shape[1], seed=0).to_csv('/media/christos/ssd/work/Infants/tmp/neotype_colors.tsv',sep='\t', index=False)
    Df.index.equals(MasterDf.index)
    Df.index = Df.index.map(dict(zip(Df.index, MasterDf[clusterCol])))
    Df = Df.groupby(Df.index).agg('sum')

    Df = relativeAbundance_from_dataframe(Df)
    #print(Df.idxmax(axis=1))
    Df.index.names = [clusterCol]
    Df.reset_index(inplace=True)
    Dfmelt = pd.melt(Df, id_vars=clusterCol, value_vars=taxaCol,
                     value_name="Relative_abundance", var_name="Taxonomy")
    Dfmelt = Dfmelt.loc[Dfmelt["Relative_abundance"] > 0]
    Dfmelt.to_csv("/media/christos/ssd/work/Infants/tmp/RelAbuNeotype.tsv", sep='\t', index=False)
    subprocess.run('Rscript plot_neotype_microbiome.R /media/christos/ssd/work/Infants/tmp/RelAbuNeotype.tsv ' +outfilename+" "+clusterCol,shell=True)

def areaplot_rel_abundance_patient_over_time(CountsDf, TaxDf, MetaDf, taxlevel='genus',outfilename='relAbuPatientOvertime.png'):
    Df = CountsDf.copy()
    Df = convert_df_to_taxlevel_from_myClassifier(Df,TaxDf, taxlevel=taxlevel)
    Df = relativeAbundance_from_dataframe(Df)
    Df = group_non_abundant_taxa(Df, TaxDf, threshold=0.05)

    MasterDf = Df.merge(MetaDf, on='SampleID', how='inner')
    MasterDf.reset_index(inplace=True)
    Dfmelt = pd.melt(MasterDf, id_vars=list(set(MasterDf.columns) - set(Df.columns)), value_vars=Df.columns,
                     value_name="Relative_abundance", var_name="Taxonomy")
    #Dfmelt = Dfmelt.loc[Dfmelt["Relative_abundance"] > 0]
    Dfmelt.to_csv("/media/christos/ssd/work/Infants/tmp/frameRelAbu.tsv", sep='\t', index=False)

    subprocess.run('Rscript plot_areaplot_relativeAbundance_patient_overtime.R /media/christos/ssd/work/Infants/tmp/frameRelAbu.tsv ' + outfilename,
        shell=True)
    os.remove('/media/christos/ssd/work/Infants/tmp/frameRelAbu.tsv')

def pcoa_plot_bodysite(CountsDf, MetaDf,  metric='braycurtis', outfilename='PCoA_bodysite'):
    plt.close('all')
    import scipy
    from skbio.stats.distance import DistanceMatrix
    from skbio.stats.ordination import pcoa
    Distance_matrix = DistanceMatrix(scipy.spatial.distance.pdist(CountsDf.values, metric=metric)).data
    Df = pd.DataFrame(Distance_matrix, index=CountsDf.index, columns=CountsDf.index)
    My_pcoa = pcoa(Df)
    pcoaDf = My_pcoa.samples[["PC1", "PC2"]]
    pcoaDf = pcoaDf.assign(SampleID=Df.index)
    pcoaDf.set_index('SampleID', inplace=True)
    pc1_prop_explained, pc2_prop_explained = My_pcoa.proportion_explained.loc["PC1"], My_pcoa.proportion_explained.loc["PC2"]
    pcoaDf = pcoaDf.merge(MetaDf, on='SampleID', how='inner')

    stat, pvalue = scipy.stats.kruskal(pcoaDf.loc[pcoaDf['SampleType'] == 'Gut', 'PC1'],
                        pcoaDf.loc[pcoaDf['SampleType'] == 'Skin', 'PC1'],
                        pcoaDf.loc[pcoaDf['SampleType']== 'Oral_cavity','PC1'])
    stat2, pvalue2 = scipy.stats.kruskal(pcoaDf.loc[pcoaDf['SampleType'] == 'Gut', 'PC2'],
                        pcoaDf.loc[pcoaDf['SampleType'] == 'Skin', 'PC2'],
                        pcoaDf.loc[pcoaDf['SampleType']== 'Oral_cavity','PC2'])
    print("Kruskal-Wallis H-test on seperation in the PC1 axis, stat={}, p-value={} .".format(stat,pvalue))
    print("Kruskal-Wallis H-test on seperation in the PC2 axis, stat={}, p-value={} .".format(stat2,pvalue2))
    pcoaDf.to_csv('/media/christos/ssd/work/Infants/tmp/pcoa_bodysite_toR.tsv', sep='\t')
    pal = {'Gut':sns.xkcd_rgb["brown"], 'Skin':sns.xkcd_rgb["peach"], 'Oral_cavity':sns.xkcd_rgb["salmon"]}
    #fig, ax = plt.subplots(figsize=(10,8))
    #ax = sns.scatterplot(x="PC1",y='PC2', hue="SampleType", data=pcoaDf, style='Timepoint', s=90, alpha=.8, palette=pal)
    #ax.set(xlabel="PC1 ({0:.2f}%)".format(pc1_prop_explained*100), ylabel="PC2 ({0:.2f}%)".format(pc2_prop_explained*100) )
    #plt.title("PCoA on Bray-Curtis Dissimilarities")

    #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    #plt.tight_layout()
    #plt.show()
    #fig.savefig(outfilename+'.png')
    #fig.savefig(outfilename + '.svg')
    #plt.close()
def pcoa_plot_weighted_unifrac(CountsDf, MetaDf, TaxDf, Treefile='/media/christos/ssd/work/Infants/qiime2/rooted_tree/tree.nwk'):
    plt.close()
    from skbio import TreeNode
    tree = TreeNode.read(Treefile)
    wu = div.beta_diversity('weighted_unifrac', CountsDf.astype(int).values, CountsDf.index, otu_ids=TaxDf.loc[CountsDf.index,'asvid'].values, tree=tree)
    pc = pcoa(wu)
    pcDf = pc.samples[['PC1','PC2']]
    pcDf = pcDf.assign(SampleID=CountsDf.index)
    pcDf.set_index('SampleID', inplace=True)
    pc1_prop_explained, pc2_prop_explained = pc.proportion_explained.loc['PC1'], pc.proportion_explained.loc['PC2']
    pcDf = pcDf.merge(MetaDf, on='SampleID', how='inner')
    pal = {'Gut':sns.xkcd_rgb["brown"], 'Skin':sns.xkcd_rgb["peach"], 'Oral_cavity':sns.xkcd_rgb["salmon"]}
    fig, ax = plt.subplots(figsize=(10,8))
    ax = sns.scatterplot(x="PC1",y='PC2', hue="SampleType", data=pcDf, style='Timepoints', s=90, alpha=.8, palette=pal)
    ax.set(xlabel="PC1 ({0:.2f}%)".format(pc1_prop_explained*100), ylabel="PC2 ({0:.2f}%)".format(pc2_prop_explained*100) )
    plt.title("PCoA on Weighted-Unifrac distance matrix")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    #plt.show()
    #fig.savefig(outfilename)

    return wu
def pcoa_plot_per_bodysite_patient(CountsDf, MetaDf, TaxDf, Taxlevel='asv', metric='braycurtis', outfilename='/media/christos/ssd/work/Infants/figures/PCoA_bodysite_patient.png'):
    Df = CountsDf.copy()
    Df = convert_df_to_taxlevel_from_myClassifier(Df,TaxDf, taxlevel=Taxlevel)
    #Df = Df.merge(MetaDf, on='SampleID', how='inner')
    Df.to_csv('/media/christos/ssd/work/Infants/tmp/Counts_topcoaR.tsv',index=True,sep='\t')
    subprocess.run("Rscript /media/christos/ssd/work/Infants/scripts/plot_pcoa.R "+outfilename,shell=True)

def pcoa_plot_date(MasterDf, CountsDf,  metric='braycurtis', outfilename='PCoA_bodysite.png'):
    """Meta Df"""
    plt.close('all')
    import scipy
    from skbio.stats.distance import DistanceMatrix
    from skbio.stats.ordination import pcoa
    Distance_matrix = DistanceMatrix(scipy.spatial.distance.pdist(CountsDf.values, metric=metric)).data
    Df = pd.DataFrame(Distance_matrix, index=CountsDf.index, columns=CountsDf.index)
    My_pcoa = pcoa(Df)
    pcoaDf = My_pcoa.samples[["PC1", "PC2"]]
    pcoaDf = pcoaDf.assign(SampleID=Df.index)
    pcoaDf.set_index('SampleID', inplace=True)
    pc1_prop_explained, pc2_prop_explained = My_pcoa.proportion_explained.loc["PC1"], My_pcoa.proportion_explained.loc["PC2"]
    pcoaDf = pcoaDf.merge(MasterDf, on='SampleID', how='inner')
    sns.set()
    fig, ax = plt.subplots(figsize=(10,8))
    cmap = sns.diverging_palette(10,220, s=90, n=pcoaDf.drop_duplicates('Time').shape[0])
    #cmap = sns.cubehelix_palette(rot=-.2, as_cmap=True)
    sns.scatterplot(x="PC1",y='PC2', hue="Time", data=pcoaDf, style='Cluster', s=90, ax=ax, palette=cmap, edgecolor='k')
    ax.set(xlabel="PC1 ({0:.2f}%)".format(pc1_prop_explained*100), ylabel="PC2 ({0:.2f}%)".format(pc2_prop_explained*100) )
    plt.title("PCoA on Bray-Curtis Dissimilarities")
    h,l = ax.get_legend_handles_labels()
    col_lgd = plt.legend(h[:-4:5], l[:-4:5], loc='lower center',
                     bbox_to_anchor=(0.5, -0.45), fancybox=True, shadow=True, ncol=4)
    #import matplotlib
    #cax, _ = matplotlib.colorbar.make_axes(ax)
    #normalize = matplotlib.colors.Normalize(vmin=pcoaDf[pcoaDf["SampleType"]=='Stool']['Time'].min())
    #cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=normalize)
    plt.gca().add_artist(col_lgd)
    plt.tight_layout()
    fig.savefig(outfilename)


def count_ASVs_per_genus(CountsDf, GenusCountsDf, TaxDf):
    """Returns a dataframe that has index the columns of the genus-counts dataframe and a count column for the count of ASVs that fall into this taxonomy"""
    NewDf = pd.DataFrame({"Genus":[], 'ASVCounts':[]})
    for genus in GenusCountsDf.columns:
        counts = len([col for col in CountsDf.columns if TaxDf.loc[col, 'genus'] == genus])
        NewDf = NewDf.append(pd.DataFrame({"Genus":[genus], "ASVCounts":[counts]}))
    NewDf.set_index('Genus', inplace=True)
    NewDf.sort_values('ASVCounts',ascending=False, inplace=True)
    return  NewDf

def pcoa_from_dataframe(CountsDf , metric='braycurtis'):
    import scipy
    from skbio.stats.distance import DistanceMatrix
    from skbio.stats.ordination import pcoa
    My_pcoa = pcoa(scipy.spatial.distance.pdist(CountsDf.values, metric=metric))
    PcoaDf = My_pcoa.samples[["PC1","PC2"]].set_index(CountsDf.index)
    pc1_prop_explained, pc2_prop_explained = My_pcoa.proportion_explained.loc["PC1"], My_pcoa.proportion_explained.loc["PC2"]
    return PcoaDf, pc1_prop_explained, pc2_prop_explained

def envfit_process_for_R(envFile):
    env= pd.read_csv(envFile, sep='\t').replace('Nan', '',regex=True).replace('.','')
    #### keep the entries that correspond to estimated values for all timepoints together (bodysite seperated)
    env = env[env['Timepoint'] == 'all']
    env["Variable"] = env['Variable'].replace({'AntiCumIndex':'Antibiotics','oxygenCumIndex':'Oxygenation', "il6WinIndex":'IL-6'})
    ColorD = {'dbCluster':'#f1c50d', 'Cluster':"#f1c50d",
        'BPD':'#d16666','ROP':'#d16666','IVH':'#d16666',
        'BW':'#3581bd','SEX':'#3581bd','GA':'#3581bd','DM':'#3581bd',
        'Age':'#e01f93',
        'Antibiotics':'#21cc8e','IL-6':'#21cc8e','Oxygenation':"#21cc8e"}
    env['Color'] = env['Variable'].map(ColorD)
    groupDict = {'Antibiotics':'MedIndex','Oxygenation':'MedIndex','IL-6':'MedIndex', 'SEX':'char', 'GA':'char','DM':'char', 'Age':'time'}
    env['grouping'] = env.Variable.map(groupDict)
    env.SampleType.replace({'Oral': 'Oral cavity'}, inplace=True)
    #env = env.melt(value_vars='r2', value_name='r2', id_vars=['Variable','SampleType'])
    env.to_csv('/media/christos/ssd/work/Infants/tmp/envFit_to_R.tsv',sep='\t', index=False)
    return env

def find_best_number_of_clusters(CountsDf, MetaDf, metric='braycurtis', exportFolder=""):
    plt.close()
    #import importlib
    #import matplotlib as mpl
    #importlib.reload(mpl); importlib.reload(plt); importlib.reload(sns)
    pcoaDf, v1, v2 = pcoa_from_dataframe(CountsDf, metric=metric)
    from sklearn.metrics import calinski_harabasz_score
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score, silhouette_samples
    import matplotlib.cm as cm
    ch_scores = []
    sil_scores = []
    for n_clusters in range(2,11):
        fig, (ax1, ax2) = plt.subplots(1,2)
        fig.set_size_inches(18,7)
        ax1.set_ylabel([0, len(pcoaDf.values) + (n_clusters+1)*10])
        kmeans = KMeans(n_clusters=n_clusters, random_state=9).fit(pcoaDf.values)
        kmean_clusters = kmeans.labels_
        ch_score = calinski_harabasz_score(pcoaDf.values, kmean_clusters)
        sample_silhoutte_values = silhouette_samples(pcoaDf.values, kmean_clusters)
        sil_score = silhouette_score(pcoaDf.values, kmean_clusters)
        sil_scores.append(sil_score)
        ch_scores.append(ch_score)
        y_lower = 10
        for i in range(n_clusters):
            ith_cluster_silhouette_values =  sample_silhoutte_values[kmean_clusters == i]
            ith_cluster_silhouette_values.sort()
            size_cluster_i = ith_cluster_silhouette_values.shape[0]
            y_upper = y_lower + size_cluster_i
            color = cm.nipy_spectral(float(i) / n_clusters)
            ax1.fill_betweenx(np.arange(y_lower, y_upper),
                          0, ith_cluster_silhouette_values,
                          facecolor=color, edgecolor=color, alpha=0.7)

            ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))
            y_lower = y_upper + 10
        # The vertical line for average silhouette score of all the values
        ax1.axvline(x=sil_score, color="red", linestyle="--")
        ax1.set_yticks([])  # Clear the yaxis labels / ticks
        ax1.set_title("The silhouette plot for the various clusters.")
        ax1.set_xlabel("The silhouette coefficient values")
        ax1.set_ylabel("Cluster label")
        colors = cm.nipy_spectral(kmean_clusters.astype(float) / n_clusters)
        ax2.scatter(pcoaDf.values[:, 0], pcoaDf.values[:, 1], marker='.', s=30, lw=0, alpha=0.7,
                c=colors, edgecolor='k')
        # Labeling the clusters
        centers = kmeans.cluster_centers_
        # Draw white circles at cluster centers
        ax2.scatter(centers[:, 0], centers[:, 1], marker='o',
        c="white", alpha=1, s=200, edgecolor='k')

        for i, c in enumerate(centers):
            ax2.scatter(c[0], c[1], marker='$%d$' % i, alpha=1,
                    s=50, edgecolor='k')
        ax2.set_title("The visualization of the clustered data.")
        ax2.set_xlabel("PC1")
        ax2.set_ylabel("PC2")
        plt.suptitle(("Silhouette analysis for KMeans clustering on sample data "
                  "with n_clusters = %d" % n_clusters),
                 fontsize=14, fontweight='bold')
        if exportFolder == "":
            fig.savefig('/media/christos/ssd/work/Infants/figures/Silhouette_figures'+str(n_clusters)+'.svg')
        else:
            fig.savefig(exportFolder+'Silhouette_figures' + str(n_clusters) + '.svg')
    #plt.show()
    plt.close()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot([i for i in range(2,11)], np.log(ch_scores)/10.0, label="Calinski-Harabz score (log /10)")
    ax.plot([i for i in range(2,11)], np.array(sil_scores), label="Sihlouette score")
    ax.legend()
    plt.xlabel("Number of Clusters")
    plt.tight_layout()
    if exportFolder == "":
        fig.savefig('/media/christos/ssd/work/Infants/figures/number-of-clusters.svg')
        fig.savefig('/media/christos/ssd/work/Infants/figures/number-of-clusters.png')
    else:
        fig.savefig(exportFolder+'number-of-clusters.svg')
        fig.savefig(exportFolder +'number-of-clusters.png')
    #plt.show()
def clustering_general(Df, MetaDf, metric='braycurtis', nclusters=4, eps=0.1, outfolder=''):
    MasterDf = Df.merge(MetaDf, on="SampleID", how='inner')
    sns.set(style="whitegrid")
    plt.rcParams.update({'axes.labelsize': 'x-large',
                         'xtick.labelsize': 'large',
                         'ytick.labelsize': 'large'})

    pcoaDf, v1, v2 = pcoa_from_dataframe(Df, metric=metric)
    import sklearn.cluster as skc
    from sklearn.cluster import DBSCAN
    from sklearn import metrics
    kmeans = skc.KMeans(n_clusters=nclusters).fit(pcoaDf.values)
    # import sklearn.preprocessing.StandardScaler as StandardScaler
    from sklearn import preprocessing
    # import sklearn.preprocessing.StandardScaler as StandardScaler

    X = preprocessing.StandardScaler().fit_transform(pcoaDf.values)
    db = DBSCAN(eps=eps).fit(X)
    dbscan_clusters = DBSCAN().fit_predict(X)
    # print (dbscan_clusters)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_
    labels_true = MasterDf["Patient"]
    ##number of clusters in labels ignoring noise if present (-1 is the 'cluster' label for noisy samples (not core)
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    clusters = kmeans.labels_
    db_clusters = labels
    MasterDf = MasterDf.assign(PC1=pcoaDf.PC1)
    MasterDf = MasterDf.assign(PC2=pcoaDf.PC2)
    MasterDf = MasterDf.assign(Cluster=clusters)  ### these are the clusters from kmeans
    MasterDf = MasterDf.assign(dbCluster=db_clusters)
    MostAbundant = pd.DataFrame({'Cluster': [], 'label': []})
    for name, cdf in MasterDf.groupby("Cluster"):
        numDf = pd.DataFrame({"Median": cdf[Df.columns].apply(np.median, axis=0),
                              "Mean": cdf[Df.columns].apply(np.mean, axis=0),
                              'MeanRelative': relativeAbundance_from_dataframe(cdf[Df.columns]).apply(np.mean,
                                                                                                            axis=0),
                              "Cluster": [name for x in Df.columns],
                              "taxon": Df.columns})
        numDf = numDf.sort_values(["MeanRelative", "Median", "Mean"], ascending=False).head(1)

        numDf['label'] = numDf['taxon']  # + ["_" for x in numDf.index] + numDf['Cluster'].map(str)
        MostAbundant = MostAbundant.append(pd.DataFrame({'Cluster': numDf['Cluster'], 'label': numDf['label']}))
        # keySpeciesDict[name] = numDf.iloc[0]["Taxa"]+'_'+str(name)
    MostAbundant.set_index('Cluster', inplace=True)

    MostAbundantDBscan = pd.DataFrame({'dbCluster': [], 'label': []})
    for name, cdf in MasterDf.groupby("dbCluster"):
        numDf = pd.DataFrame({"Median": cdf[Df.columns].apply(np.median, axis=0),
                              "Mean": cdf[Df.columns].apply(np.mean, axis=0),
                              'MeanRelative': relativeAbundance_from_dataframe(cdf[Df.columns]).apply(np.mean,
                                                                                                            axis=0),
                              "dbCluster": [name for x in Df.columns],
                              "taxon": Df.columns})
        # numDf = numDf.assign(Percent=numDf["Counts"]*100/numDf["Counts"].sum())
        numDf = numDf.sort_values(['MeanRelative', "Median", "Mean"], ascending=False).head(1)
        numDf['label'] = numDf['taxon']  # + ["_" for x in numDf.index] + numDf['dbCluster'].map(str)
        MostAbundantDBscan = MostAbundantDBscan.append(
            pd.DataFrame({'dbCluster': numDf['dbCluster'], 'label': numDf['label']}))
    MostAbundantDBscan.set_index('dbCluster', inplace=True)

    centroids = kmeans.cluster_centers_
    pcoaDf = pcoaDf.assign(Cluster=clusters)  ### these are the clusters from kmeans

    ### DBSCAN plotting
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)
    unique_labels = set(labels)  ### dbscan labels
    colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]
    for k, col in zip(unique_labels, colors):
        if k == -1:
            ## Black used for noise
            col = [0, 0, 0, 1]
        class_member_mask = (labels == k)
        xy = X[class_member_mask & core_samples_mask]
        db_centroid_x, db_centroid_y = np.mean(xy[:, 0]), np.mean(xy[:, 1])

        ax.scatter(xy[:, 0], xy[:, 1], marker='o', color=tuple(col), edgecolor='k', label=k, s=90)
        xy = X[class_member_mask & ~core_samples_mask]

        ax.scatter(xy[:, 0], xy[:, 1], marker='o', color=tuple(col), edgecolor='k', s=30)
        if k != -1:
            ax.annotate(MostAbundantDBscan.loc[k, 'label'], (db_centroid_x+0.05, db_centroid_y+0.05))
    handles, labels = ax.get_legend_handles_labels()
    legend1 = ax.legend(handles, labels, title=r'dbscan Cluster', fontsize=14, frameon=True, bbox_to_anchor=(1.301, 1))
    plt.title('PCoA-BrayCurtis: Estimated number of clusters: %d' % n_clusters_)
    fig.text(0.03, 0.5, "PC2 [{0:.0%}]".format(v2), rotation='vertical', va='center', fontsize=16)
    fig.text(0.4, 0.03, "PC1 [{0:.0%} ]".format(v1), va='center', fontsize=16)
    # fig.tight_layout()
    plt.subplots_adjust(right=0.8)
    plt.show()
    #fig.savefig(outfolder + "/Dbscan_clusters.png")
    #fig.savefig(outfolder + "/Dbscan_clusters.svg")

    # plt.close()

    fig = plt.figure(figsize=(10, 8))
    sizefactor = 2
    # pal = dict(Stool=sns.xkcd_rgb["brown"], Skin=sns.xkcd_rgb["peach"], Oral=sns.xkcd_rgb["salmon"])
    ax1 = fig.add_subplot(111)
    grouped = pcoaDf.groupby("Cluster")
    pal = sns.color_palette("Set2", grouped.ngroups)
    Index = 0
    for name, group in grouped:
        for x, y in zip(group["PC1"], group["PC2"]):
            ax1.plot([centroids[int(name)][0], x], [centroids[int(name)][1], y], color='k', alpha=0.2)
        ax1.scatter(group["PC1"], group["PC2"], color=pal[Index], s=90, alpha=0.8, marker='o', edgecolor='k',
                    label=name)
        ax1.scatter(centroids[int(name)][0], centroids[int(name)][1], s=100, color="white", marker="^", edgecolor='k',
                    linewidths=1)
        # ax1.annotate(keySpeciesDict[name], (centroids[int(name)][0],centroids[int(name)][1]))
        ax1.annotate(MostAbundant.loc[name, 'label'], (centroids[int(name)][0]+0.01, centroids[int(name)][1]+0.01))
        Index += 1
    handles, labels = ax1.get_legend_handles_labels()
    legend1 = ax1.legend(handles, labels, title=r'kmeans Cluster', fontsize=14, frameon=True, bbox_to_anchor=(1.305, 1))
    legend1.get_title().set_fontsize(14)
    fig.text(0.03, 0.5, "PC2 [{0:.0%}]".format(v2), rotation='vertical', va='center', fontsize=16)
    fig.text(0.4, 0.03, "PC1 [{0:.0%}]".format(v1), va='center', fontsize=16)
    # fig.suptitle(title)
    # fig.tight_layout()
    plt.subplots_adjust(right=0.8)
    plt.show()
    #fig.savefig(outfolder + "/Kmeans_clusters.png")
    #fig.savefig(outfolder + '/kmeans_clusters.svg')

    # plt.close()
    MasterDf["ClusterASV"] = MasterDf["Cluster"].apply(lambda x: MostAbundant.loc[x, 'label'])
    MasterDf['dbClusterASV'] = MasterDf['dbCluster'].apply(
        lambda x: MostAbundantDBscan.loc[x, 'label'] if x != -1 else 'Unclustered')
    print('PC1 axis variability: {0:.0%}'.format(v1))
    print('PC2 axis variability: {0:.0%}'.format(v2))
    return MasterDf



def clustering(CountsDf, MetaDf, TaxDf, metric='braycurtis', nclusters=2, taxlevel=None, eps=0.1, outfolder=''):
    ### Returns a pandas Series with the cluster information, in the same order as in MetaDf (Cluster col for kmeans and dbCluster for dbscan
    #plt.close("all")

    if taxlevel is not None:
        CountsDf = convert_df_to_taxlevel_from_myClassifier(CountsDf, TaxDf, taxlevel=taxlevel)
    MasterDf = CountsDf.merge(MetaDf, on="SampleID", how='inner')
    sns.set(style="whitegrid")
    plt.rcParams.update({'axes.labelsize': 'x-large',
                         'xtick.labelsize': 'large',
                         'ytick.labelsize': 'large'})


    pcoaDf, v1, v2 = pcoa_from_dataframe(CountsDf, metric=metric)
    import sklearn.cluster as skc
    from sklearn.cluster import DBSCAN
    from sklearn import metrics
    kmeans = skc.KMeans(n_clusters=nclusters).fit(pcoaDf.values)
    #import sklearn.preprocessing.StandardScaler as StandardScaler
    from sklearn import preprocessing
    #import sklearn.preprocessing.StandardScaler as StandardScaler

    X = preprocessing.StandardScaler().fit_transform(pcoaDf.values)
    db = DBSCAN(eps=eps).fit(X)
    dbscan_clusters = DBSCAN().fit_predict(X)
    #print (dbscan_clusters)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_
    labels_true = MasterDf["Patient"]
    ##number of clusters in labels ignoring noise if present (-1 is the 'cluster' label for noisy samples (not core)
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    #print('Estimated number of clusters: %d' % n_clusters_)
    #print("Homogeneity: %0.3f" % metrics.homogeneity_score(labels_true, labels))
    #print("Completeness: %0.3f" % metrics.completeness_score(labels_true, labels))
    #print("V-measure: %0.3f" % metrics.v_measure_score(labels_true, labels))
    #print("Adjusted Rand Index: %0.3f"% metrics.adjusted_rand_score(labels_true, labels))
    #print("Adjusted Mutual Information: %0.3f"% metrics.adjusted_mutual_info_score(labels_true, labels))
    #print("Silhouette Coefficient: %0.3f"% metrics.silhouette_score(pcoaDf.values, labels))
    clusters = kmeans.labels_
    #db_clusters = labels
    MasterDf = MasterDf.assign(PC1= pcoaDf.PC1)
    MasterDf = MasterDf.assign(PC2 = pcoaDf.PC2)
    MasterDf = MasterDf.assign(Cluster=clusters) ### these are the clusters from kmeans
    #MasterDf = MasterDf.assign(dbCluster=db_clusters)
    keySpeciesDict ={}
    #MostAbundant=MasterDf.groupby('Cluster')[CountsDf.columns].mean().idxmax(axis=1).to_frame().rename(columns={0:'taxon'})
    #MostAbundant = MasterDf.groupby('Cluster')[CountsDf.columns].agg(['mean','median']).sort_values(['median','mean'], ascending=False).head(1)#idxmax(axis=1).to_frame().rename(columns={0: 'taxon'})
    #print(MostAbundant)
    #MostAbundant['taxon'] = MostAbundant['taxon'].map(dict(zip(TaxDf.index, TaxDf['genus'])))
    #MostAbundant['label'] = MostAbundant['taxon'] +["_" for x in MostAbundant.index]+MostAbundant.index.map(str)
    MostAbundant = pd.DataFrame({'Cluster': [], 'label': []})
    for name, cdf in MasterDf.groupby("Cluster"):
        numDf = pd.DataFrame({"Median": cdf[CountsDf.columns].apply(np.median, axis=0),
                              "Mean": cdf[CountsDf.columns].apply(np.mean, axis=0),
                              'MeanRelative': relativeAbundance_from_dataframe(cdf[CountsDf.columns]).apply(np.mean,axis=0),
                              "Cluster": [name for x in CountsDf.columns],
                              "taxon": CountsDf.columns})
        # numDf = numDf.assign(Percent=numDf["Counts"]*100/numDf["Counts"].sum())
        #print(count_ASVs_per_genus(cdf[CountsDf.columns], convert_df_to_taxlevel_from_myClassifier(cdf[CountsDf.columns],TaxDf,taxlevel='genus'), TaxDf))
        numDf = numDf.sort_values(["MeanRelative", "Median", "Mean"], ascending=False).head(1)

        numDf['label'] = numDf['taxon']# + ["_" for x in numDf.index] + numDf['Cluster'].map(str)
        MostAbundant = MostAbundant.append(pd.DataFrame({'Cluster': numDf['Cluster'], 'label': numDf['label']}))
        # keySpeciesDict[name] = numDf.iloc[0]["Taxa"]+'_'+str(name)
    MostAbundant.set_index('Cluster',inplace=True)
    centroids = kmeans.cluster_centers_
    pcoaDf = pcoaDf.assign(Cluster=clusters)  ### these are the clusters from kmeans

    #MostAbundantDBscan = MasterDf.groupby('dbCluster')[CountsDf.columns].mean().idxmax(axis=1).to_frame().rename(columns={0:'taxon'})
    #MostAbundantDBscan['taxon'] = MostAbundantDBscan['taxon'].map(dict(zip(TaxDf.index, TaxDf['genus'])))
    #MostAbundantDBscan['label'] = MostAbundantDBscan['taxon']+["_" for x in MostAbundantDBscan.index] + MostAbundantDBscan.index.map(str)
    #print(MostAbundant.columns)
    # MostAbundantDBscan = pd.DataFrame({'dbCluster': [], 'label': []})
    # for name, cdf in MasterDf.groupby("dbCluster"):
    #     numDf = pd.DataFrame({"Median": cdf[CountsDf.columns].apply(np.median, axis=0),
    #                           "Mean": cdf[CountsDf.columns].apply(np.mean, axis=0),
    #                           'MeanRelative': relativeAbundance_from_dataframe(cdf[CountsDf.columns]).apply(np.mean,axis=0),
    #                           "dbCluster": [name for x in CountsDf.columns],
    #                           "taxon": CountsDf.columns})
    #     # numDf = numDf.assign(Percent=numDf["Counts"]*100/numDf["Counts"].sum())
    #     numDf = numDf.sort_values(['MeanRelative', "Median", "Mean"], ascending=False).head(1)
    #     numDf['label'] = numDf['taxon'] #+ ["_" for x in numDf.index] + numDf['dbCluster'].map(str)
    #     MostAbundantDBscan = MostAbundantDBscan.append(pd.DataFrame({'dbCluster': numDf['dbCluster'], 'label': numDf['label']}))
    # MostAbundantDBscan.set_index('dbCluster',inplace=True)
    #

    #
    # ### DBSCAN plotting
    # fig = plt.figure(figsize=(10,8))
    # ax = fig.add_subplot(111)
    # unique_labels = set(labels) ### dbscan labels
    # colors = [plt.cm.Spectral(each) for each in np.linspace(0,1, len(unique_labels))]
    # for k, col in zip(unique_labels, colors):
    #     if k == -1:
    #         ## Black used for noise
    #         col = [0,0,0,1]
    #     class_member_mask = (labels == k)
    #     xy = X[class_member_mask & core_samples_mask]
    #     db_centroid_x, db_centroid_y =  np.mean(xy[:,0]), np.mean(xy[:,1])
    #
    #     ax.scatter(xy[:,0], xy[:,1],marker='o', color=tuple(col),edgecolor='k', label=k, s=90)
    #     xy = X[class_member_mask & ~core_samples_mask]
    #
    #     ax.scatter(xy[:,0], xy[:,1],marker= 'o', color=tuple(col), edgecolor='k',s=30)
    #     if k != -1:
    #         ax.annotate(MostAbundantDBscan.loc[k, 'label'], (db_centroid_x, db_centroid_y))
    # handles, labels = ax.get_legend_handles_labels()
    # legend1 = ax.legend(handles, labels, title=r'dbscan Cluster', fontsize=14,frameon=True, bbox_to_anchor=(1.301,1))
    # plt.title('PCoA-BrayCurtis: Estimated number of clusters: %d' % n_clusters_)
    # fig.text(0.03, 0.5, "PC2 [{0:.0%}]".format(v2) ,rotation='vertical',va='center', fontsize=16)
    # fig.text(0.4,0.03, "PC1 [{0:.0%} ]".format(v1), va='center',fontsize=16)
    # #fig.tight_layout()
    # plt.subplots_adjust(right=0.8)
    # #fig.savefig(outfolder+"Dbscan_clusters.png")
    # #fig.savefig(outfolder+"Dbscan_clusters.svg")
    #
    # #plt.close()

    fig = plt.figure(figsize=(10,8))
    sizefactor = 2
    #pal = dict(Stool=sns.xkcd_rgb["brown"], Skin=sns.xkcd_rgb["peach"], Oral=sns.xkcd_rgb["salmon"])
    ax1 = fig.add_subplot(111)
    grouped = pcoaDf.groupby("Cluster")
    pal = sns.color_palette("Set2", grouped.ngroups)
    Index = 0
    for name, group in grouped:
        for x, y in zip(group["PC1"], group["PC2"]):
            ax1.plot([centroids[int(name)][0], x], [centroids[int(name)][1], y], color='k', alpha=0.2)
        ax1.scatter(group["PC1"],group["PC2"], color=pal[Index],s=90,alpha=0.8, marker='o', edgecolor='k', label=name)
        ax1.scatter(centroids[int(name)][0],centroids[int(name)][1],s=100, color="white", marker="^",edgecolor='k', linewidths=1)
        #ax1.annotate(keySpeciesDict[name], (centroids[int(name)][0],centroids[int(name)][1]))
        ax1.annotate(MostAbundant.loc[name,'label'], (centroids[int(name)][0],centroids[int(name)][1]) )
        Index += 1
    handles, labels = ax1.get_legend_handles_labels()
    legend1 = ax1.legend(handles, labels, title=r'kmeans Cluster', fontsize=14,frameon=True, bbox_to_anchor=(1.305,1))
    legend1.get_title().set_fontsize(14)
    fig.text(0.03, 0.5, "PC2 [{0:.0%}]".format(v2) ,rotation='vertical',va='center', fontsize=16)
    fig.text(0.4,0.03, "PC1 [{0:.0%}]".format(v1) , va='center',fontsize=16)
    #fig.suptitle(title)
    #fig.tight_layout()
    plt.subplots_adjust(right=0.8)
    fig.savefig(outfolder+"Kmeans_clusters.png")
    fig.savefig(outfolder+'kmeans_clusters.svg')

    #plt.close()
    MasterDf["ClusterASV"] = MasterDf["Cluster"].apply(lambda x:MostAbundant.loc[x, 'label'])
    #MasterDf['dbClusterASV'] = MasterDf['dbCluster'].apply(lambda x:MostAbundantDBscan.loc[x, 'label'] if x!= -1 else 'Unclustered')
    print('PC1 axis variability: {0:.0%}'.format(v1))
    print('PC2 axis variability: {0:.0%}'.format(v2))
    return MasterDf


def plot_asv_timeseries_for_patients(CountsDf, MetaDf, TaxDf, valueCol="",typeCol="SampleType", taxlevel='asv',taxonomic_label=None,groupCol="Patient", timeCol="Age"):
    #maximum_value = MetaDf[valueCol].max()
    sns.set()
    Df = convert_df_to_taxlevel_from_myClassifier(CountsDf, TaxDf, taxlevel=taxlevel)
    Df = relativeAbundance_from_dataframe(Df)
    MasterDf = Df.merge(MetaDf, on="SampleID", how='inner')
    MasterDf.sort_values(by=[groupCol, timeCol], inplace=True)
    pal = {'Gut':sns.xkcd_rgb["brown"], 'Skin':sns.xkcd_rgb["peach"], 'Oral_cavity':sns.xkcd_rgb["salmon"]}

    g = sns.FacetGrid(MasterDf, col=groupCol, hue=typeCol, col_wrap=4, size=2, palette=pal, sharex=False, sharey=True)
    g = (g.map(plt.plot, timeCol, valueCol, marker=".", ms=15).set(xticks=[1,3,6,9,12,16]).set_titles("{col_name}").add_legend())
    if taxonomic_label is not None:
        #TaxDf.set_index('otu', inplace=True)
        g.set_axis_labels("Age", "Relative abundance (%)\n"+TaxDf.set_index('asv').loc[valueCol, taxonomic_label])
    #plt.setp(g.get_legend().get_texts(), fontsize='22') # for legend text
    #plt.setp(g.get_legend().get_title(), fontsize='32') # for legend title
    #g.set_axis_labels("Age", "Richness")
    #g.text(0.02, 0.5, "Age", va='center', rotation='vertical', fontsize=14, ha='center')
    #g.text(0.03, 0.5, "Richness", va='center', rotation='vertical', fontsize=10, ha='center')
    g.savefig('/media/christos/ssd/work/Infants/figures/8.'+valueCol+"_trajectories.png", bbox_inches='tight')
def random_forest_analysis(MasterDf, MetaDf, TaxDf, taxlabel='genus', predictCol='Cluster', outfilename='10.png'):
    ## MasterDf is merge CountsDf, metaDf and has extra info from the clustering.
    """Random forest analysis"""
    from sklearn.model_selection import train_test_split
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.preprocessing import StandardScaler
    """We need to separate data into the features: data we use for the prediction and target: the data we want to predict."""
    taxCols = [x for x in MasterDf.columns if x not in MetaDf.columns and x != 'Cluster']
    features = MasterDf[taxCols]
    target = np.array(MasterDf[predictCol])
    #train_features, test_features, train_labels, test_labels = train_test_split(features, target, test_size=0.2, random_state=95)
    features = StandardScaler().fit_transform(features)

    #test_features = StandardScaler(test_features)
    #target_onehotEncode = pd.get_dummies(masterdf["Cluster"])
    #target = np.array(target_onehotEncode)
    clf = RandomForestClassifier(n_estimators=1000, random_state=95)
    """ Train the model. Remove one OTU at a time from the training dataset and evaluate the accuracy of the trained classifier"""
    accuracies = []
    #for feature in taxCols:
    #    newfeatures = features.drop(columns=feature).values.astype(float)
        #newfeatures = newfeatures.values.astype(float)
    #    clf.fit(newfeatures,target)
    #    accuracies.append(clf.score(newfeatures, target))
    ##plot accuracy values
    #ind = np.arange(len(accuracies))
    #plt.bar(ind, accuracies, orientation='vertical')
    #plt.xticks(ind, taxCols, rotation='vertical')
    #plt.ylabel("Mean Accuracy");plt.xlabel("ASV")
    #plt.title("Mean Accuracy removing ASV")
    #plt.close()
    """ Now train classifier using the entire training dataset."""
    clf.fit(features, target)
    ## Make predictions
    predictions = clf.predict(features)
    ## i can also predict the probability for each class (cluster number)
    ## clf.predict_proba(features)
    errors = sum([1.0 for x,y in zip(target, predictions) if x!=y])
    accuracy = 100 - (100* (errors/len(predictions)))
    print ("This model's accuracy {}%".format(accuracy))
    importances = list(clf.feature_importances_)
    print ("This is the score (mean accurary): {}".format(clf.score(features, target)))
    print (importances)
    feature_importances = sorted([(feature,round(importance,2)) for feature,importance in zip(taxCols, importances) if importance >= 0.01] , key= lambda x:x[1],reverse=True)
    [print("ASV: {:20} Importance: {}".format(*pair)) for pair in feature_importances if pair[1] > 0]
    ##plot importance values
    plt.close()
    fig, ax = plt.subplots()
    ind = np.arange(len(feature_importances))
    plt.bar(ind, [x[1] for x in feature_importances], orientation='vertical')
    plt.xticks(ind,[TaxDf.loc[x[0],taxlabel] for x in feature_importances], rotation='vertical')
    plt.ylabel("Importance");plt.xlabel("Variable")
    plt.title("Variable Importance")
    fig.savefig(outfilename,  bbox_inches='tight')
    
def plot_cluster_info_for_patients(MasterDf, CountsDf, outfilename):
    """MasterDf should contain a column Cluster informed from the clustering performed previously."""
    Df = MasterDf[MasterDf.columns.difference(CountsDf.columns)]
    Df.reset_index(inplace=True)
    cluster_abundance_dict = dict(zip(Df['dbClusterASV'].unique(), [x*100 for x in range(len(Df['dbCluster'].unique()))]))
    Df['Cluster_abundance'] = Df['dbClusterASV'].map(cluster_abundance_dict)
    Df.to_csv('/media/christos/ssd/work/Infants/tmp/ClusterInfo.tsv', sep='\t', index=False)
    subprocess.run('Rscript /media/christos/ssd/work/Infants/scripts/plot_cluster_per_patient.R '
                   '/media/christos/ssd/work/Infants/tmp/ClusterInfo.tsv '+outfilename, shell=True)
    #os.remove("/media/christos/ssd/work/Infants/tmp/toR.tsv")

def plot_richness_age(MetaDf,valueCol="Richness",bodysite="", imageName=''):
    #sns.set()
    pal = dict(Gut=sns.xkcd_rgb["brown"], Skin=sns.xkcd_rgb["peach"], Oral=sns.xkcd_rgb["salmon"])
    g = (sns.jointplot("Age", valueCol, xlim=(0,MetaDf["Age"].max()),data=MetaDf[MetaDf["SampleType"]==bodysite], kind="reg", space=0, color=pal[bodysite], height=7, label='big').set_axis_labels("Age",valueCol, fontsize=18))

    #g.set_xlabel("Age", fontsize=18)
    #g.set_ylabel("Richness", fontsize=18)
    g.savefig("/media/christos/ssd/work/Infants/figures/8."+bodysite+".png", figsize=(8,8))
def plot_diversity_boxplots_for_groups(MetaDf, groupcolumn, valuecolumn, imageName="boxplotsgroups.svg",show=False):
    plt.close()
    plt.style.use('seaborn-whitegrid')
    from matplotlib import colors as mcolors
    colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
    colors = [c for c in colors.keys()]
    grouped = MetaDf.groupby(by=groupcolumn)
    np.random.seed(96)
    colors = [colors[np.random.randint(0, len(colors))] for x in range(grouped.ngroups)]
    fig, ax = plt.subplots()
    medianprops = dict(linestyle='-', linewidth=2.5, color='firebrick')
    meanpointprops = dict(marker='D', markeredgecolor='black',markerfacecolor='firebrick')
    for name, group in grouped:
        bplot = ax.boxplot([group[valuecolumn] for name,group in grouped], labels =[name for name, group in grouped],medianprops=medianprops,meanprops=meanpointprops, patch_artist=True)

    for bp in (bplot):
        for patch, color in zip(bplot['boxes'], colors):
            patch.set_facecolor(color)
    fig.text(0.03, 0.5, valuecolumn, va='center', rotation='vertical', fontsize=12)
    fig.text(0.45, 0.03, groupcolumn, va='center', fontsize=12)
    plt.show()
    #fig.savefig(imageName,bbox_inches='tight')
    plt.close()
def plot_alpha_diversity_age_scatter(MetaDf, alphaName='Richness'):
    sns.set(style="darkgrid")
    pal = {'Gut':sns.xkcd_rgb["brown"], 'Skin':sns.xkcd_rgb["peach"], 'Oral_cavity':sns.xkcd_rgb["salmon"]}

    for name, group in MetaDf.groupby('SampleType'):
        g = sns.jointplot("Age", alphaName, data=group, kind="reg",xlim=(0, 16), height=7, color=pal[name])

        #rsquare = scipy.stats.pearson(group['Age'], group[alphaName])[0]
        g= g.annotate(scipy.stats.pearsonr)
        #g = sns.JointGrid(x='Age', y=alphaName, data=group)
        #g = g.plot_joint(sns.regplot, color=pal[name])
        #g = g.plot_marginals(sns.distplot, color=pal[name])
        #g = g.annotate(scipy.stats.pearsonr)
        g.savefig("/media/christos/ssd/work/Infants/figures/7."+name+"_"+alphaName+"-age.svg")

def distance_between_sampletypes_overtime(CountsDf, MetaDf, metric='braycurtis'):
    import scipy
    from skbio.stats.distance import DistanceMatrix
    import itertools
    MetaDf = MetaDf.loc[[i for i in MetaDf.index if i in CountsDf.index],:]
    Distance_matrix = DistanceMatrix(scipy.spatial.distance.pdist(CountsDf.values, metric=metric)).data
    Df = pd.DataFrame(Distance_matrix, index=CountsDf.index, columns=CountsDf.index)
    pairnames = [a + "-" + b for a, b in itertools.combinations(MetaDf["SampleType"].sort_values().unique(), 2)]
    NewDf = pd.DataFrame()
    for name, group in MetaDf.groupby(["Age","Patient"]):
        seen_pairs = []
        for sample1, sample2 in itertools.combinations(group.index, 2):
            #get_samples = [(x, y) for x, y in zip(group.index, group["SampleType"])]
            if sample1+sample2 not in seen_pairs: ### to avoid the reversed pairs..
                seen_pairs.append(sample2+sample1)
                distance = Df.loc[sample1, sample2]
                if MetaDf.loc[sample1, "SampleType"] == MetaDf.loc[sample2, "SampleType"]:
                    comparison = ' Distance within body site'
                else:
                    comparison = ' Distance between body sites'
                NewDf = NewDf.append(pd.DataFrame({"Comparison":[comparison],"Age":[name[0]], "Patient":[name[1]],"Distance":[distance],'Pair':[ MetaDf.loc[sample1, "SampleType"]+"-"+ MetaDf.loc[sample2, "SampleType"]]}),sort=True)
    plt.close('all')
    NewDf.sort_values(by=['Patient',"Age"], inplace=True)

    g = sns.FacetGrid(NewDf, col="Patient", hue='Pair',col_wrap=4, height=2, sharex=False,sharey=True)
    g = (g.map(plt.plot, "Age", "Distance",marker=".", ms=15).set(xticks=[1,3,6,9,12,16]).set_titles("{col_name}").add_legend())
    #plt.suptitle(title)
    g.savefig("/media/christos/ssd/work/Infants/figures/9.Distance_between_body_sites_overtime_per_patient.png")

def run_permanova(MetaDf, CountsDf, distance_metric='', factor=''):
    """Permutational Multivariate Analysis of Variance (PERMANOVA) is a non-parametric method that tests whether two or more groups of objects (e.g., samples) are significantly different
    based on a categorical factor. It is conceptually similar to ANOVA except that it operates on a distance matrix, which allows for multivariate analysis. PERMANOVA computes a pseudo-F statistic.
    Statistical significance is assessed via a permutation test. The assignment of objects to groups (grouping) is randomly permuted a number of times (controlled via permutations).
    A pseudo-F statistic is computed for each permutation and the p-value is the proportion of permuted pseudo-F statisics that are equal to or greater than the original (unpermuted) pseudo-F statistic.
    """
    from skbio.stats.distance import permanova
    import scipy
    from skbio.stats.distance import DistanceMatrix
    Df = CountsDf.copy()
    MasterDf = MetaDf.merge(Df, on='SampleID', how='inner')
 #### do i need to use some transformation ?? root, clr?

    Distance_matrix = DistanceMatrix(scipy.spatial.distance.pdist(MasterDf[Df.columns].values, metric=distance_metric), ids=MasterDf.index)
    #print("Running PERMANOVA...")
    #factors =["infantID","Age","SampleType" ]
    ## create new combinational factors
    #PatientplusAge = MasterDf['DM']+MasterDf['Age'].astype(str)
    results = pd.DataFrame()
    for factor in MasterDf[MetaDf.columns]:
        #results = permanova(Distance_matrix, MasterDf, column=factor, permutations=1000)
        try:
            perm = permanova(Distance_matrix, MasterDf, column=factor, permutations=1000)
            thisDf  = pd.DataFrame(dict(zip(["F-statistic", "p-value", 'factor'],[[perm['test statistic']], [perm['p-value']], [factor]])))
            results = results.append(thisDf)
        except:
            print("Exception: All values in the grouping vector are the same/ or all different. Column:'"+factor+"' could not be used as a factor!! ")
        #print ("PERMANOVA test statistic: "+str(results['test statistic'])+" for "+factor+",  p-value: "+str(results['p-value']))
    #print(results.sort_values('F-statistic',ascending=False).to_string())
    return results.sort_values(['p-value','F-statistic'])

def run_cca(MetaDf, CountsDf, TaxDf):
    """Canonical correspondence analysis is a multivariate ordination technique.
    It appeared in community ecology [R88] and relates community composition to the variation in the environment (or in other factors).
    It works from data on abundances or counts of samples and constraints variables, and outputs ordination axes that maximize sample separation among species."""
    from skbio.stats.ordination import cca
    print("Running CCA(canonical correspondece analysis)...")
    Df = CountsDf.copy()
    MasterDf = Df.merge(MetaDf, on='SampleID', how='inner')
    constrains = MasterDf[["GA"]]
    CCA = cca(MasterDf[Df.columns], constrains)
    cca_samples = pd.DataFrame(CCA.samples[["CCA1","CCA2"]],index=MasterDf.index)
    cca_otus = pd.DataFrame(CCA.features[["CCA1","CCA2"]], index=Df.columns)
    cca_constraints_axes = CCA.biplot_scores[["CCA1","CCA2"]]
    print(cca_constraints_axes)

def differential_abundance_ancom(CountsDf, MetaDf, TaxDf, Taxlevel='asv', grouping=""):
    """Compares only two groups"""
    Df =CountsDf.copy().apply(lambda x: x+1) ### it only accepts positive values, so i create a  pseudocount
    ThisTax = TaxDf.loc[Df.columns]  ## select only the otus in the givent count table
    Df = convert_df_to_taxlevel_from_myClassifier(Df,TaxDf,taxlevel=Taxlevel)
    #NewTaxIndex  = TaxDf[Taxlevel]
    #ThisTax.reindex(Df.columns)
    MasterDf = Df.merge(MetaDf, on='SampleID', how='inner')
    #MetaDf = MetaDf.assign(Grouping=[np.nan for i in MetaDf.index])
    #Grouping_dict = {x:i for i,x in enumerate((MetaDf['BPD'].astype(str)+MetaDf['Patient']).unique())}
    #MetaDf['Grouping'] = (MetaDf["TWIN"].astype(str)+MetaDf['Patient']).map(Grouping_dict)
    Grouping = MasterDf[grouping]

    from skbio.stats.composition import ancom
    AncomDf, PercentileDf = ancom(MasterDf[Df.columns], Grouping)
    print(AncomDf.loc[AncomDf['Reject null hypothesis'] == True].sort_values('W', ascending=False)['W'].to_string())
    #print(ThisTax.loc[AncomDf.loc[AncomDf['Reject null hypothesis']== True].index])
    print(PercentileDf[50.0].loc[AncomDf.loc[AncomDf['Reject null hypothesis']== True].index].to_string())
    #print(PercentileDf.loc['5497318e515a8c328a68f95975d9c7d4'].to_string())
def distance_on_correlation_matrix(CountsDf, MetaDf, TaxDf, Taxlevel='asv'):
    """For the hierarchical clustering of the correlation matrix, in qiime (gneiss tutorial) they use ward's linkage method"""


    from scipy.spatial.distance import pdist
    import scipy.cluster.hierarchy as sch
    Df = CountsDf.copy().apply(lambda x:x+1)
    Df = convert_df_to_taxlevel_from_myClassifier(Df, TaxDf, taxlevel=Taxlevel)

    MasterDf = Df.merge(MetaDf, on='SampleID', how='inner')
    corr = MasterDf[Df.columns].corr()
    D = pdist(corr)
    D[D==0] = 10**-20
    L = sch.linkage(D, method='ward')

    print(np.where(L == 0)[0])

    #print(ilr(L))

    fig, ax = plt.subplots(figsize=(14,16))
    dn = sch.dendrogram(L, orientation='right', above_threshold_color='#bcbddc', labels=Df.columns)
    locs, labels = plt.yticks()
    plt.tight_layout()
    fig.savefig('/media/christos/ssd/work/Infants/figures/11.Species_corellation_dendogram.png')
    C = sch.fcluster(L, 0.5 , 'distance')
    Columns = [Df.columns.tolist()[i] for i in list((np.argsort(C)))]
    for column in corr.columns:
        if abs(corr[column].min()) > 0.5:
            print(column, corr[column].idxmin(), corr[column].min())
    Df = Df.reindex(Columns, axis=1)
    import matplotlib as mpl
    #mpl.style.use('ggplot')
    fig, ax = plt.subplots(figsize=(12,16))
    corr = Df.corr()
    cax = ax.matshow(corr, cmap='RdYlGn')
    plt.xticks(range(len(corr.columns)), corr.columns, rotation=90, size=6)
    plt.yticks(range(len(corr.columns)), corr.columns, size=6)
    cbar = fig.colorbar(cax, aspect=40, shrink=.6, spacing='uniform')
    plt.tight_layout()
    fig.savefig('/media/christos/ssd/work/Infants/figures/11.Species_correlation_heatmap.png')

    sns.set(font_scale=0.6)
    pal = {'Gut': sns.xkcd_rgb["brown"], 'Skin': sns.xkcd_rgb["peach"], "Oral_cavity": sns.xkcd_rgb["salmon"]}
    g = sns.clustermap(corr, method='ward', cmap='vlag_r')
    #row_colors = MasterDf.iloc[g.dendrogram_row.reordered_ind]['SampleType'].map(pal)
    #print(MasterDf.iloc[g.dendrogram_row.reordered_ind]['SampleType'])
    #print(row_colors)
    g.savefig('/media/christos/ssd/work/Infants/figures/11.Species_correlation_heatmap2.png', figsize=(14,16))

def run_gneiss(CountsDf, MetaDf, TaxDf):
    Df = CountsDf.T
    Df.to_csv("/media/christos/ssd/work/Infants/qiime2/gneiss-feature-table.tsv", sep='\t', index_label='#OTUID')

    MetaDf.loc[Df.columns].to_csv('/media/christos/ssd/work/Infants/qiime2/gneiss-meta-table.tsv', sep='\t', index_label='#SampleID')
    TaxDf.loc[Df.index].to_csv('/media/christos/ssd/work/Infants/qiime2/gneiss-taxa-table.tsv', sep='\t', index_label='#OTUID')
    subprocess.run('biom convert -i /media/christos/ssd/work/Infants/qiime2/gneiss-feature-table.tsv -o /media/christos/ssd/work/Infants/qiime2/gneiss-feature-table.biom --table-type="OTU table" --to-json',
        shell=True)


def check_for_cross_contamination(cross_tab, tax, out_dir=''):
    """CrossTabulated.txt file from Craig's cross-contamination analysis, taxonomy dataframe for my ASV dataset
    Returns a list of genera to be removed from current ASV dataset and exports a file in the taxonomy dir.
    The list includes genera that have been detected in higher abundance in samples from other datasets in the same sequencing run"""
    cr = pd.read_csv(cross_tab ,sep='\t', index_col='OTUnumber')
    ms = cr.filter(regex="SAMPLE")
    cr['MyMaxCount'] = ms.max(axis=1).values
    cr = cr[["crossedDataset", "highestSample", "highestSamplecounts", "MyMaxCount", "genus_id"]]
    crossGenus = cr.loc[cr['MyMaxCount'] < cr["highestSamplecounts"], 'genus_id'].unique()
    #print(set([x for x in cr['genus_id'].values if x in cr.loc[cr['MyMaxCount'] * 2 < cr["highestSamplecounts"], 'genus_id'].values]))
    myGenus = tax.genus.unique()
    removeGenus = []
    for x in myGenus:
        for y in crossGenus:
            if re.match(re.compile(r"(" + x + ")" + ".*"), y):
                if x not in removeGenus:
                    removeGenus.append(x)
    print(removeGenus)
    cr.sort_values(["MyMaxCount", "highestSamplecounts"], ascending=False).to_csv(out_dir+"cross.csv")
    return removeGenus



def species_prevalence(CountsDf, MetaDf, TaxDf, Taxlevel='genus', exportDir=""):

    Df = convert_df_to_taxlevel_from_myClassifier(CountsDf.copy(),TaxDf, taxlevel=Taxlevel)
    Master = Df.merge(MetaDf,how='inner',on='SampleID')
    PF = pd.DataFrame({'Taxon':Df.columns})
    print('Number of {} in the dataset initially: {}'.format(Taxlevel, PF.shape[0]))
    for name, group in Master.groupby('SampleType'):
        name  = name.replace(" ","_")
        #thisDf = group[Df.columns].astype(bool).sum(axis=0).to_frame().reset_index().rename(columns={0:"Prevalence_"+name, "index":'otuid'})
        #print(thisDf)
        PF = PF.merge(group[Df.columns].astype(bool).sum(axis=0).to_frame().reset_index().rename(columns={0:name, "index":'Taxon'}),on='Taxon',
                suffixes=("","_"+name), how='outer')
    #PF['Total'] = PF['Gut'] + PF['Oral_cavity'] + PF['Skin']
    #PS = PF.copy()
    #print(MetaDf.SampleType.unique())
    ### Estimate prevalence as percentance of samples present
    #PF[['Gut','Oral_cavity','Skin']] = PF[['Gut','Oral_cavity','Skin']].div(PF['Total'], axis=0)
    ### Estimate prevelance as percentage of each body site category seperately

    #PS['Gut'] = (PS.Gut / MetaDf.loc[MetaDf.SampleType == 'Gut'].shape[0]).round(2)
    #PS['Skin'] = (PS.Skin / MetaDf.loc[MetaDf.SampleType == 'Skin'].shape[0]).round(2)
    #PS['Oral_cavity'] = (PS.Oral_cavity / MetaDf.loc[MetaDf.SampleType == 'Oral_cavity'].shape[0]).round(2)
    #PF['otuid'] = PF['otuid'].apply(lambda x: TaxDf.loc[x,'otu'])
    #print(PF.loc[PF['otuid']=='ffc36e27c82042664a16bcd4d380b286'].to_string())
    #PF = PF.sort_values('Total',ascending=False).head(50) ### select the first most prevelant now, then maybe extract specific.
    PF = PF.melt(id_vars='Taxon', var_name='Bodysite', value_vars=['Gut','Oral_cavity','Skin'])
    #PS = PS.melt(id_vars='Taxon', var_name='Bodysite', value_vars=['Gut','Oral_cavity','Skin'])
    PF.Bodysite.replace({'Oral_cavity':'Oral cavity'},inplace=True)
    os.chdir(exportDir)
    PF.to_csv('prev_toR.tsv', sep='\t', index=None)
    #PS.to_csv('/media/christos/ssd/work/Infants/tmp/prev_sep_toR.tsv',sep='\t', index=None)
    #subprocess.run('Rscript /media/christos/ssd/work/Infants/scripts/plot_prevalence.R',shell=True)
    #subprocess.run('Rscript /media/christos/ssd/work/Infants/scripts/plot_prevalence_separately.R',shell=True)
    return (PF)
    #print(PF.loc[PF['otuid']=='ffc36e27c82042664a16bcd4d380b286'].to_string())
    #sns.barplot(x=TaxDf.loc[Df.columns,'otu'], y=Df.astype(bool).sum(axis=0))
    #plt.show()

def species_abundance_distribution(CountsDf, MetaDf, TaxDf, Taxlevel='asv', Taxlabel='asv'):
    Df = convert_df_to_taxlevel_from_myClassifier(CountsDf.copy(), TaxDf, taxlevel=Taxlevel)
    Master = Df.merge(MetaDf,how='inner',on='SampleID')
    PF = pd.DataFrame({"Taxon":Df.columns})

    for name, group in Master.groupby('SampleType'):
        name  = name.replace(" ","_")

        thisDf = group[Df.columns].mean(axis=0).to_frame().reset_index().rename(columns={0:"Mean_"+name, "index":"Taxon"})
        PF = PF.merge(thisDf,on="Taxon",suffixes=("","_"+name), how='outer')

    PF['Total'] = (PF['Mean_Gut'] + PF['Mean_Oral_cavity'] + PF['Mean_Skin']) / 3
    #PF[Taxlevel] = PF[Taxlevel].apply(lambda x: TaxDf.loc[x,Taxlabel])
    #print(PF.loc[PF['otuid']=='ffc36e27c82042664a16bcd4d380b286'].to_string())
    PF = PF.sort_values('Total',ascending=False).head(50) ### select the first most prevelant now, then maybe extract specific.
    PF = PF.melt(id_vars="Taxon",value_vars=['Mean_Gut','Mean_Oral_cavity','Mean_Skin'])
    PF.to_csv('/media/christos/ssd/work/Infants/tmp/mean_toR.tsv', sep='\t', index=None)
    subprocess.run('Rscript /media/christos/ssd/work/Infants/scripts/plot_mean_abundance.R',shell=True)
def abundance_heatmap(MasterDf, CountsDf, TaxDf, Taxlevel='asv'):
    #Meta = MasterDf[ MasterDf.columns.difference(CountsDf.columns)]
    #Df  = convert_df_to_taxlevel_from_myClassifier(CountsDf, TaxDf, taxlevel=Taxlevel).to_csv('/media/christos/ssd/work/Infants/tmp/heatmap_counts.tsv',sep='\t', index_label='SampleID')
    #Tax  = TaxDf.drop_duplicates('genus').set_index('genus').reindex(Df.columns).to_csv('/media/christos/ssd/work/Infants/tmp/heatmap_tax.tsv', sep='\t', index_label=t)
    #Meta.to_csv('/media/christos/ssd/work/Infants/tmp/heatmap_meta.tsv', sep='\t', index_label='SampleID')
    subprocess.run("Rscript /media/christos/ssd/work/Infant/scripts/plot_abundance_heatmap "
                   "/media/christos/ssd/work/Infants/tmp/CountTable.tsv "
                   "/media/christos/ssd/work/Infants/tmp/TaxonomyTable.tsv "
                   "/media/christos/ssd/work/infants/tmp/MetaTable.tsv", shell=True)

def bayes_per_bodysite( CountsDf, MetaDf, TaxDf, Taxlevel='asv', exportDir=""):
    Df = CountsDf.copy()
    MetaDf.SampleType = MetaDf.SampleType.str.replace('_', ' ')
    if Taxlevel is not None:
        Df = convert_df_to_taxlevel_from_myClassifier(Df,TaxDf, taxlevel=Taxlevel)
    MasterDf = Df.merge(MetaDf,how='inner',on='SampleID')
    bayes = pd.DataFrame({'ASV':[], 'Posterior':[], 'value':[]})
    for asv in Df.columns:
        for A, B  in itertools.permutations(MasterDf['SampleType'].unique(),2): # for each combination of body sites
            posterior = A+" | "+B
            NA = 0
            NB = 0
            Nall = 0
            Nnone = 0
            ## find populations
            for name, group in MasterDf.groupby(['Patient', 'Age']): ## for each Patient and timepoint
                if(all(x in group["SampleType"].unique() for x in [A,B])): ## if there is no missing sample, so that i compare equal number of pairs
                    if group.loc[group['SampleType'] == A, asv].values[0] > 0 and group.loc[group['SampleType'] == B, asv].values[0] > 0:
                    #if group[asv].all():
                        Nall+= 1

                    elif group.loc[group['SampleType'] == A, asv].values[0] > 0: ## implies that B = 0
                        NA += 1
                    elif group.loc[group['SampleType']== B, asv].values[0] > 0: ## implies that A = 0
                        NB += 1
                    #elif group.loc[group['SampleType'] == A, asv].item() == 0 and group.loc[group['SampleType'] == B, asv].item() == 0:
                    else:
                        Nnone += 1
            ## estimate probabilities
            Ntotal = NA + NB + Nall + Nnone
            if Nall+NA > 0: ## there must be this asv (A) to begin with
                PBonA = Nall / (Nall+ NA) ## likelihood
                PA = (Nall + NA) / Ntotal ## prior
                PB = (Nall + NB) / Ntotal ## evidence
                if PB > 0:
                    PAonB = (PBonA * PA) / PB     ## posterior
                else:
                    PAonB = np.nan ## in case there are no reads in any sample for B, PAonB cannot be defined
            else:
                PAonB = np.nan
            asv_bayes = pd.DataFrame({'ASV':[asv],"Posterior":[posterior], 'value':[PAonB]})
            bayes = bayes.append(asv_bayes)
    bayes = bayes.dropna().sort_values(['value', 'ASV'],ascending=False)
    bayes[['Site', 'Condition']]= bayes['Posterior'].str.split("|", expand=True).rename(columns={0:'Site',1:'Condition'})
    os.chdir(exportDir)
    bayes.to_csv('bayes.tsv', sep='\t', index=False)
    return  bayes

def prob_of_sharing(Df, MetaDf):
    MetaDf.SampleType = MetaDf.SampleType.str.replace('_', ' ')
    MasterDf = Df.merge(MetaDf,how='inner',on='SampleID')
    shared = pd.DataFrame({'ASV':[], 'Sites':[], 'value':[]})
    for asv in Df.columns:
        for A, B  in itertools.combinations(MasterDf['SampleType'].unique(),2): # for each combination of body sites
            Sites = A+" | "+B
            AndB = 0
            all = 0
            ## find populations
            for name, group in MasterDf.groupby(['Patient', 'Age']): ## for each Patient and timepoint
                #if(all(x in group["SampleType"].unique() for x in [A,B])): ## if there is no missing sample, so that i compare equal number of pairs
                if A in group.SampleType.values and B in group.SampleType.values:
                    if group.loc[group['SampleType'] == A, asv].values[0] > 0 and group.loc[group['SampleType'] == B, asv].values[0] > 0:
                        AndB += 1
                    all += 1
            prob_shared = AndB / all
            asv_shared = pd.DataFrame({'ASV':[asv],"Sites":[Sites], 'value':[prob_shared]})
            shared = shared.append(asv_shared)
    #shared = shared.dropna().sort_values(['value', 'ASV'],ascending=False)
    shared = shared.loc[shared.value > 0]
    #shared[['Site', 'Condition']]= shared['Sites'].str.split("|", expand=True).rename(columns={0:'Site',1:'Condition'})
    #shared.to_csv('/media/christos/ssd/work/Infants/tmp/bayes.tsv', sep='\t', index=False)
    return  shared


def odds_ratio(Df, MetaDf):
    MetaDf.SampleType = MetaDf.SampleType.str.replace('_', ' ')
    MasterDf = Df.merge(MetaDf,how='inner',on='SampleID')
    odds_df = pd.DataFrame({'OddsRatio': [], "p-value": [], "Sites": []})

    for A, B in itertools.combinations(MasterDf['SampleType'].unique(), 2):  # for each combination of body sites
        ### make a contigeny table for each ASV's presense in body site A and B --> numpy array
        ct = np.array([[0, 0], [0, 0]])  ## [[AandB, AnotB],[BnotA, notAnotB]]
        for asv in Df.columns:
            for n, g in MasterDf.groupby(['Patient', 'Timepoint']):  ## for each patient and timepoint (happenning in parallel)
                if A in g.SampleType.values and B in g.SampleType.values:
                    if g.loc[g['SampleType'] == A, asv].values[0] > 0  and g.loc[g['SampleType'] == B, asv].values[0] > 0:
                        ct[0,0] += 1
                    elif g.loc[g['SampleType'] == A, asv].values[0] > 0 and g.loc[g['SampleType'] == B, asv].values[0] == 0:
                        ct[0,1] += 1
                    elif g.loc[g['SampleType'] == A, asv].values[0] == 0 and g.loc[g['SampleType'] == B, asv].values[0] > 0:
                        ct[1,0] += 1
                    elif g.loc[g['SampleType'] == A, asv].values[0] == 0 and g.loc[g['SampleType'] == B, asv].values[0] == 0:
                        ct[1,1] += 1
        oddsratio, pvalue = scipy.stats.fisher_exact(ct)
    odds_df = odds_df.append(pd.DataFrame({'OddsRatio':[oddsratio],"p-value":[pvalue],"Sites":[A+"-"+B]}))
    odds_df.dropna(inplace=True)
    return odds_df

def cluster_transition_prob_per_bodysite(MasterDf, clusterNameColumn='dbClusterASV'):
    from IPython.display import display, HTML, Image, SVG

    ## probability that we transit to a different cluster given we are in another.
    def pairwise(iterable):
        "s -> (s0,s1), (s1,s2), (s2, s3), ..."
        a, b = itertools.tee(iterable)
        next(b, None)
        return zip(a, b)
    transition = pd.DataFrame()
    for b , bg in MasterDf.groupby('SampleType'):
        print(b)
        for c1, c2 in itertools.product(MasterDf[clusterNameColumn].unique(), repeat=2): ## all possible cluster combinations (sequences.forward and reverse are different)
            Nt1 = 0
            Nt2 = 0
            Nt1t2 = 0
            Nnone = 0
            for n , g in bg.groupby('Patient'):
                if g.shape[0] > 1: ## there must be at least two timepoints
                    g = g.sort_values('Timepoint').reset_index()
                    for t1, t2 in pairwise(g.index):
                        if g.loc[t1, clusterNameColumn] == c1 and g.loc[t2, clusterNameColumn] == c2:
                            Nt1t2 += 1
                        elif g.loc[t2, clusterNameColumn] == c2:
                            Nt2 += 1
                        elif g.loc[t1, clusterNameColumn] == c1:
                            Nt1 += 1
                        else:
                            Nnone += 1
            Ntotal = Nt1 + Nt2 + Nt1t2 + Nnone
            if Nt2 + Nt1t2 > 0:
                Pt1ont2 = Nt1t2 / (Nt1t2 + Nt2) ## likelihood
                Pt2 = (Nt1t2 + Nt2) / Ntotal ## prior
                Pt1 = (Nt1t2 + Nt1) / Ntotal ## evidence
                if Pt1 > 0:
                    Pt2ont1 = (Pt1ont2 * Pt2) / Pt1
                else:

                    Pt2ont1 = np.nan
            else:
                Pt2ont1 = np.nan
            transition = transition.append(pd.DataFrame({"From":[c1], 'To':[c2],'Posterior':[Pt2ont1], 'SampleType':[b]}),ignore_index=True)
    return  transition.sort_values(['SampleType','From','To'])
def cluster_transition_prob(MasterDf, clusterNameColumn='dbClusterASV'):
    from IPython.display import display, HTML, Image, SVG
    def pairwise(iterable):
        "s -> (s0,s1), (s1,s2), (s2, s3), ..."
        a, b = itertools.tee(iterable)
        next(b, None)
        return zip(a, b)
    transition = pd.DataFrame()
    for c1, c2 in itertools.product(MasterDf[clusterNameColumn].unique(), repeat=2): ## all possible cluster combinations (sequences.forward and reverse are different)
        Nt1 = 0
        Nt2 = 0
        Nt1t2 = 0
        Nnone = 0
        for n , g in MasterDf.groupby(['Patient','SampleType']):
            if g.shape[0] > 1: ## there must be at least two timepoints
                g = g.sort_values('Timepoint').reset_index()
                for t1, t2 in pairwise(g.index):
                    if g.loc[t1, clusterNameColumn] == c1 and g.loc[t2, clusterNameColumn] == c2:
                        Nt1t2 += 1
                    elif g.loc[t2, clusterNameColumn] == c2:
                        Nt2 += 1
                    elif g.loc[t1, clusterNameColumn] == c1:
                        Nt1 += 1
                    else:
                        Nnone += 1
        Ntotal = Nt1 + Nt2 + Nt1t2 + Nnone
        if Nt2 + Nt1t2 > 0:
            Pt1ont2 = Nt1t2 / (Nt1t2 + Nt2) ## likelihood
            Pt2 = (Nt1t2 + Nt2) / Ntotal ## prior
            Pt1 = (Nt1t2 + Nt1) / Ntotal ## evidence
            if Pt1 > 0:
                Pt2ont1 = (Pt1ont2 * Pt2) / Pt1
            else:
                Pt2ont1 = np.nan
        else:
            Pt2ont1 = np.nan
        transition = transition.append(pd.DataFrame({"From":[c1], 'To':[c2],'Posterior':[Pt2ont1]}),ignore_index=True)
    return transition




def plot_bayes(BayesDf):

    subprocess.run('Rscript plot_bodysite_detection_probabilities.R', shell=True)
def make_bodysite_covariance_matrix(CountsDf, MetaDf, Iterations=100, AsvPrev=20):
   import scipy.stats as stats

   def make_null_distr(Df, IterNum):
       """Returns a list with permutations matrices"""
       ## shuffle the OTU columns:
       permutations = []
       for i in range(IterNum):
           ## instead of countsDf use maybe the g
           cntMatrix = Df.values.astype(float)
           cntMatrixT = cntMatrix.T
           for col in cntMatrixT:
               col = np.random.shuffle(col)
           shuffledMatrix = cntMatrixT.T
           shuffledspearmanMatrix, shuffle_p_value = stats.spearmanr(shuffledMatrix.astype(float))
           permutations.append(shuffledMatrix)
       return permutations

   merged = CountsDf.merge(MetaDf, on='SampleID', how='inner')
   #merged = pd.merge(CountsDf, MetaDf, right_index=True, left_index=True)
   newDf = pd.DataFrame({"Source":[], "Target":[],"Correlation":[], "p-value":[], "Sign":[], "BodySite":[]})
   for n, g in merged.groupby("SampleType"):
       """Select ASVs based on their prevelance"""
       g = g[CountsDf.columns].loc[:, g.astype(bool).sum(axis=0).ge(AsvPrev)] ## i remove any absent OTUs from this bodysite
       print("Number of ASVs in the dataset: {}".format(g.shape[1]))
       OtuColumns = g.columns

       ## here maybe only use ASVs that are in at least 20% of samples.
       spearmanMatrix, p_value = stats.spearmanr(g.values.astype(float))
       spearmanMatrix = np.tril(spearmanMatrix, k=-1)
       permutedPvalues = make_null_distr(g, Iterations)

       for r in range(len(spearmanMatrix)):
           for c in range(len(spearmanMatrix[r])):
               source = OtuColumns[c]
       #        #source = tax.loc[source,"genus_id"]
               target = OtuColumns[r]
       #        #target = tax.loc[target,"genus_id"]
               corr = spearmanMatrix[r,c]
               pvalue  = p_value[r,c]

               ## find new empirical p-value..
               ### from all permuted matrices extract the respective corr value for this patricular cell and sort all values (null distribution)
               nulldistr  = sorted([per[r,c] for per in permutedPvalues])
               myRank = 1
               for i in range(len(nulldistr)): ## go thougth the sorted corr values from null, starting from the lowest value and stop if corr value greater than a given value
                   if corr > nulldistr[i]:
                       myRank += 1

               if len(nulldistr)/2 > myRank: ### check from both sides of the null distribution
                   empiricalPvalue = myRank / len(nulldistr)
               else:
                   empiricalPvalue = (len(nulldistr) - myRank) / myRank

               if corr < 0:
                   sign = "+"
               else:
                   sign = "-"

               #corr = abs(corr)
               if abs(corr) > 0.1:
                   thisDf  = pd.DataFrame({"Source":[source], "Target":[target],"Correlation":[np.round(corr,2)], "p-value":[empiricalPvalue], "Sign":[sign], "BodySite":[n]})
                   newDf = newDf.append(thisDf, sort=True)

   #newDf = newDf.assign(Index=[i for i in range(len(network.index))])
   #newDf.set_index("Index", inplace=True)
   #import statsmodels.sandbox.stats.multicomp as multicomp
   #newDf = newDf.assign(CorPvalue=multicomp.multipletests(newDf["p-value"], alpha=0.05, method='bonferroni')[1])
   return newDf

def transition_prob_from_master(master):
    import itertools
    def pairwise(iterable):
        "s -> (s0,s1), (s1,s2), (s2, s3), ..."
        a, b = itertools.tee(iterable)
        next(b, None)
        return zip(a, b)
    print('\n Estimate transition probabilities...')
    tr = pd.DataFrame({'site':[],'From':[],'To':[], "transition":[]})
    for n , g in master.groupby(["SampleType","Patient"]): ### for each body site and Patient
        g = g.sort_values('Timepoint').reset_index()
        #for c1, c2 in itertools.product(master["ClusterASV"].unique(), repeat=2): ### for all cluster combinations
        for t1,t2 in pairwise(g.index):
            s1 = g.loc[t1]
            s2 = g.loc[t2]
            print("{}:  ({}){} --> ({}){}".format(s1.SampleType, s1.Timepoint,s1.ClusterASV,s2.Timepoint, s2.ClusterASV))
            c1 = g.loc[t1,'ClusterASV']
            c2 = g.loc[t2,'ClusterASV']
            tr = tr.append(pd.DataFrame({'site':[n[0]],'From':[c1],'To':[c2], "transition":[c1+"-->"+c2]}))
    pr = pd.DataFrame()
    for n, g in tr.groupby(["site",'From']):
        pr =pr.append(pd.DataFrame({'site':[n[0] for i in g.To.value_counts().index],'From':[n[1] for i in g.To.value_counts().index],'To':g.To.value_counts().index, 'Prob':(g['To'].value_counts() / g['To'].value_counts().sum()).values}))
    pr['From'] = pr['From'].apply(lambda x: re.match(re.compile(r'(ASV_\d+:|)(\w).*'), x).group(2) + "C")
    pr['To'] = pr['To'].apply(lambda x: re.match(re.compile(r'(ASV_\d+:|)(\w).*'), x).group(2) + "C")
    pr.sort_values(["site",'From','To'],inplace=True)
    return pr

def cluster_freq_per_bodysite(MasterDf):
    df = MasterDf
    df['SampleCountPerSite'] = df.groupby(['SampleType'])['ClusterASV'].transform('count')
    vc = df.groupby(['SampleType', "ClusterASV"]).size().reset_index().rename(columns={0: "ClusterCountPerSite"})
    df = df.drop_duplicates(['SampleType', 'ClusterASV'])[["SampleType", 'SampleCountPerSite', "ClusterASV"]].merge(vc, on=['SampleType','ClusterASV'],how='left')
    df['ClusterFreqPerSite'] = (df.ClusterCountPerSite / df.SampleCountPerSite).apply(lambda x: round(x, 3))
    df['ClusterASV'] = df['ClusterASV'].apply(lambda x: re.match(re.compile(r'(ASV_\d+:|)(\w).*'), x).group(2) + "C")
    return(df)



def make_cluster_covariance_matrix(CountsDf, MasterDf, Iterations=100, AsvPrev=20):
   import scipy.stats as stats

   def make_null_distr(Df, IterNum):
       """Returns a list with permutations matrices"""
       ## shuffle the OTU columns:
       permutations = []
       for i in range(IterNum):
           ## instead of countsDf use maybe the g
           cntMatrix = Df.values.astype(float)
           cntMatrixT = cntMatrix.T
           for col in cntMatrixT:
               col = np.random.shuffle(col)
           shuffledMatrix = cntMatrixT.T
           shuffledspearmanMatrix, shuffle_p_value = stats.spearmanr(shuffledMatrix.astype(float))
           permutations.append(shuffledMatrix)
       return permutations

   merged = CountsDf.merge(MasterDf[['Cluster','dbCluster','ClusterASV','dbClusterASV']], on='SampleID',how='inner')
   newDf = pd.DataFrame({"Source":[], "Target":[],"Correlation":[], "p-value":[], "Sign":[], "BodySite":[]})
   for n, g in merged.groupby("dbClusterASV"):

       """Select ASVs based on their prevelance"""
       g = g[CountsDf.columns].loc[:, g.astype(bool).sum(axis=0).ge(AsvPrev)] ## i remove any absent OTUs from this bodysite
       print("Number of ASVs in the dataset: {}".format(g.shape[1]))
       OtuColumns = g.columns

       ## here maybe only use ASVs that are in at least 20% of samples.
       spearmanMatrix, p_value = stats.spearmanr(g.values.astype(float))
       spearmanMatrix = np.tril(spearmanMatrix, k=-1)
       permutedPvalues = make_null_distr(g, Iterations)

       for r in range(len(spearmanMatrix)):
           for c in range(len(spearmanMatrix[r])):
               source = OtuColumns[c]
       #        #source = tax.loc[source,"genus_id"]
               target = OtuColumns[r]
       #        #target = tax.loc[target,"genus_id"]
               corr = spearmanMatrix[r,c]
               pvalue  = p_value[r,c]

               ## find new empirical p-value..
               ### from all permuted matrices extract the respective corr value for this patricular cell and sort all values (null distribution)
               nulldistr  = sorted([per[r,c] for per in permutedPvalues])
               myRank = 1
               for i in range(len(nulldistr)): ## go thougth the sorted corr values from null, starting from the lowest value and stop if corr value greater than a given value
                   if corr > nulldistr[i]:
                       myRank += 1

               if len(nulldistr)/2 > myRank: ### check from both sides of the null distribution
                   empiricalPvalue = myRank / len(nulldistr)
               else:
                   empiricalPvalue = (len(nulldistr) - myRank) / myRank

               if corr < 0:
                   sign = "+"
               else:
                   sign = "-"

               if abs(corr) >= 0.1:
                   thisDf  = pd.DataFrame({"Source":[source], "Target":[target],"Correlation":[np.round(corr,2)], "p-value":[empiricalPvalue], "Sign":[sign], "dbCluster":[n]})
                   newDf = newDf.append(thisDf)

   #newDf = newDf.assign(Index=[i for i in range(len(network.index))])
   #newDf.set_index("Index", inplace=True)
   #import statsmodels.sandbox.stats.multicomp as multicomp
   #newDf = newDf.assign(CorPvalue=multicomp.multipletests(newDf["p-value"], alpha=0.05, method='bonferroni')[1])
   return newDf.loc[newDf['p-value'] <= 0.05]


def prepare_fasta_for_iqtree(fastaFile, tax):
    fasta = list(SeqIO.parse(fastaFile, 'fasta'))
    for rec in fasta:
        rec.id = tax.loc[rec.id, "ASVg"]
        rec.description = ""
    SeqIO.write(fasta, fastaFile, 'fasta')


def prepare_counts_for_ecol_inference(Df):
    Df = Df.rename(columns={x:x.replace(":","_") for x in Df.columns})
    Df = Df.rename(columns={x:x.replace("/","_") for x in Df.columns})
    return Df