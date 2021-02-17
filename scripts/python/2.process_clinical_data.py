import sys, os, re, argparse, subprocess
import pandas as pd
import numpy as np

#from functions import *

parser=argparse.ArgumentParser(description="This script will process clinical data..")
parser.add_argument("-i", dest="rootDir", help="Project directory - (should contain tmp and clinical_data directories)", required=True)
args=parser.parse_args()

rootDir = args.rootDir
if rootDir.endswith("/") is False:
    rootDir += '/'

tmpDir = rootDir+'tmp/'
cliDir = rootDir+'clinical_data/'

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
os.chdir(cliDir)
KHAZ_list = pd.read_excel('patient_list.xlsx',engine='openpyxl')
pdict = {k:v for k,v in zip(KHAZ_list['khaz'], KHAZ_list['patient'])}
asv = pd.read_csv(tmpDir+'CountTable.tsv', sep='\t', index_col="SampleID")
meta = pd.read_csv(tmpDir+'MetaTable.tsv', sep='\t', index_col='SampleID')
meta['Date'] = pd.to_datetime(meta.Date)



"""                                              Nutrition data                                         """
bila = pd.read_excel('Bilanzen.xlsx', engine='openpyxl').rename(columns={'KHAZ':'Patient'})
nut = pd.read_csv('nutrition_data_table.tsv', sep='\t',parse_dates=['Date'],dayfirst=True)
bila['Patient'] = bila['Patient'].map(pdict)
### milk composition according to Michaelsen etal (1990), Macronutrient (g/dL) and energy (kcal/dL)
milk = pd.DataFrame({'Group':['protein','fat','lactose','energy'],
                    "Median":[.9, 3.6, 7.2, 67],
                    '-2std':[.6,1.8,6.4,17],
                    '+2std':[1.4,8.9,7.6,117]})
milk.set_index('Group', inplace=True)
bila['Protein'] = bila.SummeEnteral * milk.loc['protein', 'Median'] / 10
bila['Fat'] = bila.SummeEnteral * milk.loc['fat', 'Median'] / 10
bila['Lactose'] = bila.SummeEnteral * milk.loc['lactose','Median'] /10
bila['Energy'] = bila.SummeEnteral * milk.loc['energy', "Median"] / 10 ## cause it is kcal/dl
bila = bila.groupby(["Patient", pd.Grouper(key='Date', freq='D')]).sum().reset_index()
bila['DoL'] = bila.groupby('Patient')['Date'].transform(lambda x: x-x.min()+pd.Timedelta(days=0)).dt.days
#bila["Timepoint"] = pd.cut(anti.Age, bins=[0,2,5,9,16, np.inf], labels=False).apply(lambda x: x+1 )
bila = bila[bila.Energy > 0]
#bila.to_csv('/media/christos/ssd/work/Infants/tmp/feeding.tsv',sep='\t', index=False)

meta.reset_index(inplace=True)
enter = meta.merge(nut, on=['Patient','Date'], how='left')[['Patient','Age','Ratio', 'Timepoint']]
#enter['DoL'] = enter.groupby('Patient')['Date'].transform(lambda x: x-x.min()+pd.Timedelta(days=1)).dt.days
enter['Enteral_feeding'] = enter.Ratio.astype(str).apply(lambda x: int(x.split(":")[1])) ### this is the percentage of enteral feeding received for this day
enter.drop(columns=['Ratio'],inplace=True)
meta = meta.merge(enter.drop_duplicates(['Patient','Age']), on=['Patient','Age','Timepoint'], how='left')
enteral_melt = pd.melt(enter, value_vars='Enteral_feeding', id_vars=['Patient',"Age", 'Timepoint'])
enteral_melt.to_csv(tmpDir+'enteral_toR.csv',index=False)




""""                                Ventilation data               """
venti = pd.read_excel('Ventilation.xlsx', engine='openpyxl').rename(columns={'KHAZ':'Patient'})
venti['Patient'] = venti['Patient'].map(pdict)
ventilation_types = ['spontan','non-invasiv', 'invasiv', 'Atem-unterst.']
venti['Date'] = pd.to_datetime(venti['Date'], format='%m/%d/%y')
venti = venti[['Patient','Date','VentilationType']]
venti['Age'] = venti.groupby('Patient')['Date'].transform(lambda x: x-x.min()+pd.Timedelta(days=0)).dt.days
#venti["Timepoint"] = pd.cut(venti.Age, bins=[0,2,5,9,16,np.inf], labels=False).apply(lambda x: x+1 )


"""                                Oxygen data                    """
oxy = pd.read_excel('oxygenLevel.xlsx', engine='openpyxl')
oxy.rename(columns = {'Value':'oxygen', 'KHAZ':'Patient'}, inplace=True)
oxy['Date'] = pd.to_datetime(oxy['DATUMZTP'],format='%m/%d/%y %H:%M:%S')
oxy['Patient'] = oxy['Patient'].map(pdict)

oxy = oxy.groupby(['Patient', pd.Grouper(key='Date', freq='D')]).agg({'oxygen':'std'}).reset_index()
respi = oxy.merge(venti, on=['Patient','Date'], how='outer')
### FiO2 is measured for a relatively short period compared to the oxygen and ventilation data.
### pcocess oxygen adjustment in the instrument (if needed to be used))
## two values on the table fiO2 and Insp. O2, have to seperate in the beginning cause they are measured with latency.
#### if fio2 is 21 (%), it means that the patient breaths by himself/herself
fio = pd.read_excel('FiO2.xlsx', engine='openpyxl')
fio.rename(columns={'DATUMZTP':"Date",'KHAZ':'Patient'}, inplace=True) ### rename columns from German
fio['Patient'] = fio['Patient'].map(pdict)

###################################################################################
fio2 = fio.loc[fio['shortLabel.1']== 'Insp. O2'] ### i extract only the Insp. 02 rows
fio2 = fio2.groupby(['Patient', pd.Grouper(key='Date', freq='D')]).agg({'WERT':'min'}).reset_index().rename(columns={'WERT':'FiO2'})
#fio2 = fio2.set_index('Date').resample('D').max().reset_index().rename(columns = {'WERT':'Insp_O2'})[['Date','KHAZ','Insp_O2']]
respi2 = respi.merge(fio2, on=['Patient', 'Date'], how='outer')
respi2['Age'] = respi2.groupby('Patient')['Date'].transform(lambda x: x-x.min()+pd.Timedelta(days=0)).dt.days
#respi2["Timepoint"] = pd.cut(respi2.Age, bins=[0,2,5,9,16, np.inf], labels=False).apply(lambda x: x+1 )
respi2.to_csv(tmpDir+'respi.tsv',sep='\t', index=False, na_rep='NaN')
##################################################################################

oxd = fio.loc[fio['shortLabel.1']== 'Insp. O2'].rename(columns={'WERT':'FiO2'})
oxd = oxd.groupby(['Patient', pd.Grouper(key='Date', freq='h')]).agg({'FiO2':'min'}).reset_index()
oxd.FiO2 = oxd.FiO2.astype(float)
def get_oxygen_days(col):
    if col.size > 12:
        return 1
    else:
        return 0

oxd = oxd.loc[oxd.FiO2 > 21].groupby(['Patient', pd.Grouper(key='Date',freq='D')]).agg({'FiO2': get_oxygen_days}).rename(columns={'FiO2':'OxygenDay'}).reset_index()
oxd['Age'] = oxd.groupby('Patient')['Date'].transform(lambda x: x-x.min()+pd.Timedelta(days=0)).dt.days


oxyI = pd.DataFrame()
DaysPrior = 2
for n, g in meta.drop_duplicates(['Patient','Date']).groupby(['Patient', 'Date']):
    Age = g.Age.head(1).values[0]
    pao = oxd.loc[oxd.Patient == n[0]] ### subset for each patient
    r = pd.date_range(start=pao.Date.min(), end=pao.Date.max())  ### define the time period for each patient
    pao = pao.set_index('Date').reindex(r).fillna(0.0).rename_axis('Date').reset_index()
    ## index that looks at the whole history (prior to the timepoint)
    oxygenCumIndex = pao.loc[pao.Date.le(n[1]), 'OxygenDay'].sum()
    ## index that looks at a defined time window before the timepoint (including the timepoint)
    oxygenIndex = pao.loc[pao.Date.between(n[1]-pd.DateOffset(days=DaysPrior),n[1]), 'OxygenDay'].sum()
    oxyI = oxyI.append(pd.DataFrame({'Patient':[n[0]], 'Age':[Age], 'oxygenCumIndex':[oxygenCumIndex],
                                    'oxygenPastIndex':[oxygenIndex]}), ignore_index=True)


"""                   Inflammation (IL-6, CPR)                    """
infe = pd.read_excel('Infektparameter.xlsx', engine='openpyxl').rename(columns={'KHAZ':'Patient', 'ZEITPKT':'Date'})
infe['Patient'] = infe['Patient'].map(pdict)
infe = infe.groupby(['Patient', 'Variable', pd.Grouper(key='Date',freq='D')]).agg({'Value':'max'}).reset_index()
infe['Age'] = infe.groupby('Patient')['Date'].transform(lambda x: x-x.min()+pd.Timedelta(days=0)).dt.days
#infe["Timepoint"] = pd.cut(infe.Age, bins=[0,2,5,9,16, np.inf], labels=False).apply(lambda x: x+1 )

il6 = infe.loc[infe.Variable == 'Interleukin-6 [pg/ml]']
#il6 = il6[il6.Timepoint != 5]
def define_inflammation(col):
    ## col contains information per day, will return the number of days that IL-6 was > 150 pg/ml
    return len(col[col > 150])
il6 = il6.groupby(['Patient','Age', pd.Grouper(key='Date',freq='D')]).agg({'Value':define_inflammation}).reset_index().rename(columns={'Value':'InflammationIndex'})

il6I = pd.DataFrame()
DaysPrior = 2
DaysAfter = DaysPrior
for n, g in meta.drop_duplicates(['Patient','Date']).groupby(['Patient', 'Date']):
    Age = g.Age.head(1).values[0]
    pao = il6.loc[il6.Patient == n[0]]
    r = pd.date_range(start=pao.Date.min(), end=pao.Date.max())
    pao = pao.set_index('Date').reindex(r).fillna(0.0).rename_axis('Date').reset_index()
    ## index that looks at the whole history (prior to the timepoint)
    #Il6CumIndex = pao.loc[pao.Date.le(n[1]), 'OxygenDay'].sum()
    ## index that looks at a defined time window before the timepoint
    il6PastIndex = pao.loc[pao.Date.between(n[1]-pd.DateOffset(days=DaysPrior),n[1]), 'InflammationIndex'].sum()
    il6WinIndex = pao.loc[pao.Date.between(n[1] - pd.DateOffset(days=DaysPrior), n[1] + pd.DateOffset(days=DaysAfter)), 'InflammationIndex'].sum()

    il6I = il6I.append(pd.DataFrame({'Patient':[n[0]], 'Age':[Age], 'il6PastIndex':[il6PastIndex],
                                    'il6WinIndex':[il6WinIndex]}), ignore_index=True)

#il6 = il6.groupby(['Patient','Timepoint']).agg({'Value':define_inflammation}).reset_index().rename(columns={'Value':'InflammationIndex'})
infe.to_csv(tmpDir+'infe.tsv', sep='\t', index=False)





"""                                Antibiotics                                        """
anti = pd.read_excel('Antiinfektiva.xlsx', engine='openpyxl').rename(columns={"ZEITPKT":"Date", "KHAZ":"Patient"})
anti['Patient'] = anti['Patient'].map(pdict)
antibiotics = ['Ampicillin','Gentamicin', 'Ofloxacin' ,
               'Piperacillin/Tazobactam','Meropenem', 'Vancomycin',
               'Penicillin G','Clarithromycin','Mupirocin',
               'Flucloxacillin', 'Fusidinsäure', 'Bacitracin/Neomycin','Linezolid']
print('Different units in the dataset: {}'.format(", ".join(anti['UOM'].unique())))
def extract_dosage_for_gtt(String):
    import re
    return float(re.match(re.compile('\w+\s+\(\w+\)\s(\d+).+'), String).group(1))

def extract_administration(String):
    ## both i.v and KI are intraveinously.
    import re
    if re.match(re.compile(r'.+(i.v|KI|Appl).+'), String):
        return re.match(re.compile(r'.+(i.v|KI|Appl).+'), String).group(1)
    else:
        return 'i.v'

anti.loc[anti['UOM'] == 'gtt', 'Dose'] = anti.loc[anti['UOM'] == 'gtt','Additional'].apply(extract_dosage_for_gtt) # is this right ?
anti = anti.groupby(["Patient",'Medication','UOM','Additional', pd.Grouper(key='Date', freq='D')]).mean().reset_index()
anti['Age'] = anti.groupby('Patient')['Date'].transform(lambda x: x-x.min()+pd.Timedelta(days=0)).dt.days
## extract only the antibiotics (antibacterial)
anti = anti.loc[anti.Medication.isin(antibiotics)] ## extract only anti bacterial medication
anti['Admin'] = anti['Additional'].apply(extract_administration)
## quick fix for vancomycin --> none (i assumed it is administered v.i), ask Lukas though
def change_vanco_label(row):
    if row.Medication == 'Vancomycin':
        if row.Admin == 'i.v':
            return 'Vancomycin_iv'
        else:
            return 'Vancomycin_ki'
    else:
        return row.Medication
anti.Medication = anti.apply(change_vanco_label, axis=1)
antibiotics = list(anti.Medication.unique())
DaysPrior = 2

antiI = pd.DataFrame()
for n, g in meta.drop_duplicates(['Patient','Date']).groupby(['Patient', 'Date']):
    Age = g.Age.head(1).values[0]
    pao = anti.loc[anti.Patient == n[0]]
    #r = pd.date_range(start=pao.Date.min(), end=pao.Date.max())
    #pao = pao.set_index('Date').reindex(r).fillna(0.0).rename_axis('Date').reset_index()
    antiCumIndex = pao.loc[pao.Date.between(n[1]-pd.DateOffset(days=DaysPrior),n[1]), "Medication"].size
    antiCumIndex = pao.loc[pao.Date.le(n[1]), "Medication"].size
    antiIndex = pao.loc[pao.Date.between(n[1]-pd.DateOffset(days=DaysPrior),n[1]), "Medication"].size
    antiI = antiI.append(pd.DataFrame({"Patient":[n[0]], "Age":[Age],
                                       "AntiCumIndex":[antiCumIndex],
                                       "AntiIndex":[antiIndex]}), ignore_index=True, sort=True)


def extract_antibiotic_index(col):
    ## total number of antibiotic entries for the window, normalized by the size of the window (days span)
    return (col.size) / col.drop_duplicates().size

anti.to_csv(tmpDir+'anti.tsv', sep='\t', index=False)



"""                        Blood gas analysis                        """
### relevant to be used, pH , CO2, HbF(hemoglubin), SBC ( standard bicarbonate)
bga = pd.read_excel('BGA.xlsx', engine='openpyxl')
bga.groupby(['KHAZ', 'Variable', pd.Grouper(key='Timepoint', freq='D')]).mean()
list(bga['Variable'].unique())


"""                        Export data                              """
meta = antiI.merge(meta, on = ['Patient','Age'], how='inner')
meta  = il6I.merge(meta, on = ["Patient","Age"], how='inner')
meta = oxyI.merge(meta, on = ['Patient','Age'], how='inner')
meta.set_index('SampleID', inplace=True)

selected_cols = ['DM','SEX','Age','GA', 'SampleType','BW','BPD', 'ROP','IVH', 'Late-onset-sepsis', 'LOS_bc_pos','il6WinIndex', 'AntiCumIndex', 'oxygenCumIndex', 'Timepoint', "Enteral_feeding"]
envfit  = meta[selected_cols]
#metadata_correlation(meta_df)
print(envfit.columns)
envfit.to_csv(tmpDir+'meta_envfit.tsv', sep='\t')
meta.to_csv(tmpDir+'meta_with_clinical.tsv', sep='\t')
C = envfit.corr().stack()
C = C[C.index.get_level_values(0) != C.index.get_level_values(1)].sort_values().drop_duplicates()
print(C)
