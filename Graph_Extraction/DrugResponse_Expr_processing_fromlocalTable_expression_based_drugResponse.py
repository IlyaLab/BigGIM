#!/usr/bin/env python
# coding: utf-8

# In[1]:



try:
    import statsmodels.stats.multitest as multi   
    print("module 'statsmodels' is installed")
except ModuleNotFoundError:
    get_ipython().system('pip install statsmodels')
    import statsmodels.stats.multitest as multi

from scipy.stats import stats

import pandas as pd
import numpy as np
import sys
Data_dir = "/Users/guangrong/Documents/GitHub_project/drug_response_kp/Data/GDSC/" 


# In[70]:


def Cohen_dist(x,y):
    n1 = len(x)
    n2 = len(y)
    s = np.sqrt(((n1 - 1)*(np.std(x))*(np.std(x)) + (n2 - 1) * (np.std(y)) * (np.std(y))) / (n1 + n2 -2))
    d = (np.mean(x) - np.mean(y)) / s
    return(d)


# In[2]:


Expr = pd.read_csv(Data_dir+"Cell_line_RMA_proc_basalExp.txt", sep = '\t')


# In[27]:


Expr.index = Expr['GENE_SYMBOLS']


# In[28]:


Expr_df = Expr.iloc[0:17736,2:1020]


# In[3]:


CL_anno = pd.read_excel(Data_dir + "Cell_Lines_Details.xlsx")


# In[5]:


dic_cancerType_cls = {}
for cancerType in set(CL_anno['Cancer Type\n(matching TCGA label)']):
    if cancerType != 'NA':
        dic_cancerType_cls[cancerType] = set()
        COSMIC_ids = CL_anno.loc[CL_anno['Cancer Type\n(matching TCGA label)'] == cancerType]['COSMIC identifier']
        for ids in COSMIC_ids:
            dic_cancerType_cls[cancerType].add('DATA.'+str(int(ids)))
    
    


# In[9]:


DR_table = pd.read_excel(Data_dir+"GDSC1_fitted_dose_response_25Feb20.xlsx")


# In[ ]:


CL_ID = []
for CL in DR_table['COSMIC_ID']:
    CL_ID.append('DATA.'+str(CL))
DR_table['CL_ID'] = CL_ID


# In[ ]:





# In[ ]:



Tumor_sele = sys.argv[1]
result_df = pd.DataFrame()
Gene_list = list(Expr_df.index)

LIHC_DR = DR_table.loc[DR_table['CL_ID'].isin(dic_cancerType_cls[Tumor_sele])]

for Drug_sele in list(set(DR_table['DRUG_NAME'])):
    LIHC_DR_Drug_sele = LIHC_DR.loc[LIHC_DR['DRUG_NAME'] == Drug_sele]

    LIHC_DR_Drug_sele_effect = LIHC_DR_Drug_sele.loc[LIHC_DR_Drug_sele['LN_IC50']<np.log(0.5)]
    LIHC_DR_Drug_sele = LIHC_DR_Drug_sele.loc[LIHC_DR_Drug_sele['CL_ID'].isin(Expr_df.columns)]

    if LIHC_DR_Drug_sele_effect.shape[0] >= 3:
        print(Drug_sele)
        
        sensitive_threshold = np.quantile(LIHC_DR_Drug_sele["LN_IC50"],0.25)
        resistance_threshold = np.quantile(LIHC_DR_Drug_sele["LN_IC50"],0.75)
        CL_sensitive  = list(LIHC_DR_Drug_sele.loc[LIHC_DR_Drug_sele['LN_IC50']<sensitive_threshold]['CL_ID'].values)
        CL_resistance = list(LIHC_DR_Drug_sele.loc[LIHC_DR_Drug_sele['LN_IC50']>resistance_threshold]['CL_ID'].values)

        Expr_dr_sen = Expr_df[CL_sensitive]
        Expr_dr_res = Expr_df[CL_resistance]

        p_list = []
        se_list = []

        for Gene in Gene_list:
            Expr_sen = Expr_dr_sen.loc[[Gene],:].values[0] #12
            Expr_res = Expr_dr_res.loc[[Gene],:].values[0]

            p_list.append(stats.ttest_ind(Expr_sen, Expr_res)[1])
            se_list.append(Cohen_dist(Expr_sen, Expr_res))

        FDR_List= multi.multipletests(p_list, alpha=0.05, method='fdr_bh', is_sorted=False)[1]

        result = pd.DataFrame({'Gene':Gene_list,
                                            'SE':se_list,
                                            "p_ttest":p_list, 
                                            "FDR":FDR_List, 
                                           })

        result_sig = result.loc[result['FDR'] < 0.25]
        result_sig.loc[:,['Drug']] = Drug_sele
        result_sig.loc[:,['disease']] = Tumor_sele

        result_df = pd.concat([result_df, result_sig])
result_df.to_csv("/Users/guangrong/Documents/GitHub_project/drug_response_kp/Results/" + Tumor_sele+"_expr_drugResponse.csv")


