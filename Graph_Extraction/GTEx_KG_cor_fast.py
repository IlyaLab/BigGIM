import pandas as pd
from scipy import stats
import time
import sys
import _thread
import numpy as np

'''
#Usage: python GTEx_KG_cor.py Breast 25000 30001 1000 &
# Breast tissue type in GTEx
# 25000: the first gene considered in the analysis
# 30000: the last gene considered in the analysis

#Example shell script
Tissue='Bladder'
python3 GTEx_KG_cor_fast.py $Tissue 0 5001 1000 &
sleep 250s
python3 GTEx_KG_cor_fast.py $Tissue 5000 10001 1000 &
sleep 350s
python3 GTEx_KG_cor_fast.py $Tissue 10000 15001 1000 &
sleep 350s
python3 GTEx_KG_cor_fast.py $Tissue 15000 20001 1000 &
sleep 450s
python3 GTEx_KG_cor_fast.py $Tissue 20000 30000 1000 &


#Description: This script is used to get the significantly highly correlated gene pairs 
for different tissue type using the GTEx dataset. 

Souce data: TPM for genes.
Filtering: Any gene with TPM count greater than 0.5 for at least 20 samples and 10% of the total samples are considered. 
Calculation: Spearman correlation is used.
Result: Only gene pairs with abosulute correlation coefficent greater than 0.5 and p-value smaller than 0.05 were outputed to the result table.
'''

def Filter_by_TPM(Expr_list):
    count = 0
    Remain = False
    for value in Expr_list:
        if value > 0.5:
            count = count + 1

    if count > 20 and (count > 0.1 * len(Expr_list)):
        Remain = True
    return(Remain)
    

def read_data_new(sample_select):
    Expr_file = "../Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"
    fin = open(Expr_file)
    count = 0
    gene_list =[]
    dic_expr = {}
    #dic_ensemble_gene = []
    num_samples_matched = 0
    for line in fin.readlines():
        count = count + 1
        
        if count == 3:
            samples = line.strip().split('\t')
            #print(samples)
            col_sele = []
            col_curr = 0
            for s in samples:

                if s in sample_select:
                    col_sele.append(col_curr)
                    num_samples_matched = num_samples_matched + 1
                col_curr = col_curr + 1
            print("Number of samples with tpm: "+ str(num_samples_matched))
        if count > 3:
            words = line.strip().split('\t')
            #ensemble = words[0].split('.')[0]
            gene = words[1]
            #dic_ensemble_gene[ensemble] = gene
            expr = []
            #check_gene = 0
            for i in col_sele:
                expr.append(float(words[i]))

            if Filter_by_TPM(expr) == True:
                dic_expr[gene] = expr
                gene_list.append(gene)
        
    fin.close()
    print("Number of genes: " + str(len(gene_list)))
    return(dic_expr, gene_list)

def rank_simple(vector):
    return sorted(range(len(vector)), key=vector.__getitem__)

def rankdata(a):
    n = len(a)
    ivec=rank_simple(a)
    svec=[a[rank] for rank in ivec]
    sumranks = 0
    dupcount = 0
    newarray = [0]*n
    for i in range(0,n):
        sumranks += i
        dupcount += 1
        if i==n-1 or svec[i] != svec[i+1]:
            averank = sumranks / float(dupcount) + 1
            for j in range(i-dupcount+1,i+1):
                newarray[ivec[j]] = averank
            sumranks = 0
            dupcount = 0
    return newarray

def rp_spearman_batch(x_ranks, y_ranks):
    dist_array = x_ranks - y_ranks
    d2 = np.einsum("ij,ij->i", dist_array, dist_array, optimize=True)
    n = x_ranks.shape[1]
    r = 1 - 6 * d2 / (n*(n*n-1))
    fac = (r + 1) * (1 - r)
    t = r * np.sqrt((n-2)/fac)
    t = np.where(fac <= 0.0, 0.0, t)
    return (r, n, t)



#def rp_spearman_batch_torch(x_ranks, y_ranks):
#    x_ranks = torch.tensor(x_ranks).cuda()
#    y_ranks = torch.tensor(y_ranks).cuda()
#    dist_array = x_ranks - y_ranks
#    d2 = torch.einsum("ij,ij->i", dist_array, dist_array)
#    n = x_ranks.shape[1]
#    r = 1 - 6 * d2 / (n*(n*n-1))
#    fac = (r + 1) * (1 - r)
#    t = r * torch.sqrt((n-2)/fac)
#    t = torch.where(fac <= 0.0, 0.0, t)
#    r = r.cpu()
#    t = t.cpu()
#    return (r, n, t)

def rp_spearman(x_rank,y_rank):
    #handling NA values
    #x_rank = rankdata(x)
    #y_rank = rankdata(y)
    n = len(x_rank)
    
    dist_array = np.array(x_rank) - np.array(y_rank)
    
    d2 = np.dot(dist_array, dist_array)
    
    #for i in range(0,len(x)):
    #    d2 = (x_rank[i] - y_rank[i]) * (x_rank[i] - y_rank[i]) + d2
        
    r = 1 - 6 * d2 /(n*(n*n -1))
    
    fac = (r + 1) * (1 - r)
    if fac > 0:
        t = r * np.sqrt((n - 2)/fac)
        
        p = stats.t.sf(abs(t), df=n-2)*2 # two sided test
        
    else:
        p = 0
        
    return(r,p)

#'''
#def measure_cor_from_states(dic_expr1,genes_for_query, tissue_type):
#    fou = open("result/"+tissue_type+str(genes_for_query)+".csv",'w')
#    #print((sample_select))
#    fou.write("Gene1,Gene2,rho_spearman,pvalue\n")
#
#    dic_result = {}
#    gene1 = genes_for_query
#    for gene2 in list(dic_expr.keys()):
#        if gene1 != gene2:
#            if gene1<gene2:
#                rho, pval = stats.spearmanr(dic_expr1[gene1], dic_expr1[gene2])
#                if pval < 0.05 :
#                    if (rho > 0.5) or (rho < -0.5):
#                        dic_result[gene1+','+ gene2] = [rho, pval]
#
#    for i in dic_result:
#        fou.write(i + "," +str(dic_result[i][0])+ ","+ str(dic_result[i][1]) +"\n")
#
#    fou.close()
#    return()
#'''

def measure_cor_spearman_batch(dic_expr_rank, genes_for_query1, genes_for_query2, tissue_type):
    #fou = open("result/"+tissue_type_spearman_cor+".csv",'w')
    ##fou.write("Gene1,Gene2,rho_spearman,pvalue\n")
    gene1_list = []
    #ensemble_list1= []
    gene2_list = []
    #ensemble_list2 = []
    r_list = []
    p_list = []
    count = 0

    g1_list = []
    g2_list = []

    def batch_calculate():
        batch_size = len(g1_list)
        if batch_size == 0:
            return
        rank1 = [dic_expr_rank[g] for g in g1_list]
        rank2 = [dic_expr_rank[g] for g in g2_list]
        t_start = time.time()
        r_arr, n, t_arr = rp_spearman_batch(np.array(rank1), np.array(rank2))
        #print(f"rp_spearman_batch batch_size = {batch_size} time = {time.time()-t_start}")
        
        t_start = time.time()
        p_arr = stats.t.sf(np.abs(t_arr), df=n-2)*2
        r_idx = np.where(np.abs(r_arr) > 0.5, True, False)
        p_idx = np.where(p_arr < 0.05, True, False)
        idx = np.logical_and(r_idx, p_idx)
        idx = np.where(idx)
        r_list.extend(np.round(r_arr[idx], 3).tolist())
        p_list.extend(p_arr[idx].tolist())
        gene1_list.extend(np.array(g1_list)[idx].tolist())
        gene2_list.extend(np.array(g2_list)[idx].tolist())
        #for i in range(r_arr.shape[0]):
        #    r = r_arr[i].item()
        #    t = t_arr[i].item()
        #    p = p_arr[i].item()
        #    if abs(r) > 0.5 and p < 0.05:
        #        r_list.append(round(r,3))
        #        p_list.append(p)
        #        gene1_list.append(g1_list[i])
        #        gene2_list.append(g2_list[i])
        #print(f"p_value batch_size = {batch_size} time = {time.time()-t_start}")
        g1_list.clear()
        g2_list.clear()

    for gene1 in genes_for_query1:
        for gene2 in genes_for_query2:
            if gene1 < gene2:
                g1_list.append(gene1)
                g2_list.append(gene2)
                if len(g1_list) >= 50000:
                    batch_calculate()
    batch_calculate()

    result = pd.DataFrame({"Gene1": gene1_list,
                            "Gene2": gene2_list,
                            "rho":r_list,
                            "pvalue":p_list,
                            "sample_size":[len(dic_expr[gene2])] * len(p_list),
                            "Tissue_type":[tissue_type] *len(p_list)
                            })

    #result.to_csv("result/"+tissue_type_spearman_cor+".csv",'w')
    return(result)


def measure_cor_spearman(dic_expr_rank, genes_for_query1, genes_for_query2, tissue_type):
    #fou = open("result/"+tissue_type_spearman_cor+".csv",'w')
    ##fou.write("Gene1,Gene2,rho_spearman,pvalue\n")
    gene1_list = []
    #ensemble_list1= []
    gene2_list = []
    #ensemble_list2 = []
    r_list = []
    p_list = []
    count = 0
    for gene1 in genes_for_query1:
        for gene2 in genes_for_query2:
            if gene1 < gene2:
                r,p = rp_spearman(dic_expr_rank[gene1],dic_expr_rank[gene2])
                if abs(r) > 0.5 and p < 0.05:
                    gene1_list.append(gene1)
                    gene2_list.append(gene2)
                    r_list.append(round(r,3))
                    p_list.append(p)
    result = pd.DataFrame({"Gene1": gene1_list, 
                            "Gene2": gene2_list,
                            "rho":r_list,
                            "pvalue":p_list,
                            "sample_size":[len(dic_expr[gene2])] * len(p_list),
                            "Tissue_type":[tissue_type] *len(p_list)
                            })

    #result.to_csv("result/"+tissue_type_spearman_cor+".csv",'w')
    return(result)

### ----------------------------------------Start main function ------------------------------
start_time = time.time()

tissue_type = sys.argv[1]

samples = pd.read_csv("../Data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", sep = "\t")
sample_select = set(samples.loc[samples['SMTS'] == tissue_type]['SAMPID'])
print("Number of samples : " + str(len(sample_select)))

print("--- %s seconds ---" % (time.time() - start_time))
start_time = time.time()

#### read data
dic_expr, gene_list_all = read_data_new(sample_select)
print("Finish reading!")

print("--- %s seconds ---" % (time.time() - start_time))
start_time = time.time()

#### rank data
dic_expr_rank = {}
for gene in list(dic_expr.keys()):
    dic_expr_rank[gene] = rankdata(dic_expr[gene])

print("Finish rank")
print("--- %s seconds ---" % (time.time() - start_time))
start_time = time.time()

start = int(sys.argv[2])
internel = int(sys.argv[4])
end = start + internel
All_gene_num = len(dic_expr_rank)
Stop_num = int(sys.argv[3]) 

#if int(sys.argv[3]) > All_gene_num:
#    Stop_num = All_gene_num

while end < Stop_num:
    if end < All_gene_num:
        gene_list1 = gene_list_all[start:end]
        result = measure_cor_spearman_batch(dic_expr_rank, gene_list1, gene_list_all, tissue_type)
        result.to_csv(tissue_type+"_spearman_cor_" + str(start)+ "_"+str(end) +".csv", index=False)
        start = end
        end = start + internel
    else:
        gene_list1 = gene_list_all[start:All_gene_num]
        result = measure_cor_spearman_batch(dic_expr_rank, gene_list1, gene_list_all, tissue_type)
        result.to_csv(tissue_type+"_spearman_cor_" + str(start)+ "_"+str(end) +".csv", index=False)
        break

    #result = measure_cor_spearman(dic_expr_rank, gene_list1, gene_list_all, tissue_type)
    #result.to_csv(tissue_type+"_spearman_cor_" + str(start)+ "_"+str(end) +".csv")
    #start = end
    #end = start + internel
 
print("--- %s seconds ---" % (time.time() - start_time))

