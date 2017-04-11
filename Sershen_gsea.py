
# coding: utf-8

# In[54]:

# Program for GSEA


# In[55]:

import numpy as np
import scipy as sp
import pandas as pd
from __future__ import division


# In[56]:

#set the number of permutations, number of columns in file
num_permut=1000
# expression data parameters
no_cols=72
no_rows=2672
num_ALL=47
#pathways data parameters
cols_to_drop=15
num_non_gene_cols=2


# In[57]:

# read in the pathways

df = pd.read_table('./data/pathways.gmt', header=None)




# In[58]:

#read in expression data, a tab delimited file, including headers
dfexpr = pd.read_csv('./data/leukemia3.txt', sep='\t',header=0)


# In[59]:

# drop the non-numerical col "gene"
dfexpr_1 = dfexpr



# In[60]:

#Create numpy arr of exper values
arrayn=dfexpr_1.iloc[:,0:no_cols].values


# In[61]:

arraym=np.reshape(arrayn,(no_rows,no_cols))
arr_ALL=arraym[:,0:num_ALL-1]
arr_AML=arraym[:,num_ALL:no_cols]


# In[62]:

#function to calculate the correlations for each gene
def run_corr (mtx):
  my_array = np.zeros(no_rows)
  
  for i in xrange(0, len(mtx[:,0])):  
      miu=np.mean(mtx[i,:])
      stdev=np.std(mtx[i,:])  
      correl=(miu)/stdev  
      my_array[i]=correl

  return my_array
    
final_STON_ALL=run_corr(arr_ALL)
final_STON_AML=run_corr(arr_AML)


# In[63]:

# add on the signal to noise ratio to df
dfexpr_ALL=dfexpr_1[0:no_rows+1][0:]
dfexpr_AML=dfexpr_1[0:no_rows+1][0:]

dfexpr_ALL["ston"]=final_STON_ALL
dfexpr_AML["ston"]=final_STON_AML


# In[64]:

# remove the rows with less than 15 genes from the pathway df
df["num_NaN"]=df.isnull().sum(axis=1)


# In[65]:

#find num cols in df
df_cols=len(df.columns)


# In[66]:

# create new df with dropped cols for pathways with lt 15 genes ex index and first 2 cols
df_filtered=df[df["num_NaN"]<df_cols-(cols_to_drop+num_non_gene_cols+1)]
df_filtered.reset_index(drop=True,inplace=True)


# In[67]:

def hitormiss(datafrtosearch, searchstr):
  return searchstr in datafrtosearch[3:].values
  


# In[68]:

df_filtered.reset_index(inplace=True)


# In[69]:

headerlist = list(df_filtered)
headerlist1 = headerlist[1:]


# In[70]:

df_filtered_pivot=df_filtered.transpose()
df_filtered_pivot = df_filtered_pivot


# In[71]:

series_gene = dfexpr_ALL["gene"].values
series_ston = dfexpr_ALL["ston"].values
dfexpr_ALL_reduced = pd.DataFrame(data=series_gene,columns=['gene'])
dfexpr_AML_reduced = pd.DataFrame(data=series_gene,columns=['gene'])
dfexpr_ALL_reduced['ston'] = dfexpr_ALL["ston"].values
dfexpr_AML_reduced['ston'] = dfexpr_AML["ston"].values


# In[72]:

for i in xrange(0,len(df_filtered_pivot.columns)):
    varname = df_filtered_pivot.iloc[1,i]  
    
    dfexpr_ALL_reduced[varname]=dfexpr_ALL_reduced["gene"].apply(lambda x: x in df_filtered_pivot.iloc[2:,i].values) #hitormiss,(df_filtered, x)
    dfexpr_AML_reduced[varname]=dfexpr_AML_reduced["gene"].apply(lambda x: x in df_filtered_pivot.iloc[2:,i].values) #hitormiss,(df_filtered, x)


# In[73]:

#sort the df in descend order by ston values
dfexpr_ALL_reduced_sorted=dfexpr_ALL_reduced.sort_values(by="ston",ascending=False)
dfexpr_AML_reduced_sorted=dfexpr_AML_reduced.sort_values(by="ston",ascending=False)


# In[74]:

#to write lambda for looping
series_pathways = df_filtered_pivot.ix[1]

def calc_ES(dframe, dframe2, series_pathways):

  ES=np.zeros(len(dframe2.columns))

  for i in xrange(0,len(dframe2.columns)):
    
    series_pathw = dframe[series_pathways[i]].values
    df_calc = pd.DataFrame(data=series_pathw,columns=[series_pathways[i]])
    
    df_calc["ston"] =  dframe["ston"]
    df_calc["TrueFalse"] = df_calc[df_calc.columns[0]]
    df_calc["forTrue"] = df_calc["TrueFalse"].apply(lambda x: 1 if (x==True) else 0)
    df_calc["forFalse"] = df_calc["TrueFalse"].apply(lambda x: 1 if (x==False) else 0)
    df_calc["stonhits"] = df_calc["forTrue"]*df_calc["ston"].abs()
    df_calc["abs_ston_hits"]=df_calc["stonhits"].abs()
    df_calc["sum_abs_ston_hits"]=df_calc["abs_ston_hits"].cumsum()
    sum_ston_hits=df_calc["sum_abs_ston_hits"].iloc[-1]
    df_calc["P_hit"]=df_calc["sum_abs_ston_hits"]/sum_ston_hits.max()
    df_calc["interim"]=df_calc["forTrue"].cumsum()
    df_calc["interim2"]=1/(2672-df_calc["interim"])
    df_calc["P_miss"]=df_calc["interim2"].cumsum()
    df_calc["P_hit-P_miss"]=df_calc["P_hit"]-df_calc["P_miss"]

    max1=df_calc["P_hit-P_miss"].max()
    min1=df_calc["P_hit-P_miss"].min()
    
    ES[i]=ES_fn(min1,max1)
  return ES
    
    
def ES_fn(min3,max4):
    if abs(min3)>max4:
        return min3
    else:
        return max4
    
ES_ALL=calc_ES(dfexpr_ALL_reduced_sorted,df_filtered_pivot,series_pathways)
ES_AML=calc_ES(dfexpr_AML_reduced_sorted,df_filtered_pivot,series_pathways)

print ES_ALL.shape


# In[75]:

# random calcs
ES_ALL_ran=[]
ES_AML_ran=[]
NES_ALL=np.zeros(len(ES_ALL))
NES_AML=np.zeros(len(ES_AML))
FDR_ALL=np.zeros(len(ES_ALL))
FDR_AML=np.zeros(len(ES_AML))

for i in xrange(0,num_permut):
  temp=arraym[:, np.random.permutation(arraym.shape[1])]
  arr_ALL_ran=temp[:,0:num_ALL-1] 
  arr_AML_ran=temp[:,num_ALL:no_cols]
  final_STON_ALL_ran=run_corr(arr_ALL_ran) 
  final_STON_AML_ran=run_corr(arr_AML_ran)
  #Calculate ALL T/F matrix
  dfexpr_ALL_reduced["ston"]=final_STON_ALL_ran
  dfexpr_ALL_reduced_sorted=dfexpr_ALL_reduced.sort_values(by="ston",ascending=False)
  ES_ALL_ran=np.append(ES_ALL_ran,calc_ES(dfexpr_ALL_reduced_sorted,df_filtered_pivot,series_pathways))
  dfexpr_AML_reduced["ston"]=final_STON_AML_ran
  dfexpr_AML_reduced_sorted=dfexpr_AML_reduced.sort_values(by="ston",ascending=False)
  ES_AML_ran=np.append(ES_AML_ran,calc_ES(dfexpr_AML_reduced_sorted,df_filtered_pivot,series_pathways))


# In[76]:

#calc NES, p_val and FDR

NES_ALL=np.zeros(len(ES_ALL))
NES_AML=np.zeros(len(ES_AML))
FDR_ALL=np.zeros(len(ES_ALL))
FDR_AML=np.zeros(len(ES_AML))

p_value_ALL=np.zeros(len(ES_ALL))
p_value_AML=np.zeros(len(ES_AML))
denominator1=np.mean([i for i in ES_ALL_ran if i > 0])
denominator2=np.mean([i for i in ES_AML_ran if i > 0])
denominator3=np.mean([i for i in ES_ALL_ran if i < 0])
denominator4=np.mean([i for i in ES_AML_ran if i < 0])
NES_ALL_pi=np.zeros(len(ES_ALL_ran))
NES_AML_pi=np.zeros(len(ES_ALL_ran))

for i in xrange(len(ES_ALL)):
  if ES_ALL[i]>0:  
    p_value_ALL[i]=(ES_ALL_ran[:]>ES_ALL[i]).sum()/len(ES_ALL_ran)
    NES_ALL[i]=ES_ALL[i]/denominator1
  else: 
    p_value_ALL[i]=(ES_ALL_ran[:]<ES_ALL[i]).sum()/len(ES_ALL_ran)
    NES_ALL[i]=ES_ALL[i]/denominator3
  if ES_AML[i]>0:  
    p_value_AML[i]=(ES_ALL_ran[:]>ES_ALL[i]).sum()/len(ES_ALL_ran)
    NES_AML[i]=ES_AML[i]/denominator2
  else: 
    p_value_AML[i]=(ES_ALL_ran[:]<ES_ALL[i]).sum()/len(ES_ALL_ran)
    NES_AML[i]=ES_AML[i]/denominator4  
    
for i in xrange(len(ES_ALL_ran)):
  if ES_ALL_ran[i]>0:  
    NES_ALL_pi[i]=ES_ALL_ran[i]/denominator1   
  else: 
    NES_ALL_pi[i]=ES_ALL_ran[i]/denominator3
  if ES_AML_ran[i]>0:  
    NES_AML_pi[i]=ES_AML_ran[i]/denominator2
  else: 
    NES_AML_pi[i]=ES_AML_ran[i]/denominator4
    

for i in xrange(len(NES_AML)):
  if NES_AML[i]>0:
    numerator=(NES_AML_pi[:]>NES_AML[i]).sum()/(NES_AML_pi[:]>0).sum()
  else:
    numerator=(NES_AML_pi[:]<NES_AML[i]).sum()/(NES_AML_pi[:]<0).sum()

  if NES_AML[i]>0:
    denominator=(NES_AML[:]>ES_AML[i]).sum()/(NES_AML[:]>0).sum()
  else:
    denominator=(NES_AML[:]<ES_AML[i]).sum()/(NES_AML[:]<0).sum()

  FDR_AML[i] = numerator / denominator

    
for i in xrange(len(NES_ALL)):
  if NES_ALL[i]>0:
    numerator=(NES_ALL_pi[:]>NES_ALL[i]).sum()/(NES_ALL_pi[:]>0).sum()
  else:
    numerator=(NES_ALL_pi[:]<NES_ALL[i]).sum()/(NES_ALL_pi[:]<0).sum()

  if NES_ALL[i]>0:
    denominator=(NES_ALL[:]>ES_ALL[i]).sum()/(NES_ALL[:]>0).sum()
  else:
    denominator=(NES_ALL[:]<ES_ALL[i]).sum()/(NES_ALL[:]<0).sum()

  FDR_ALL[i] = numerator / denominator


# In[77]:

#print out the results
ALL_results=pd.DataFrame(data=series_pathways, columns=['pathway'])
ALL_results['ES']=ES_ALL
ALL_results['p_val']=p_value_ALL
ALL_results['NES']=NES_ALL
ALL_results['FDR']=FDR_ALL
ALL_results_sorted=ALL_results.sort_values(by="NES",ascending=False)
ALL_results_sorted.to_csv('ALL_results.csv')

AML_results=pd.DataFrame(data=series_pathways, columns=['pathway'])
AML_results['ES']=ES_AML
AML_results['p_val']=p_value_AML
AML_results['NES']=NES_AML
AML_results['FDR']=FDR_AML
AML_results_sorted=AML_results.sort_values(by="NES",ascending=False)
AML_results_sorted.to_csv('AML_results.csv')

