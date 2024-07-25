import pandas as pd
import streamlit as st
import seaborn as sns
import numpy as np
import plotly.express as px
import scipy.stats as stats
from statsmodels.stats.multitest import fdrcorrection
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import dash_bio

#optional transformations
def LogTransformMatrix(user_matrix_df):
   #log transform values of user-uploaded matrix
   geneNames = user_matrix_df.iloc[:, 0]
   temp = user_matrix_df.iloc[:, 1:]
   log_temp = np.log(temp)
   log_matrix = pd.concat([geneNames, log_temp], axis=1)
   return log_matrix

def convertToZscore(sig_sample_df):
    #calculate average and standard deviation for each signature
    for i, row in enumerate(sig_sample_df.values):
        avg = row.mean()
        sd = row.std()
    #convert each value in signature-sample matrix to a z-score
        for j, value in enumerate(row):
          z = (value - avg) / sd
          sig_sample_df.iloc[i, j] = z
    return sig_sample_df

def DiffFromMean(user_matrix_df):
   #Transform values to difference from mean (for each gene) 
   #Note: not sure if this logic is correct
   print(user_matrix_df)
   geneNames = user_matrix_df.iloc[:, 0]
   temp = user_matrix_df.iloc[:, 1:]
   for i, row in enumerate(temp.values):
        avg = row.mean()
        for j, value in enumerate(row):
          new = (value-avg)
          temp.iloc[i, j] = new
   transformed_matrix = pd.concat([geneNames, temp], axis=1)
   print(transformed_matrix)
   return transformed_matrix

def ConvertToRanks(user_matrix_df):
   #Convert matrix values to ranks
   geneNames = user_matrix_df.iloc[:, 0]
   temp = user_matrix_df.iloc[:, 1:]
   ranked_df = temp.rank(axis=1)
   ranked_matrix = pd.concat([geneNames, ranked_df], axis=1)
   return ranked_matrix

def constructHeatMapvalueMatrix():
  signature_dict = SignatureToGeneSymbols("SaVanT", "", "")
  gene_to_sample_value_dict = GeneSymbolsToSampleValue([])
  signature_to_sample_sum = {}

  numberOfSamples = len(next(iter(gene_to_sample_value_dict.values())))

  for signature in signature_dict:
    runningSumPerSample = []
    for sample in range(numberOfSamples):
      runningSum = 0
      length = 0
      for gene in signature_dict[signature]:
        if gene in gene_to_sample_value_dict:
          runningSum += gene_to_sample_value_dict[gene][sample]
        else:
          runningSum += 0
        length += 1
      runningSumPerSample.append(float(runningSum / length))
    signature_to_sample_sum[signature] = runningSumPerSample
  print(len(signature_to_sample_sum))
  for key, value in list(signature_to_sample_sum.items())[:5]:
        print(f"{key}: {value}")
  heatMapDF = pd.DataFrame(signature_to_sample_sum)
  heatMapDF = heatMapDF.transpose() #rotate heatmap
  ax = sns.heatmap(heatMapDF, cmap='coolwarm', annot=False, fmt=".2f", cbar = 1)
  ax.xaxis.tick_top() #moves y-axis to top
  st.set_option('deprecation.showPyplotGlobalUse', False) #gets rid of Pyplot warning
  st.pyplot()

def SignatureToGeneSymbols(group, category, selected):
  if group == "Enrichr":
    signature_dict = {}
    # converts signature matrix to a dataframe making it easier to work with
    for i in range(len(category)):
      signature_matrix_path = 'files/Enrichr/' + category[i] + '.txt'
      with open(signature_matrix_path, 'r') as file:
        lines = file.readlines()
        data = [line.strip().split('\t') for line in lines]
        matrix_df = pd.DataFrame(data)

      #creates dictionary out of each signature set and then adds this to big dictionary of signatures
      temp_dict = matrix_df.set_index(0).transpose().to_dict('list')
      signature_dict.update(temp_dict)
    print("number of sigs: ", len(signature_dict.keys()))
  else: 
    signature_matrix_path = 'files/SaVanT_Signatures_Release01.tab.txt'
    matrix_df = pd.read_csv(signature_matrix_path, delimiter='\t', header=None)#, nrows=20)
  # drop null values
  
  # takes in dataframe and converts it into a hashmap that maps the signature to its corresponding genes
    signature_dict = matrix_df.set_index(0).transpose().to_dict('list')
  return signature_dict

def GeneSymbolsToSampleValue(options):
  gene_matrix_path = 'files/SaVanT_ExampleMatrix.txt'
  delimeter = '\t'
  gene_df = pd.read_csv(gene_matrix_path, delimiter=delimeter, header=None, skiprows=[0,1])

  #if user selects to log-transform data
  if "logtransform" in options:
    gene_df = LogTransformMatrix(gene_df)

  #if user selects to convert matrix values to ranks
  if "ranks" in options:
     gene_df = ConvertToRanks(gene_df)

  #if user selects to transform to difference from mean
  if "delta" in options:
     gene_df = DiffFromMean(gene_df)
  
  gene_to_sample_value_dict = gene_df.set_index(0).transpose().to_dict('list')
  return gene_to_sample_value_dict

def constructHeatMapFromCategory(group, category, signature, options):
  signature_dict = SignatureToGeneSymbols(group, category, signature)
  gene_to_sample_value_dict = GeneSymbolsToSampleValue(options)
  signature_to_sample_sum = {}
  group_list= sampleToGroup()

  #if user does not select specific signatures, all the signatures in the selected category will be used
  if group == 'Enrichr' and signature == []:
    for i in range(len(category)):
        signature_matrix_path = 'files/Enrichr/' + category[i] + '.txt'
        with open(signature_matrix_path, 'r') as file:
            lines = file.readlines()
            data = [line.strip().split('\t') for line in lines]
            matrix_df = pd.DataFrame(data)
            for value in matrix_df.iloc[:, 0]:
              signature.append(value)

  for sig in signature:
     sampleaverages=[]
     for sample in range(7):
        sum =0 
        length = 0
        for gene in signature_dict[sig]:
           if gene in gene_to_sample_value_dict:
              sum += gene_to_sample_value_dict[gene][sample]
              length+=1
        sampleaverages.append(float(sum/length))  #for each gene in a signature, calculate average, and do this for every sample
     signature_to_sample_sum[sig] = sampleaverages   #key is signature, value is array of each sample's avg

  
  #print('Signature to sample sum: ', list(signature_to_sample_sum.values())[0])
  #print('Signature to sample sum: ', list(signature_to_sample_sum.values())[1])
  heatMapDF = pd.DataFrame(signature_to_sample_sum, index=[1,2,3,4,5,6,7]) #index is there to fix a dataframe error
  heatMapDF = heatMapDF.transpose() #rotate heatmap


  pVals=anovaTest(group_list, signature_to_sample_sum) #PVALS change based on order selected rn
 
  heatMapColors = ['rgb(237, 229, 207)', 'rgb(224, 194, 162)', 'rgb(211, 156, 131)', 'rgb(193, 118, 111)', 'rgb(166, 84, 97)', 'rgb(129, 55, 83)', 'rgb(84, 31, 63)'] #same as px.Brwnyl

  #if user selects to convert to zscore
  if "zscores" in options:
    convertToZscore(heatMapDF)
    heatMapColors = [[0, '#0000FF'],[0.5, '#FFFFFF'],[1.0, '#FF0000']]
  
  if "threshold" in options:
     #go through the og heatmap, find which labels to remove, update and remove from copy
     #need to remove pVals
     count= 0
     labels_to_remove = []
     while count < len(pVals):
        #print(heatMapDF.index)
        if pVals[count][0] > 0.5:
           labels_to_remove.append(heatMapDF.index[count])
           print(labels_to_remove)
        count+=1
     print(pVals)
     DF_updated= heatMapDF.drop(labels=labels_to_remove)
     fig = px.imshow(DF_updated, color_continuous_scale= heatMapColors)
     fig.update_layout(margin=dict(l=300,r=100,b=100,t=100,pad=4))
     fig.update_traces(text=pVals)
     fig.update_traces(hovertemplate='Signature: %{y}<br>Sample: %{x}<br>Avg Exp: %{z}<br>P Value: %{text}<extra></extra>')
     fig.show()

  if "cluster" in options:
     num_sigs = len(signature_dict)
     columns = list(heatMapDF.columns.values)
     rows = list(heatMapDF.index)
     fig2 = dash_bio.Clustergram(data = heatMapDF, row_labels=rows, column_labels=columns, color_map= heatMapColors, height = num_sigs*12, width = 850, center_values = False)
     fig2.update_traces(text=pVals)
     fig2.update_traces(hovertemplate='Signature: %{y}<br>Sample: %{x}<br>Avg Exp: %{z}<br>P Value: %{text}<extra></extra>')
     fig2.show()
  elif "threshold" not in options:
     fig = px.imshow(heatMapDF, color_continuous_scale= heatMapColors)
     fig.update_layout(margin=dict(l=300,r=100,b=100,t=100,pad=4))
     fig.update_traces(text=pVals)
     fig.update_traces(hovertemplate='Signature: %{y}<br>Sample: %{x}<br>Avg Exp: %{z}<br>P Value: %{text}<extra></extra>')
     fig.show()



def anovaTest(group_list, sigsample_dict):
   group1_avgs = []
   group2_avgs = []
   pVals=[]
   corrected_pVals=[]
   for i in sigsample_dict:
    sigValues = sigsample_dict[i]
    for i in range(len(group_list)):
      if group_list[i] == 1:
        group1_avgs.append(sigValues[i])
      else:
        group2_avgs.append(sigValues[i])
    p= stats.f_oneway(group1_avgs, group2_avgs).pvalue
    pVals.append(p)
   rejectedarr, correctedarr= fdrcorrection(pVals, alpha=0.05)
   for i in correctedarr:
      corrected_pVals.append([i])
   return corrected_pVals

def sampleToGroup(): #returns an array of group numbers, that correspond to sample numbers
   group_path = 'files/SaVanT_ExampleMatrix.txt'
   delimeter = '\t'
   groups = pd.read_csv(group_path, delimiter=delimeter, header=None, skiprows=[0], nrows=1)
   group_list = list(groups.iloc[0].values)
   ga = group_list.pop(0)#skips savantgroup column placeholder
   return group_list