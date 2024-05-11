import pandas as pd
import streamlit as st
import seaborn as sns
import numpy as np
import plotly.express as px
import scipy.stats as stats
import matplotlib.pyplot as plt
import plotly.graph_objects as go

def constructHeatMapvalueMatrix():
  signature_dict = SignatureToGeneSymbols("SaVanT", "", "")
  gene_to_sample_value_dict = GeneSymbolsToSampleValue()
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
  # converts signature matrix to a dataframe making it easier to work with
  if group == "Enrichr":
    #creates dictionary out of each signature set and then adds this to big dictionary of signatures
    signature_dict = {}
    for i in range(len(category)):
      signature_matrix_path = 'files/Enrichr/' + category[i] + '.txt'
      with open(signature_matrix_path, 'r') as file:
        lines = file.readlines()
        data = [line.strip().split('\t') for line in lines]
        matrix_df = pd.DataFrame(data)
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


def GeneSymbolsToSampleValue():
  gene_matrix_path = 'files/SaVanT_ExampleMatrix.txt'
  delimeter = '\t'
  gene_df = pd.read_csv(gene_matrix_path, delimiter=delimeter, header=None, skiprows=[0,1])
  gene_to_sample_value_dict = gene_df.set_index(0).transpose().to_dict('list')
  return gene_to_sample_value_dict

def convertToZscore(sig_sample_df):
    #Calculate mean and standard deviation of sig_sample matrix
    avg = np.mean(sig_sample_df.values)
    sd = np.std(sig_sample_df.values)

    #transform matrix to z-scores 
    for i, row in enumerate(sig_sample_df.values):
        for j, value in enumerate(row):
          z = (value - avg) / sd
          sig_sample_df.iloc[i, j] = z

    return sig_sample_df


def constructHeatMapFromCategory(group, category, signature, zscores):
  signature_dict = SignatureToGeneSymbols(group, category, signature)
  gene_to_sample_value_dict = GeneSymbolsToSampleValue()
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
 # color= sns.color_palette("dark:seagreen", "ch:light=.5", as_cmap=True)

  pVals=anovaTest(group_list, signature_to_sample_sum)
  print(pVals)
 
  #if user selects to convert to zscore
  if zscores:
    convertToZscore(heatMapDF)

  fig = px.imshow(heatMapDF, color_continuous_scale="Brwnyl")
  fig.update_layout(margin=dict(l=300,r=100,b=100,t=100,pad=4))
  fig.update_traces(text=pVals)
  fig.update_traces(hovertemplate='Signature: %{y}<br>Sample: %{x}<br>Avg Exp: %{z}<br>P Value: %{text}<extra></extra>')
  fig.show()

  # Create subplot for additional row or col of info
  """fig = make_subplots(
    rows=2, cols=1,  # 2 rows, 1 column
    shared_xaxes=True,
    vertical_spacing=0.1,  # Adjust vertical spacing between subplots
    row_heights=[0.8, 0.2]
)
  additional_row_trace = go.Scatter(x=[1, 2, 3, 4], y=[1, 2, 3, 4], mode='markers', marker=dict(color='red', size=10))
  heatmap_trace = go.Heatmap(z=heatMapDF)
  fig.add_trace(heatmap_trace, row=1, col=1)
  fig.add_trace(additional_row_trace, row=2, col=1)
  fig.show()
  """

def anovaTest(group_list, sigsample_dict):
   group1_avgs = []
   group2_avgs = []
   pVals=[]
   for i in sigsample_dict:
    sigValues = sigsample_dict[i]
    for i in range(len(group_list)):
      if group_list[i] == 1:
        group1_avgs.append(sigValues[i])
      else:
        group2_avgs.append(sigValues[i])
    p= [stats.f_oneway(group1_avgs, group2_avgs).pvalue]
    pVals.extend([p])
   return pVals

def sampleToGroup(): #returns an array of group numbers, that correspond to sample numbers
   group_path = 'files/SaVanT_ExampleMatrix.txt'
   delimeter = '\t'
   groups = pd.read_csv(group_path, delimiter=delimeter, header=None, skiprows=[0], nrows=1)
   group_list = list(groups.iloc[0].values)
   ga = group_list.pop(0)#skips savantgroup column placeholder
   return group_list