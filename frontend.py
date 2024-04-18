import streamlit as st
import requests
import json

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
import scipy.stats as stats


def SignatureToGeneSymbols(species, category, selected):
  # converts signature matrix to a dataframe making it easier to work with
  if species == "Enrichr":
    #signature_matrix_path = 'files/Enrichr/' + category + '.txt'
    signature_matrix_path = 'files/Enrichr/Achilles_fitness_decrease.txt'
    with open(signature_matrix_path, 'r') as file:
      lines = file.readlines()
    data = [line.strip().split('\t') for line in lines]
    matrix_df = pd.DataFrame(data)
  else: 
    signature_matrix_path = 'files/SaVanT_Signatures_Release01.tab.txt'
    matrix_df = pd.read_csv(signature_matrix_path, delimiter='/t', header=None)#, nrows=20)
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


def constructHeatMapvalueMatrix():
  signature_dict = SignatureToGeneSymbols()
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
  st.pyplot()

def constructHeatMapFromCategory(species, category, signature):
  signature_dict = SignatureToGeneSymbols(species, category, signature)
  gene_to_sample_value_dict = GeneSymbolsToSampleValue()
  signature_to_sample_sum = {}


  for sig in signature:
     sampleaverages=[]
     for sample in range(7):
        sum =0 
        length = 0
        for gene in signature_dict[sig]:
           if gene in gene_to_sample_value_dict:
              sum += gene_to_sample_value_dict[gene][sample]
           else:
              sum += 0
           length += 1
        sampleaverages.append(float(sum/length))  #for each gene in a signature, calculate average, and do this for every sample
     signature_to_sample_sum[sig] = sampleaverages   #key is signature, value is array of each sample's avg

  
  print('Signature to sample sum: ', list(signature_to_sample_sum.values())[0])
  heatMapDF = pd.DataFrame(signature_to_sample_sum, index=[1,2,3,4,5,6,7]) #index is there to fix a dataframe error
  heatMapDF = heatMapDF.transpose() #rotate heatmap
 # color= sns.color_palette("dark:seagreen", "ch:light=.5", as_cmap=True)
  
  
  hm = sns.heatmap(heatMapDF, annot=True, fmt=".2f", cbar = 1, cmap="YlGnBu", linewidths=0.3)
  ax2= hm.twiny()
  hm.xaxis.tick_top() #moves y-axis to top
  ax2.set_xlim(hm.get_xlim())
  ax2.set_xlabel("Group Number")
  ax2.xaxis.set_label_position("bottom")
  ax2.xaxis.tick_bottom()
  #go into samples, look at row of groups
  #map the groups to signature value array
  group_list= sampleToGroup()
  #ax2.xaxis.set_ticks(group_list)
  ax2.set_xticklabels(group_list)    
  for label in ax2.get_xticklabels():   #works assuming there are 2 groups
    if label.get_text() == '1':
        label.set_color('red')
    elif label.get_text() == '2':
        label.set_color('green')
  
  st.pyplot()
  
  #fig.canvas.mpl_connect('motion_notify_event')
  sigValues = list(signature_to_sample_sum.values())[0]

  group1_avgs = []
  group2_avgs = []

  
  for i in range(len(group_list)):
     if group_list[i] == 1:
        group1_avgs.append(sigValues[i])
     else:
        group2_avgs.append(sigValues[i])

  x= stats.f_oneway(group1_avgs, group2_avgs)
  print('Anova test', x)

def anovaTest(group_list, sigsample_dict):
   #for every sig in dictionary, assign samples to groups and conduct a test
   result_dict= {} #dictionary of tuples, one for group, one for sigvalue
   #maybe create a tuple in original heatmap construct?
   #for value, group in zip(group_list, sigsample_dict array)
   #for value in sigsample_dict:
      
   return
  
def sampleToGroup(): #returns an array of group numbers, that correspond to sample numbers
   group_path = 'files/SaVanT_ExampleMatrix.txt'
   delimeter = '\t'
   groups = pd.read_csv(group_path, delimiter=delimeter, header=None, skiprows=[0], nrows=1)
   group_list = list(groups.iloc[0].values)
   ga = group_list.pop(0)#skips savantgroup column placeholder
   return group_list

def main():
    # Sidebar
    with st.sidebar:
        
        # Upload File Var
        st.title("User Upload Matrix")
        uploaded_file = st.file_uploader('Upload a Gene Expression Matrix: ', type=['txt'])
        add_file = st.button('Submit Matrix')
        upload_success = False

        # User upload 
        if add_file:
          if uploaded_file and add_file and not upload_success:
              matrix_content = uploaded_file.read()
              matrix_data = {'matrix': matrix_content.decode('utf-8')}
              headers = {'Content-Type': 'application/json'}  # Set content type to JSON
              response = requests.post('http://127.0.0.1:8000/upload_matrix/', data=json.dumps(matrix_data), headers=headers)

              if response.status_code == 200 and response.json().get('status') == 'success':
                st.success('Matrix uploaded successfully!')
                upload_success = True
          else:
              st.error('Failed to upload matrix.')
        
        # Choose Ranked Signature
        st.title('Select / Upload Signatures')
        st.title('Ranked Signaturess')
        
        select_dict = {
            "Mouse": ["Mouse Body Atlas", "ImmGen"],
            "Human": ["Skin Samples & Diseases ('SkinDB')", "Swindell ('WRS') Cell Types", "Th Cell Data", "Brain Samples", "Human Pertubation", "Macrophage Activation", "Human Body Atlas", "Primary Cell Atlas (Curated)", "Human Monocyte Subsets", "GTEx Tissues"],
            "Enrichr": ["Achilles_fitness_decrease"]
        }

        species = st.selectbox("Choose Species", options=select_dict.keys())
        category = st.multiselect("Choose category", options=select_dict[species])

        if species == "Human":
            human_category2_dict = {
                "Skin Samples & Diseases ('SkinDB')": ['Acne', 'Acute wound (0h after injury)', 'Allergic contact dermatitis'],
                "Swindell ('WRS') Cell Types": ['WRS_B_cell', 'WRS_CD138+Plasma_Cell', 'WRS_CD34+cell'],
                "Th Cell Data": ['TH_Th17', 'TH_Th1_Harvard'], 
                "Brain Samples": ['Astrocytes', 'Cortical neurons'],
                "Human Pertubation": ['MacCyto_adPBMC_IL4_6h', 'MacCyto_adPBMC_IL4_24h'],
                "Macrophage Activation": ['MA_B', 'MA_DC_imm'],
                "Human Body Atlas": ['HBA_721_B_lymphoblasts', 'HBA_Adipocyte'],
                "Primary Cell Atlas (Curated)": ['HPCA_Adipocytes', 'HPCA_B_cells'],
                "Human Monocyte Subsets": ['Classical Monocytes: CD14++CD16-', 'Intermediate Monocytes: CD14++CD16+'],
                "GTEx Tissues": ['GTEx adipose - subcutaenous', 'GTEx adipose - visceral (omentum)']
            }
            signatures= []
            for sig in category:
                  subcategories = human_category2_dict[sig]
                  signatures.extend(subcategories)
            signatures_selected = st.multiselect('Choose a signature', options=signatures)
        elif species == "Mouse":
            mouse_category_2_dict = {
                "Mouse Body Atlas": ['MBA_3T3-L1', 'MBA_adipose_brown'],
                "ImmGen": ['Stem Cells', 'B Cells']
            }
            signatures= []
            for sig in category:
                  subcategories = mouse_category_2_dict[sig]
                  signatures.extend(subcategories)
            signatures_selected = st.multiselect('Choose a signature', options=signatures)
        else:
            Enrichr_category_2_dict = {
                "Achilles_fitness_decrease": ['22RV1-prostate', '697-haematopoietic and lymphoid tissue'],
            }
            signatures= []
            for sig in category:
                  subcategories = Enrichr_category_2_dict[sig]
                  signatures.extend(subcategories)
            signatures_selected = st.multiselect('Choose a signature', options=signatures)

  
    

    # Main App Contents
    st.title("SaVanT (Signature Visualization Tool)")
    st.text("Visualize molecular signatures in the context of gene expression matrices")
    if st.button("Generate Test Heatmap"):
      constructHeatMapvalueMatrix()
    if st.button("Generate Heatmap"):
        st.text("test")
        constructHeatMapFromCategory(species, category, signatures_selected)
    else:
            st.text("Upload a matrix or choose one from the drop down menu...")
            st.text("Example: ")
            st.video("https://www.youtube.com/watch?v=bVZ5Ki7aR4o")
       


if __name__ == "__main__":
    main()