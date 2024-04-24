import streamlit as st
import requests
import json

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import plotly.graph_objects as go
#import plotly.io as pio
import plotly.express as px
from plotly.subplots import make_subplots


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

def constructHeatMapFromCategory(group, category, signature):
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
           else:
              sum += 0
           length += 1
        sampleaverages.append(float(sum/length))  #for each gene in a signature, calculate average, and do this for every sample
     signature_to_sample_sum[sig] = sampleaverages   #key is signature, value is array of each sample's avg

  
  print('Signature to sample sum: ', list(signature_to_sample_sum.values())[0])
  heatMapDF = pd.DataFrame(signature_to_sample_sum, index=[1,2,3,4,5,6,7]) #index is there to fix a dataframe error
  heatMapDF = heatMapDF.transpose() #rotate heatmap
 # color= sns.color_palette("dark:seagreen", "ch:light=.5", as_cmap=True)
  
  """fig = go.Figure(data=go.Heatmap(
        z=heatMapDF.values,  # Pass DataFrame values
        x=heatMapDF.columns,  # Use DataFrame columns for x-axis
        y=heatMapDF.index,  # Use DataFrame index for y-axis
        hoverongaps=False,
        colorscale="blues" 
    ))"""
  #pio.show(fig)

  fig = px.imshow(heatMapDF, color_continuous_scale="Brwnyl")
  fig.update_layout(margin=dict(l=300,r=100,b=100,t=100,pad=4))
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
        
        selectAll = False
        # Choose Ranked Signature
        st.title('Select / Upload Signatures')
        st.title('Choose Signatures:')
        
        select_dict = {
            "Enrichr": ['ARCHS4_Cell-lines', 'ARCHS4_IDG_Coexp', 'ARCHS4_Kinases_Coexp', 'ARCHS4_TFs_Coexp', 'ARCHS4_Tissues', 'Achilles_fitness_decrease', 'Achilles_fitness_increase', 'Aging_Perturbations_from_GEO_down', 'Aging_Perturbations_from_GEO_up', 'Allen_Brain_Atlas_10x_scRNA_2021', 'Allen_Brain_Atlas_down', 'Allen_Brain_Atlas_up', 'Azimuth_2023', 'Azimuth_Cell_Types_2021', 'BioCarta_2013', 'BioCarta_2015', 'BioCarta_2016', 'BioPlanet_2019', 'BioPlex_2017', 'CCLE_Proteomics_2020', 'CORUM', 'COVID-19_Related_Gene_Sets', 'COVID-19_Related_Gene_Sets_2021', 'Cancer_Cell_Line_Encyclopedia', 'CellMarker_Augmented_2021', 'ChEA_2013', 'ChEA_2015', 'ChEA_2016', 'ChEA_2022', 'Chromosome_Location', 'Chromosome_Location_hg19', 'ClinVar_2019', 'DSigDB', 'Data_Acquisition_Method_Most_Popular_Genes', 'DepMap_WG_CRISPR_Screens_Broad_CellLines_2019', 'DepMap_WG_CRISPR_Screens_Sanger_CellLines_2019', 'Descartes_Cell_Types_and_Tissue_2021', 'Diabetes_Perturbations_GEO_2022', 'DisGeNET', 'Disease_Perturbations_from_GEO_down', 'Disease_Perturbations_from_GEO_up', 'Disease_Signatures_from_GEO_down_2014', 'Disease_Signatures_from_GEO_up_2014', 'DrugMatrix', 'Drug_Perturbations_from_GEO_2014', 'Drug_Perturbations_from_GEO_down', 'Drug_Perturbations_from_GEO_up'],
            "SaVanT signatures": ["Mouse Body Atlas", "ImmGen", "Skin Samples & Diseases ('SkinDB')", "Swindell ('WRS') Cell Types", "Th Cell Data", "Brain Samples", "Human Pertubation", "Macrophage Activation", "Human Body Atlas", "Primary Cell Atlas (Curated)", "Human Monocyte Subsets", "GTEx Tissues"]
        }

        group = st.selectbox("Choose Group", options=select_dict.keys())
        category = st.multiselect("Choose category", options=select_dict[group])

        if group == "SaVanT signatures":
           savant_category2_dict = {
                "Skin Samples & Diseases ('SkinDB')": ['Acne', 'Acute wound (0h after injury)', 'Allergic contact dermatitis'],
                "Swindell ('WRS') Cell Types": ['WRS_B_cell', 'WRS_CD138+Plasma_Cell', 'WRS_CD34+cell'],
                "Th Cell Data": ['TH_Th17', 'TH_Th1_Harvard'], 
                "Brain Samples": ['Astrocytes', 'Cortical neurons'],
                "Human Pertubation": ['MacCyto_adPBMC_IL4_6h', 'MacCyto_adPBMC_IL4_24h'],
                "Macrophage Activation": ['MA_B', 'MA_DC_imm'],
                "Human Body Atlas": ['HBA_721_B_lymphoblasts', 'HBA_Adipocyte'],
                "Primary Cell Atlas (Curated)": ['HPCA_Adipocytes', 'HPCA_B_cells'],
                "Human Monocyte Subsets": ['Classical Monocytes: CD14++CD16-', 'Intermediate Monocytes: CD14++CD16+'],
                "GTEx Tissues": ['GTEx adipose - subcutaenous', 'GTEx adipose - visceral (omentum)'],
                "Mouse Body Atlas": ['MBA_3T3-L1', 'MBA_adipose_brown'],
                "ImmGen": ['Stem Cells', 'B Cells']
            }
           signatures= []
           for sig in category:
                  subcategories = savant_category2_dict[sig]
                  signatures.extend(subcategories)
           signatures_selected = st.multiselect('Choose a signature', options=signatures)
        elif group == "Enrichr":
            Enrichr_category_2_dict = {
                'Achilles_fitness_decrease': ['22RV1-prostate', '697-haematopoietic and lymphoid tissue', '786O-kidney', 'A1207-central nervous system', 'A172-central nervous system', 'A204-soft tissue', 'A2058-skin', 'A549-lung', 'A673-bone', 'ACHN-kidney', 'AGS-stomach', 'AM38-central nervous system', 'AML193-haematopoietic and lymphoid tissue', 'ASPC1-pancreas', 'BT20-breast', 'BT474-breast', 'BXPC3-pancreas', 'C2BBE1-large intestine', 'C32-skin', 'CADOES1-bone', 'CAL120-breast', 'CAL51-breast', 'CALU1-lung', 'CAOV3-ovary', 'CAOV4-ovary', 'CAS1-central nervous system', 'CFPAC1-pancreas', 'CH157MN-central nervous system', 'COLO205-large intestine', 'COLO704-ovary', 'COLO741-skin', 'COLO783-skin', 'CORL23-lung', 'COV362-ovary', 'COV434-ovary', 'COV504-ovary', 'COV644-ovary', 'DBTRG05MG-central nervous system', 'DKMG-central nervous system', 'DLD1-large intestine', 'EFE184-endometrium', 'EFM19-breast', 'EFO21-ovary', 'EFO27-ovary', 'EW8-bone', 'EWS502-bone', 'F36P-haematopoietic and lymphoid tissue', 'GB1-central nervous system', 'GP2D-large intestine', 'HCC1187-breast', 'HCC1395-breast', 'HCC1954-breast', 'HCC2218-breast', 'HCC2814-lung', 'HCC364-lung', 'HCC44-lung', 'HCC70-breast', 'HCC827-lung', 'HCC827GR5-lung', 'HCT116-large intestine', 'HEC1A-endometrium', 'HEYA8-ovary', 'HL60-haematopoietic and lymphoid tissue', 'HLF-liver', 'HNT34-haematopoietic and lymphoid tissue', 'HPAC-pancreas', 'HPAFII-pancreas', 'HS683-central nervous system', 'HS766T-pancreas', 'HS944T-skin'],
                'Achilles_fitness_increase': ['22RV1-prostate', '697-haematopoietic and lymphoid tissue', '786O-kidney', 'A1207-central nervous system'],
                'ARCHS4_Cell-lines': [],
                'ARCHS4_IDG_Coexp': [],
                'ARCHS4_Kinases_Coexp': [],
                'ARCHS4_TFs_Coexp': [],
                'ARCHS4_Tissues': [],
                'Aging_Perturbations_from_GEO_down': [],
                'Aging_Perturbations_from_GEO_up': [],
                'Allen_Brain_Atlas_10x_scRNA_2021': [],
                'Allen_Brain_Atlas_down': [],
                'Allen_Brain_Atlas_up': [],
                'Azimuth_2023': [],
                'Azimuth_Cell_Types_2021': [],
                'BioCarta_2013': [],
                'BioCarta_2015': [],
                'BioCarta_2016': [],
                'BioPlanet_2019': [],
                'BioPlex_2017': [],
                'CCLE_Proteomics_2020': [],
                'CORUM': [],
                'COVID-19_Related_Gene_Sets': [],
                'COVID-19_Related_Gene_Sets_2021': [],
                'Cancer_Cell_Line_Encyclopedia': [],
                'CellMarker_Augmented_2021': [],
                'ChEA_2013': [],
                'ChEA_2015': [],
                'ChEA_2016': [],
                'ChEA_2022': [],
                'Chromosome_Location': [],
                'Chromosome_Location_hg19': [],
                'ClinVar_2019': [],
                'DSigDB': [],
                'Data_Acquisition_Method_Most_Popular_Genes': [],
                'DepMap_WG_CRISPR_Screens_Broad_CellLines_2019': [],
                'DepMap_WG_CRISPR_Screens_Sanger_CellLines_2019': [],
                'Descartes_Cell_Types_and_Tissue_2021': [],
                'Diabetes_Perturbations_GEO_2022': [],
                'DisGeNET': [],
                'Disease_Perturbations_from_GEO_down': [],
                'Disease_Perturbations_from_GEO_up': [],
                'Disease_Signatures_from_GEO_down_2014': [],
                'Disease_Signatures_from_GEO_up_2014': [],
                'DrugMatrix': [],
                'Drug_Perturbations_from_GEO_2014': [],
                'Drug_Perturbations_from_GEO_down': [],
                'Drug_Perturbations_from_GEO_up': [],
            }
            signatures= []
            for sig in category:
                  subcategories = Enrichr_category_2_dict[sig]
                  signatures.extend(subcategories)
            signatures_selected = st.multiselect('Choose a signature', options=signatures)
        st.title('Or Select All Signatures:')
        if st.checkbox('Select All'):
           selectAll = True
  
    

    # Main App Contents
    st.title("SaVanT (Signature Visualization Tool)")
    st.text("Visualize molecular signatures in the context of gene expression matrices")
    if st.button("Generate Test Heatmap"):
      constructHeatMapvalueMatrix()
    if st.button("Generate Heatmap"):
        #st.text("test")
        if selectAll:
           #constructHeatMapFromCategory('All', '', '')
           constructHeatMapvalueMatrix() #revise
        else:
          constructHeatMapFromCategory(group, category, signatures_selected)
    else:
            st.text("Upload a matrix or choose one from the drop down menu...")
            st.text("Example: ")
            st.video("https://www.youtube.com/watch?v=bVZ5Ki7aR4o")
       


if __name__ == "__main__":
    main()