import streamlit as st
import requests
import json

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#testing
def SignatureToGeneSymbols():
  # converts signature matrix to a dataframe making it easier to work with
  signautre_matrix_path = 'files/SaVanT_Signatures_Release01.tab.txt'
  delimiter = '\t'
  matrix_df = pd.read_csv(signautre_matrix_path, delimiter=delimiter, header=None)
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
  ax = sns.heatmap(heatMapDF, cmap='coolwarm', annot=False, fmt=".2f", cbar = False)
  ax.xaxis.tick_top() #moves y-axis to top
  st.pyplot()


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
            "Human": ["Skin Samples & Diseases ('SkinDB')", "Swindell ('WRS') Cell Types", "Th Cell Data", "Brain Samples", "Human Pertubation", "Macrophage Activation", "Human Body Atlas", "Primary Cell Atlas (Curated)", "Human Monocyte Subsets", "GTEx Tissues"]
        }

        species = st.selectbox("Choose Species", options=select_dict.keys())
        category = st.selectbox("Choose category", options=select_dict[species])

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
            signature = st.selectbox('Choose a signature', options=human_category2_dict[category])
        else:
            mouse_category_2_dict = {
                "Mouse Body Atlas": ['MBA_3T3-L1', 'MBA_adipose_brown'],
                "ImmGen": ['Stem Cells', 'B Cells']
            }
            signature = st.selectbox('Choose a signature', options=mouse_category_2_dict[category])


        

    # Main App Contents
    st.title("SaVanT (Signature Visualization Tool)")
    st.text("Visualize molecular signatures in the context of gene expression matrices")
    if st.button("Generate Heatmap"):
        st.text("test")
        constructHeatMapvalueMatrix()
    else:
            st.text("Upload a matrix or choose one from the drop down menu...")
            st.text("Example: ")
            st.video("https://www.youtube.com/watch?v=bVZ5Ki7aR4o")
       


if __name__ == "__main__":
    main()