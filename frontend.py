import streamlit as st
import requests
import json

import backend



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
        st.title('Optional Transformations')
        options = []
        if st.checkbox('Log-transform matrix'):
            options.append("logtransform")
        if st.checkbox('Convert to z-scores'):
           options.append("zscores")

           
  
    

    # Main App Contents
    st.title("SaVanT (Signature Visualization Tool)")
    st.text("Visualize molecular signatures in the context of gene expression matrices")
    if st.button("Generate Test Heatmap"):
      backend.constructHeatMapvalueMatrix()
    if st.button("Generate Heatmap"):
        #st.text("test")
        if selectAll:
           #constructHeatMapFromCategory('All', '', '')
           backend.constructHeatMapvalueMatrix() #revise
        else:
          backend.constructHeatMapFromCategory(group, category, signatures_selected, options)
    else:
            st.text("Upload a matrix or choose one from the drop down menu...")
            st.text("Example: ")
            st.video("https://www.youtube.com/watch?v=bVZ5Ki7aR4o")
       


if __name__ == "__main__":
    main()