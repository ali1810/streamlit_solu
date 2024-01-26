import streamlit as st
from PIL import Image
from stmol import showmol
import py3Dmol
import re
import numpy as np 
from rdkit.Chem import AllChem
import pubchempy as pcp
import streamlit as st
import pickle
from PIL import Image
import pandas as pd
from rdkit import Chem
#from rdkit.Chem import Draw
import xgboost
from sklearn.neural_network import MLPRegressor
from sklearn.ensemble import RandomForestRegressor
import base64
import pickle
import numpy as np
import pandas as pd
from rdkit import Chem,DataStructs
from rdkit.Chem import MolFromSmiles, Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import Crippen
import streamlit as st
import base64
from streamlit_shap import st_shap
import shap
from xgboost import XGBRegressor
import xgboost as xgb
from urllib.request import urlopen
from PIL import Image
#from flask.wrappers import Request
#import threading
import requests

def page1():
    col1, col2, col3 = st.columns([10,2,11.5])

    with col1:
	    st.header("   2 D Structure of the smiles  ")
    with col2:
	    st.write("")
    with col3:
            st.header(" 3 D Structure  of the smiles")
            st.write("""Use mouse pointer to rotate the structure""")

    st.sidebar.write('**Type SMILES below**')
    SMILES = st.sidebar.text_input('then press predict button', value ="CC(=O)OC1=CC=CC=C1C(=O)O")
    prop=pcp.get_properties([ 'MolecularWeight'], SMILES, 'smiles')
    x = list(map(lambda x: x["CID"], prop))
    y=x[0]
    #print(y)
    x = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/PNG?image_size=400x300"
    url=(x % y)
#print(url)
    img = Image.open(urlopen(url))
    mol = Chem.MolFromSmiles(SMILES)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    mblock = Chem.MolToMolBlock(mol)	
    xyzview = py3Dmol.view(width=400,height=300)
    #xyzview = py3Dmol.view(query=′pdb:1A2C′)
    xyzview.addModel(mblock,'mol')
    xyzview.setStyle({'model': -1}, {"cartoon": {'color': 'spectrum'}})
    #bcolor = st.sidebar.color_picker('Pick background Color', '#0C0C0B')
    style = st.sidebar.selectbox('Chemical structure',['stick','ball-and-stick','line','cross','sphere'])
#spin = st.sidebar.checkbox('Spin', value = False)
    spin = st.sidebar.checkbox('Animation', value = True)
    xyzview.spin(True)
    if spin:
      xyzview.spin(True)
    else:
      xyzview.spin(False)
    #xyzview.setStyle({'sphere':{}})
    xyzview.setBackgroundColor('#EAE5E5')
    xyzview.zoomTo()
    xyzview.setStyle({style:{'color':'spectrum'}})
    #showmol(xyzview,height=300,width=400) 


	
    # img=smiles_to_img(smiles)
#st.write("a logo and text next to eachother")
    col1, mid, col2 = st.columns([15,0.5,15])
    with col1:
           st.image(img, use_column_width=False)
    with col2:
           showmol(xyzview,height=300,width=400) 
           #render_mol(blk)
    if st.sidebar.button('Predict'):
	 # define the rdkit moleculer object
         mol1 = Chem.MolFromSmiles(SMILES)
    
    # calculate the log octanol/water partition descriptor
         single_MolLogP = Descriptors.MolLogP(mol1)
    
    # calculate the molecular weight descriptor
         single_MolWt   = Descriptors.MolWt(mol1)
    
    # calculate of the number of rotatable bonds descriptor
         single_NumRotatableBonds = Descriptors.NumRotatableBonds(mol1)
    
    # calculate the aromatic proportion descriptor
         single_AP = getAromaticProportion(mol1)

    # Calculate ring count 
         single_RC= Descriptors.RingCount(mol1)

    # Calculate TPSA 
         single_TPSA=Descriptors.TPSA(mol1)

    # Calculate H Donors  
         single_Hdonors=Lipinski.NumHDonors(mol1)

    # Calculate saturated Rings 
         single_SR= Lipinski.NumSaturatedRings(mol1) 

    # Calculate Aliphatic rings 
         single_AR =Lipinski.NumAliphaticRings(mol1)
    
    # Calculate Hydrogen Acceptors 
         single_HA = Lipinski.NumHAcceptors(mol1)

    # Calculate Heteroatoms
         single_Heter = Lipinski.NumHeteroatoms(mol1)

         single_Max_Partial_Charge =  Descriptors.MaxPartialCharge(mol1)
         single_FP_density =  Descriptors.FpDensityMorgan1(mol1)
         single_num_valence_electrons = Descriptors.NumValenceElectrons(mol1)
         single_NHOH_count = Lipinski.NHOHCount(mol1)
         single_SP3_frac = Lipinski.FractionCSP3(mol1)
         single_SP_bonds = len(mol1.GetSubstructMatches(Chem.MolFromSmarts('[^1]')))
        
    

    # put the descriptors in a list
        rows = np.array([single_MolLogP, single_MolWt, single_NumRotatableBonds, single_AP,single_RC,
        single_TPSA,single_Hdonors,single_SR,single_AR,
        single_HA,single_Heter,single_Max_Partial_Charge,single_FP_density,single_num_valence_electrons,
        single_NHOH_count,single_SP3_frac,single_SP_bonds])
    
    # add the list to a pandas dataframe
    #single_df = pd.DataFrame(single_list).T
        baseData = np.vstack([rows])
    # rename the header columns of the dataframe
    
    #columnNames = ["MolLogP", "MolWt", "NumRotatableBonds", "AromaticProportion","Ring_Count","TPSA","H_donors","Saturated_Rings","AliphaticRings","H_Acceptors","Heteroatoms"]
        columnNames = ["MolP","MolWt", 
                   "NumRotatableBonds", "AromaticProportion"
                  ,"Ring_Count","TPSA","H_donors", "Saturated_Rings","AliphaticRings","H_Acceptors","Heteroatoms","Max_Partial_Charge",
                  "valence_electrons","FP_density","NHOH_count","SP3_frac","SP_bonds"]
 
        generated_descriptors1 = pd.DataFrame(data=baseData, columns=columnNames)
      
	    
        #generated_descriptors1= predictSingle(smiles)
        mol = Chem.MolFromSmiles(smiles)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=512)
        arr = np.zeros((0,), dtype=np.int8)
        arr1=DataStructs.ConvertToNumpyArray(fp,arr)
        arr2 =pd.DataFrame(arr)
        array = arr2.T
        trained_model= xgb.Booster()
        trained_model.load_model('models/model_xgb_95 2.bin')
        df3=pd.concat([array,generated_descriptors1],axis=1)
        df3 = xgb.DMatrix(df3)
        pred_rf1 = trained_model.predict(df3)
        pred_rf1 =  (pred_rf1-0.30)
        pred_rf2 =  np.round(pred_rf1,2)	
        mol_liter1   =  10**pred_rf1
        mol_liter2   = np.round(mol_liter1,2)
        c_name    =smiles_iupac(smiles)
        mol = Chem.MolFromSmiles(smiles)
        MolWt1     = generated_descriptors1["MolWt"]
        Gram_liter1  =(10**pred_rf1)*MolWt1
        P_sol1 =smiles_to_sol(smiles) ## calling function to get the solubility from <pubchem
        data = dict(IUPAC_Name=c_name,SMILES=smiles, Predicted_LogS=pred_rf2, 
        Mol_Liter=mol_liter2,Gram_Liter=Gram_liter1,Experiment_Solubility_PubChem=P_sol1)
        df = pd.DataFrame(data, index=[0])
        st.header('Predicted LogS values for single smiles')
        st.table(df.style.format({"Predicted_LogS": "{:.2f}","Mol_Liter":"{:.2f}","Gram_Liter":"{:.2f}"}))
    #df
    #st.write('Good Morning') #displayed when the button is clicked
        st.header('Computed molecular descriptors')
        generated_descriptors1 # Skips the dummy first item
    else:
        st.write('Note for users - 1>Enter Single smiles and click on predict button') #displayed when the button is unclicked
	
    
    #def smiles_to_img(SMILES):
    

    ## Read SMILES input
    #SMILES_input = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    #\nCC(=O)OC1=CC=CC=C1C(=O)O"
    #SMILES_input = " "

    
    #SMILES = SMILES.split('\n')
   #col1, col2, col3 = st.columns([10,2,11.5]) 
	
    #print(y)
        
    #return img 
	
	
	
	#st.title("AqSolPred: Online Solubility Prediction Tool")
    
    #st.markdown("<h1 style='text-align: left;position: fixed;  top: 0; width: 100%; color: blue; margin-top: -1; padding-top: -1;'>AqSolPred: Online Solubility Prediction Tool</h1>", unsafe_allow_html=True)
    #st.markdown("<h1 style='font-size: 34px ;position: fixed;  top: 0.1; width: 100%; color: blue; margin-top: 0.5; padding-top: 0.5;'>AqSolPred: Online Solubility Prediction Tool</h1>", unsafe_allow_html=True)
 
	
    #st.markdown("<h1 style='text-align: left; color: blue; position: fixed;margin-top: -1; padding-top: -1; width: 100%;
    #'>AqSolPred: Online Solubility Prediction Tool</h1>", unsafe_allow_html=True)
	
    #st.markdown("<h1 style='text-align: center; color: blue;margin-top: 0; padding-top: 0;>AqSolPred: Online Solubility Prediction Tool</h1>", unsafe_allow_html=True)
    #st.title("AqSolPred: Online Solubility Prediction Tool")
    #st.write("This is the content of Page 1.")        
def page2():
    #st.title("Page 2")
    #st.write("This is the page for project details")
    image = Image.open('Flow2.jpeg')
    col1, col2, col3 = st.columns([0.001,2.0,0.5])
    with col1:
	    st.write("")
    with col2:
            st.image(image, use_column_width=2)
    with col3:	
            st.write("")	
def page3():
    #st.title("Page 3")
    #st.write("This is the content of Institute and contact details .")
    st.write("""
# About AqSolPred1
Sol Prediction  is an  accurate solubility prediction model that consists consensus of 3 ML algorithms (Neural Nets, Random Forest, and XGBoost). 
AqSolPred is developed using a quality-oriented data selection method described in [1] and trained on AqSolDB [2] largest publicly available aqueous solubility dataset.
AqSolPred showed a top-performance (R2 - 0.95 )

If you are using the predictions from Sol Pred on your work, 
 please cite these papers: [1, 2]

[1] Sorkun, M. C., Koelman, J.M.V.A. & Er, S. (2021). [Pushing the limits of solubility prediction via quality-oriented data selection](https://www.cell.com/iscience/fulltext/S2589-0042(20)31158-5), iScience, 24(1), 101961.

[2] Sorkun, M. C., Khetan, A., & Er, S. (2019).  [AqSolDB, a curated reference set of aqueous solubility and 2D descriptors for a diverse set of compounds](https://www.nature.com/articles/s41597-019-0151-1). Scientific data, 6(1), 1-8.

[3] Huuskonen, J. Estimation of aqueous solubility for a diverse set of organic compounds based on molecular topology. Journal of Chemical Informationand Computer Sciences 40, 773–777 (2000).
Special thanks: 

This web app is developed based on the tutorials and the template of [DataProfessor's repository](https://github.com/dataprofessor/code/tree/master/streamlit/part7). 

** For any feedback or suggestion please write me -- mushtaq.ali@kit.edu                                                                                         
**Contact over Linkdin :** [Mushtaq Ali](www.linkedin.com/in/mushtaq-ali/)	""")

def main():
    st.sidebar.title("Navigation")
    selected_page = st.sidebar.radio("Visit for ", ["Solubility Prediction", "Project Detail", "Contact Detail"])
    
    #st.sidebar.[element_name]
	
    #selected_option = st.selectbox("Choose an option", ["Option 1", "Option 2", "Option 3"])
	

    if selected_page ==   "Solubility Prediction":
        page1()
    elif selected_page == "Project Detail":
        page2()
    elif selected_page == "Contact Detail":
        page3()
if __name__ == "__main__":
    main()
