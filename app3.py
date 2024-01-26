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
    smiles = st.sidebar.text_input('then press predict button', value ="CC(=O)OC1=CC=CC=C1C(=O)O")
    prop=pcp.get_properties([ 'MolecularWeight'], SMILES, 'smiles')
    x = list(map(lambda x: x["CID"], prop))
    y=x[0]
    #print(y)
    x = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/PNG?image_size=400x300"
    url=(x % y)
#print(url)
    img = Image.open(urlopen(url))
    # img=smiles_to_img(smiles)
#st.write("a logo and text next to eachother")
    col1, mid, col2 = st.columns([15,0.5,15])
    with col1:
           st.image(img, use_column_width=False)
    #with col2:
     #      blk=makeblock(smiles)
      #     render_mol(blk)
	
    
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
    st.write("This is the page for project details")
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
    st.write("This is the content of Institute and contact details .")

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
