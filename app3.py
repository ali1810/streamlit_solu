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
import streamlit as st
#from streamlit.components.v1 import ComponentMeta

import streamlit as st

# Function to display the target page
#if selected_page ==   "Solubility Prediction":
 #       page1()
  #  elif selected_page == "Project Detail":
   #     page2()
   # elif selected_page == "Contact Detail":
    #    page3()


#def target_page():
 #   st.write("This is the target page. You can put any content here.")

# Display a clickable heading



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

   



def page2():
    #st.title("Page 2")
    #st.write("This is the page for project details")
    image = Image.open('flow5.jpeg')
    #image_path = 'flow5.jpeg'
    col1, mid, col2 = st.columns([15,0.5,5])
    with col1:
           st.image(image, use_column_width=False)
    with col2:
           #showmol(xyzview,height=300,width=400)
	    st.write("")
    col1, col2, col3 = st.columns([0.001,2,1])
    #with col1:
	#    st.write("")
    #with col2:
	#    st.image(image, use_column_width=2)
    #with col3:	
     #       st.write("")
    #st.markdown(
     #   f'<div style="display: flex; justify-content: center;">'
      #  f'<img src="{image_path}" style="width: 50%; height: auto;">'
      #  f'</div>',
       # unsafe_allow_html=True
    #)	
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
    #st.sidebar.title("Navigation")
    #selected_page = st.sidebar.radio("Visit for ", ["Solubility Prediction", "Project Detail", "Contact Detail"])
    #selected_page = 
    #st.sidebar.button(["Solubility Prediction", "Project Detail", "Contact Detail"])
    
    button1_clicked = st.sidebar.button("Solubility Prediction 1")
    button2_clicked = st.sidebar.button("Project details")
    button3_clicked = st.sidebar.button("Contact Details")
    
    if button1_clicked:
         page1()	    
           
    if button2_clicked:
         page2()
    if button3_clicked:
         page3()



	
    #selected_option = st.selectbox("Choose an option", ["Option 1", "Option 2", "Option 3"])
   # if st.sidebar.button("Solubility Prediction "):
    #    page1()
    #if st.sidebar.button("Project details "):
     #    page2()
    #if st.sidebar.button("Contact Details"):
     #    page3()	

    #if selected_page ==   "Solubility Prediction":
     #   page1()
    #elif selected_page == "Project Detail":
     #   page2()
    #elif selected_page == "Contact Detail":
     #   page3()
if __name__ == "__main__":
    main()
