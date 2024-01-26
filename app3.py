import streamlit as st
from PIL import Image



def smiles_to_img(SMILES):
    prop=pcp.get_properties([ 'MolecularWeight'], SMILES, 'smiles')
    x = list(map(lambda x: x["CID"], prop))
    y=x[0]
    #print(y)
    x = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/PNG?image_size=400x300"
    url=(x % y)
#print(url)
    img = Image.open(urlopen(url))
    return img 

def render_mol(xyz):
    xyzview = py3Dmol.view(width=400,height=300)
    #xyzview = py3Dmol.view(query=′pdb:1A2C′)
    xyzview.addModel(xyz,'mol')
    xyzview.setStyle({'model': -1}, {"cartoon": {'color': 'spectrum'}})
    #bcolor = st.sidebar.color_picker('Pick background Color', '#0C0C0B')
    style = st.sidebar.selectbox('Chemical structure',['stick','line','cross','sphere'])
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
    showmol(xyzview,height=300,width=400) 
col1, mid, col2 = st.columns([15,0.5,15])
with col1:
        st.image(img, use_column_width=False)
with col2:
        blk=makeblock(smiles)
        render_mol(blk)


def page1():
	
	
    
    #def smiles_to_img(SMILES):
         st.sidebar.write('**Type SMILES below**')

    ## Read SMILES input
    #SMILES_input = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    #\nCC(=O)OC1=CC=CC=C1C(=O)O"
    #SMILES_input = " "
         smiles1 = st.sidebar.text_input('then press predict button', value ="CC(=O)OC1=CC=CC=C1C(=O)O")
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
    st.title("Page 2")
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
    st.title("Page 3")
    st.write("This is the content of Institute and contact details .")

def main():
    st.sidebar.title("Navigation")
    selected_page = st.sidebar.radio("Go to", ["Page 1", "Page 2", "Page 3"])
    #st.sidebar.[element_name]
	
    #selected_option = st.selectbox("Choose an option", ["Option 1", "Option 2", "Option 3"])
	

    if selected_page ==   "Page 1":
        page1()
    elif selected_page == "Page 2":
        page2()
    elif selected_page == "Page 3":
        page3()
if __name__ == "__main__":
    main()
