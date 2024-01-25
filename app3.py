import streamlit as st
from PIL import Image

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


def page1():
    #st.title("AqSolPred: Online Solubility Prediction Tool")
    
    #st.markdown("<h1 style='text-align: left;position: fixed;  top: 0; width: 100%; color: blue; margin-top: -1; padding-top: -1;'>AqSolPred: Online Solubility Prediction Tool</h1>", unsafe_allow_html=True)
    #st.markdown("<h1 style='font-size: 34px ;position: fixed;  top: 0.1; width: 100%; color: blue; margin-top: 0.5; padding-top: 0.5;'>AqSolPred: Online Solubility Prediction Tool</h1>", unsafe_allow_html=True)
 
	
    #st.markdown("<h1 style='text-align: left; color: blue; position: fixed;margin-top: -1; padding-top: -1; width: 100%;
    #'>AqSolPred: Online Solubility Prediction Tool</h1>", unsafe_allow_html=True)
	
    #st.markdown("<h1 style='text-align: center; color: blue;margin-top: 0; padding-top: 0;>AqSolPred: Online Solubility Prediction Tool</h1>", unsafe_allow_html=True)
    #st.title("AqSolPred: Online Solubility Prediction Tool")
    #st.write("This is the content of Page 1.")

    #st.write("This is the content of Page 1.")
    #st.set_page_config(page_title="AqSolPred: Online Solubility Prediction Tool",layout="wide")
    #st.write("""# Solibility Prediction on Aqueous Solvent """)
    
    image = Image.open('Flow2.jpeg')
    col1, col2, col3 = st.columns([0.001,2.0,0.5])
    with col1:
	    st.write("")
    with col2:
            st.image(image, use_column_width=8)
    with col3:	
            st.write("")
def page2():
    st.title("Page 2")
    st.write("This is the content of Page 2.")
def page3():
    st.title("Page 3")
    st.write("This is the content of Page 3.")

def main():
    st.sidebar.title("Navigation")
    selected_page = st.sidebar.radio("Go to", ["Page 1", "Page 2", "Page 3"])

    if selected_page == "Page 1":
        page1()
    elif selected_page == "Page 2":
        page2()
    elif selected_page == "Page 3":
        page3()
if __name__ == "__main__":
    main()
