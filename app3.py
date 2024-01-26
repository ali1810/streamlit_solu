import streamlit as st
from PIL import Image



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
	smiles1 = st.sidebar.text_input('then press predict button', value ="CC(=O)OC1=CC=CC=C1C(=O)O")
	
    
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
