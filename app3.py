import streamlit as st

def page1():
    st.title("AqSolPred: Online Solubility Prediction Tool")
    #st.write("This is the content of Page 1.")
    #st.set_page_config(page_title="AqSolPred: Online Solubility Prediction Tool",layout="wide")
    st.write("""# Solibility Prediction on Aqueous Solvent """)
    image = Image.open('Flow.jpeg')
    col1, col2, col3 = st.columns([0.5,2.0,0.5])
    image = Image.open('Flow.jpeg')
    col1, col2, col3 = st.columns([0.5,2.0,0.5])
    with col1:
	st.write("")

    with col2:
	
        st.image(image, use_column_width=6)
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