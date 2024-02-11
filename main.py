import streamlit as st
from streamlit_option_menu import option_menu

import contact,home,project

st.set_page_config(page_title="Solubility Prediction")



class MultiApp:

    def __init__(self):
        self.apps = []

    def add_app(self, title, func):

        self.apps.append({
            "title": title,
            "function": func
        })

    def run():
        # app = st.sidebar(
        with st.sidebar:        
            app = option_menu(
                menu_title='Solubility Prediction ',
                options=['Home','Project','Contact'],
                icons=['house-fill','person-circle','info-circle-fill'],
                menu_icon='chat-text-fill',
                default_index=1,
                styles={
                    "container": {"padding": "5!important","background-color":'white'},
        "icon": {"color": "black", "font-size": "23px"}, 
        "nav-link": {"color":"black","font-size": "20px", "text-align": "left", "margin":"1px", "--hover-color": "blue"},
        "nav-link-selected": {"background-color": "#02ab21"},}
                
                )

        
        if app == "Home":
            home.app()
        if app == "Project":
            project.app()    
        if app == "Contact":
            contact.app()        
        
          
             
    run()            
         
