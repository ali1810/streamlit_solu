import streamlit as st

from streamlit-option-menu import option_menu

import contact,home,project

st.set_page_config(page_titel="Solubility Prediction",)

class MultiApp:
  def _init_(self):
    self.apps=[]
  def add_app(self,titel,function):
    seld.apps.append({ "titel":titel,"function":function})

  def run():
    with st.sidebar:
      app=option_menu(menu_titel='solubility prediction', options=['Home','Project','Contact'],
                      icons=['house-fill','person-circle','chat-fill'],
                      menu_icon ='chat-text-fill',
                      default_index=1, 
                      styles= {"container":{"padding":"5!important" ,"background-color":"black"},
                                            "icon":{"color":"white","font-size":"23px"},
                                             "nav-link": {"color":"white","font-size": "20px", "text-align": "left", "margin":"0px", "--hover-color": "blue"},
                                             "nav-link-selected": {"background-color": "#02ab21"},}
                

                                            )
    if app == "Home":
          home.app()
    if app == "Project":
          project.app()    
    if app == "Contact":
          contact.app()        
        #if app == 'Your Posts':
         #   your.app()
        #if app == 'about':
         #   about.app()    
             
          
             
    run()            
         
