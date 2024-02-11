import streamlit as st 

from streamlit_option_menu import option_menu

import contact,home,project

st.set_page_config(page_titel="Solubility Prediction",)

class MultiApp:
  def _init_(self):
    self.apps=[]
  def add_app(self,titel,function):
    seld.apps.append({ "titel":titel,"function":function})

  def run():
    with st.sidebar:
      app=option_menu(menu_titel='solubility prediction', options=['home','project','contact'],
                      icons=['house-fill','person-circle','contact-fill']
                      menu_icon ='chat-text-fill',
                      default_index=1, 
                      styles= {"container":{"padding":"5!important" ,"background-color":'black'
                                            "icon":{"color":"white","font-size","20px","text-align":"left",
                                                    "nav-link_selected":{"background-color":"#02ab21"}}
                                            )
    
