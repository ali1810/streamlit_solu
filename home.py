

def app():
	
      st.write('**Type SMILES below**')
      SMILES = st.text_input('then press predict button', value ="CC(=O)OC1=CC=CC=C1C(=O)O")
#style = st.selectbox('Chemical structure',['stick','ball-and-stick','line','cross','sphere'])
#spin = st.checkbox('Spin', value = False)	
      #col1, col2, col3 = st.columns([10,2,11.5])
      #with col1:	
	#  st.write("   2 D Structure of the smiles  ")
      #with col2:
	#  st.write("")
      #with col3:
       #   st.write(" 3 D Structure  of the smiles")
 #       st.write("""Use mouse pointer to rotate the structure""")

      prop=pcp.get_properties([ 'MolecularWeight'], SMILES, 'smiles')
      x = list(map(lambda x: x["CID"], prop)) 
      y=x[0]
    #print(y) 
      x = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/PNG?image_size=300x300"
      url=(x % y)
#print(url)
      img = Image.open(urlopen(url))
      mol = Chem.MolFromSmiles(SMILES)
      mol = Chem.AddHs(mol)
      AllChem.EmbedMolecule(mol)
      mblock = Chem.MolToMolBlock(mol)
      xyzview = py3Dmol.view(width=300,height=300)
    #xyzview = py3Dmol.view(query=′pdb:1A2C′)
      xyzview.addModel(mblock,'mol')
      xyzview.setStyle({'model': -1}, {"cartoon": {'color': 'spectrum'}})
      style = 'stick'
#bcolor = st.sidebar.color_picker('Pick background Color', '#0C0C0B')
#xyzview.spin(True)
#if spin:
 #   xyzview.spin(True)
#else:
 #    xyzview.spin(False)
    #xyzview.setStyle({'sphere':{}})
      xyzview.setBackgroundColor('#EAE5E5')
      xyzview.zoomTo()
      xyzview.setStyle({style:{'color':'spectrum'}})
    


      col1, mid, col2 = st.columns([15,2.5,15])
#col1, mid, col2 = st.columns([15,0.5,15])
      with col1:
          st.image(img, use_column_width=False)
      with col2:
          showmol(xyzview,height=300,width=400) 
      #if st.button("predict"):
     # if st.button("Predict"):
	#      st.write("This is some content that should remain on the page.")        	      
          #page1()
      #        st.write("work in Progress") 
#def page1():      
