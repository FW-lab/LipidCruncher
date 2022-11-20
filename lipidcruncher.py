import streamlit as st
from multiapp import MultiApp
from PIL import Image
from apps import lipidsearch4, lipidsearch5, lipidxplorer # import modules 

app = MultiApp()

st.markdown("""
# LipidCruncher: a web application for processing, visualizing, and analyzing lipidomics data
The advancements in biomedical engineering over the past decade has enabled the collection of lipidomics data at a scale that was not possible before. 
However, in the absence of a lipid analytics tool, biologists are forced to analyze their data manually using tools such as Excel. 
That is a tedious process with a high possibility of error. 
We present LipidCruncher, a web-based tool for processing, visualizing, and analyzing lipidomics data from LipidSearch and LipidXplorer datasets. 
LipidCruncher comes with two sub-modules. The first submodule allows the biologists to clean, filter, and explore their data, run quality checks on it, 
and detect anomalies. The second submodule provides more in-depth analysis such as volcano plots for hypothesis testing, 
saturation profile analysis and lipidomics pathway visualization.  

""")

image = Image.open('fig1.png')
st.image(image)

# Add all your application here
app.add_app("LipidSearch 4.1 Module", lipidsearch4.app)
app.add_app("LipidSearch 5.0 Module", lipidsearch5.app)
app.add_app("LipidXplorer Module", lipidxplorer.app)
# The main app
app.run()
