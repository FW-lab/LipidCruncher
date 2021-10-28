import streamlit as st
from multiapp import MultiApp
from apps import home, quality_control, analysis # import modules 

app = MultiApp()

st.markdown("""
# LipidCruncher

In the context of lipidomics, the output of mass spectrometry is the relative abundance of the lipid species that make up the sample under study. 
This output is in the form of a spectrum in which the peaks represent an identified lipid species and the area underneath each peak reperesnts the relative 
abundance of the corresponding lipid species. There are two pieces of software that turn this spectrum into a lipidomics dataset: LipidSearch and LipidXplorer.
LipidCruncher is a web app that allows the user to perform lipidomics analysis on the LipidSearch and LipidXplorer datasets. 

""")

# Add all your application here
app.add_app("Home", home.app)
app.add_app("Quality Control Module", quality_control.app)
app.add_app("Analysis & Hypothesis Testing Module", analysis.app)
# The main app
app.run()
