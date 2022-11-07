import streamlit as st
from multiapp import MultiApp
from apps import lipidsearch4, lipidsearch5, lipidxplorer # import modules 

app = MultiApp()

st.markdown("""
# LipidCruncher

In the context of lipidomics, the output of mass spectrometry is the relative abundance of the lipid species that make up the sample under study. 
This output is in the form of a spectrum in which the peaks represent an identified lipid species and the area underneath each peak reperesnts the relative 
abundance of the corresponding lipid species. There are two pieces of software that turn this spectrum into a lipidomics dataset: LipidSearch and LipidXplorer.
LipidCruncher is a web app that allows the user to perform lipidomics analysis on the LipidSearch and LipidXplorer datasets. 

![Lipidomics Workflow](fig1.png)

""")

# Add all your application here
app.add_app("LipidSearch 4.1 Module", lipidsearch4.app)
app.add_app("LipidSearch 5.0 Module", lipidsearch5.app)
app.add_app("LipidXplorer Module", lipidxplorer.app)
# The main app
app.run()
