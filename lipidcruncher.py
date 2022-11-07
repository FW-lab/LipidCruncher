import streamlit as st
from multiapp import MultiApp
from PIL import Image
from apps import lipidsearch4, lipidsearch5, lipidxplorer # import modules 

app = MultiApp()

st.markdown("""
# LipidCruncher: a web application for performing general lipidomics 
One of our missions at Farese/Walther lab is making transformative discoveries in lipid metabolism. 
We want to understand how our bodies store and channel energy. We are hypothesis driven. 
A typical study in our lab is as follows: the lab biologists come up with a hypothesis (driven by deep expertise and intuition built over many years),
then they design/run an experiment to test the hypothesis. 
They usually have a control group where nothing is changed and one or more experimental groups where something is changed 
(i.g. a certain gene is knocked out or is over expressed). They prepare multiple replicates for each group and run them though mass spectrometer. 
One of the things that they are interested in is whether there is a significant difference in the lipid profile between 
the control group and the experimental group(s).
In the context of lipidomics, the output of mass spectrometry is the relative abundance of the lipid species that make up the sample under study. 
This output is in the form of a spectrum in which the peaks represent an identified lipid species and the area underneath each peak represents 
the relative abundance of the corresponding lipid species. 
We use two pieces of software that turn this spectrum into a lipidomics dataset: LipidSearch and LipidXplorer. 
We have built LipidCruncher which is a web-based tool that allows the biologists to perform lipidomics analysis on 
their LipidSearch and LipidXplorer datasets.  

""")

image = Image.open('fig1.png')
st.image(image)

# Add all your application here
app.add_app("LipidSearch 4.1 Module", lipidsearch4.app)
app.add_app("LipidSearch 5.0 Module", lipidsearch5.app)
app.add_app("LipidXplorer Module", lipidxplorer.app)
# The main app
app.run()
