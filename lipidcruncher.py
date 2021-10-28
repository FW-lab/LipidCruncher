import streamlit as st
from multiapp import MultiApp
from apps import home, quality_control, analysis # import your app modules here

app = MultiApp()

st.markdown("""
# LipidCruncher

LipidCruncher is a web app built to perform lipidomics analysis on LipidSearch and LipidXplorer datasets. 

""")

# Add all your application here
app.add_app("Home", home.app)
app.add_app("Quality Control Module", quality_control.app)
app.add_app("Analysis & Hypothesis Testing Module", analysis.app)
# The main app
app.run()
