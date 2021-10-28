import streamlit as st
from PIL import Image 

def app():
    
    st.header('Analysis & Hypothesis Testing Module')

    img = Image.open("under_construction.png")
    
    st.image(img)

