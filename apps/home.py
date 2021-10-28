import streamlit as st

def app():
    st.header('Home')

    st.markdown("""
               
               Currently, Lipid Lab only supports LipidSearch datasets.
               LipidXplorer capability will be added soon.
               
               """)
               
    st.info("""
            
            Note: LipidSearch datasets must have the following columns: 'Rej', 'LipidMolec', 'Class', 'Calc Mass', 'BaseRt', 'MainArea[c]', 'MainArea[s1]', 
                    ..., 'MainArea[sN]', APvalue[s1], ..., APvalue[sN], MainGrade[s1], ..., MainGrade[sN].
                    
            N is the total number of the samples.
            
            """)
            
    st.info("""
            
            Note: LipidSearch datasets come with a short description of the samples at the beginning of the dataset. 
            Make sure to remove that description before uploading the dataset. Follow the example below. 
            
            """)

