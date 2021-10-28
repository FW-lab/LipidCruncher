import streamlit as st

def app():
    
    st.header('Home')

    st.markdown("""
               
               This is the Homepage of LipidCruncher where you can find some general instructions before starting your analysis.
               Start by selecting the type of your dataset:
               
               """)
               
    dataset_type = st.radio('Select the type of your dataset', ['LipidSearch', 'LipidXplorer'], 0)
    
    if dataset_type == 'LipidSearch':
               
        st.info("""
            
                LipidSearch datasets must have the following columns: 
                
                *Rej*: one of the built-in filtering mechanisms of LipidSearch which either takes 0 (i.e. accepted) or 1 (i.e. rejected)
                    
                *LipidMolec*: the class that the lipid species belong to and its number of carbon atoms and double bonds  
                    
                *Class*: the class that the lipid species belong to
                
                *Calc Mass*: the calculated mass of the lipid species
                
                *BaseRt*: the retention time of the lipid species in the chromatography method 
                
                *MainArea[sc]*: Area Under the Curve (AUC) representing the relative abundance of the lipid species in the 
                control (i.e. internal standards) sample which is used for calibration purposes  
                
                *MainArea[s1], ..., MainArea[sN]*: Area Under the Curve (AUC) representing the relative abundance of the lipid species
                in samples s1 to sN where N stands for the total number of the sampels  
                    
                *APvalue[s1], ..., APvalue[sN]*: p-value representing the statistical significance of the identification of the lipid sepecies
                in samples s1 to sN
                    
                *MainGrade[s1], ..., MainGrade[sN]*: a letter representing the quality of the identification of the lipid species in samples s1 to sN 
                (letters 'A' and 'B' stand for high quality and letters 'C' and 'D' stand for low quality)
            
                """)
            
        st.info("""
            
                LipidCruncher allows the user to filter the data based on MainGrade and APvalue columns.
                
                *MainGrade filter*: LipidCruncher assigns score 1 to letters A and B and 0 to letters C and D
                and asks the user to set the minimum required MainGrade score. For example, if the user sets the MainGrade score equal to 4, it means 
                each lipid species need to have at least 4 A or B grades across all N samples in order to pass through the filter.
                
                *APvalue filter*: LipidCruncher sets the p-value threshold at 0.001 and asks the user to set the minimum required APvalue score.
                For example, if the user sets the APvalue score at 4, it means each lipid species need to have at least 4 identifications 
                across all N samples with a p-value smaller than the threshold in order to pass through the filter.  
                    
                """)
                
        st.info("""
            
                A missing value does not necessarily mean zero abundance. It often means "no information available".
                If a datapoint has too many missing values, it can not be useful in the process of statistical inference. 
                LipidCruncher allows the user to apply filters to the data to remove datapoints with too many missing values.  
                    
                """)
                
    elif dataset_type == 'LipidXplorer':
        
        st.info("""
            
                LipidXplorer datasets must have the following columns: 
                    
                *SPECIES*: the class that the lipid species belong to and its number of carbon atoms and double bonds  
                    
                *CLASS*: the class that the lipid species belong to
                
                *MASS*: the calculated mass of the lipid species
                
                *FORMULA*: the chemical formula of the lipid species 
                
                *ERROR*: the difference between theoritical and calculated mass
                
                *PRECURSORINTENSITY[s1], ..., PRECURSORINTENSITY[sN]*: Area Under the Curve (AUC) representing the relative abundance of the lipid species
                in samples s1 to sN where N stands for the total number of the sampels  
            
                """)
                
        st.info("""
                
                The only built-in filtering mechanism in LipidXplorer is the ERROR column. If ERROR > 5, the lipid species is removed, 
                otherwise it is kept.
                
                """)
                
        st.info("""
            
                A missing value does not necessarily mean zero abundance. It often means "no information available".
                If a datapoint has too many missing values, it can not be useful in the process of statistical inference. 
                LipidCruncher allows the user to apply filters to the data to remove datapoints with too many missing values.
                    
                """)

