import csv
import streamlit as st
import pandas as pd
import numpy as np
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, HoverTool
import matplotlib.pyplot as plt
import seaborn as sns 
from sklearn.preprocessing import scale # Data scaling
from sklearn import decomposition # PCA
from scipy import stats

##########################################################################################################################################
# functions used in the main code of the app 

# building LipidXplorer dataframe (only positive or negative mode)
def build_single_lipidxplorer_df(lipid_xplorer):
        
        with open(lipid_xplorer.name, "wb") as f:
            
            f.write(lipid_xplorer.getbuffer()) # saves the file in the app directory
        
        with open(lipid_xplorer.name, newline='') as f:
            
            reader = csv.reader(f)
            
            row1 = next(reader)
            
            if 'SPECIES' in row1[0]:
                
                df = pd.read_csv(lipid_xplorer.name, error_bad_lines=False)
                
                df, intsta_df = extract_internal_standards(df)
                
                intsta_df.dropna(how='all', inplace=True) # removes empty rows
                
                intsta_df['CLASS'] = intsta_df['SPECIES'].apply(lambda x: x.split(' ')[0])
                
                if ('FORMULA' and 'MASS' and 'ERROR') in df.columns:
                
                    return df, intsta_df
                
                else: 
                    
                   st.sidebar.error('This is not a valid LipidXplorer dataset!')
                
                   return None, None 
            
            else:
                
                st.sidebar.error('This is not a valid LipidXplorer dataset!')
                
                return None, None
            
            
            
            
            
# function to extract the internal standards dataframe from the uploaded dataset 
def extract_internal_standards(df):
    
        mol_lst = df['SPECIES'].values.tolist()
        
        intsta_lst = []
        
        for ele in mol_lst:
            
            ele_split = str(ele).split(' ')
            
            if (len(ele_split) > 1) and (str(ele_split[1]) == 'd7' or str(ele_split[1]) == 'd9'):
                
                index = mol_lst.index(ele)
                
                intsta_lst.append(index)
                
        intsta_df = df.iloc[intsta_lst, :]
        
        df.drop(intsta_lst,0,inplace=True) # removes the internal standard rows from the dataset
            
        return df, intsta_df 
    
    
    
    
    
 # building LipidXplorer dataframe (merging positive and negatve modes)
def build_merged_lipidxplorer_df(pos_df, pos_intsta_df, neg_df, neg_intsta_df):
        
        st.sidebar.info("""
                        
                        *Merging process*: In the positive-mode dataset, we only keep the lipids that belong to 
                        the following classes: Cer, HexCer, DAG, CE, LPC, LPC-O, PC, PC-O, SM and TG.
                        In the negative-mode dataset, we only keep the lipids that belong to 
                        the following classes: LPA, LPE, PE, PE-O, PA, PI, PG, PS. Then, we merge the two datasets.
                        
                        """)
        
        pos_class = ['Cer', 'HexCer', 'DAG', 'CE', 'LPC', 'LPC-O', 'PC', 'PC O-', 'PC-O', 'SM', 'TG']
        
        neg_class = ['LPA', 'LPE', 'PE', 'PE-O', 'PA', 'PI', 'PG', 'PS']

        pos_cols = ['SPECIES','FORMULA', 'MASS', 'ERROR', 'CLASS']
        
        neg_cols = ['SPECIES','FORMULA', 'MASS', 'ERROR', 'CLASS']
        
        # if both datasets are uploaded 
        
        if (pos_df is not None) and (neg_df is not None):
            
            # dealing with classes that are detected but not included in pos_class + neg_class
            
            all_class = pos_class + neg_class
            
            pos_df_unique_class = pos_df['CLASS'].unique()
            
            pos_class_rest = [x for x in pos_df_unique_class if (x not in all_class) and (str(x) != 'nan')]
            
            if len(pos_class_rest) > 0:
            
                st.sidebar.markdown('The following additional class(es) are detected in the positive mode. \
                                    Use the drop down menu to add additional classes to the list of the positive mode classes if necessary.')
            
                for lipid_class in pos_class_rest:
                    
                    st.sidebar.write(lipid_class)
                    
                pos_class_add = st.sidebar.multiselect('Add classes to the list of positive mode classes', pos_class_rest)
                
                pos_class = pos_class + pos_class_add
                
            neg_df_unique_class = neg_df['CLASS'].unique()
            
            neg_class_rest = [x for x in neg_df_unique_class if (x not in all_class) and (x not in pos_class) and (str(x) != 'nan')]
            
            if len(neg_class_rest) > 0:
            
                st.sidebar.markdown('The following additional class(es) are detected in the negative mode. \
                                    Use the drop down menu to add additional classes to the list of the negative mode classes.')
            
                for lipid_class in neg_class_rest:
                    
                    st.sidebar.write(lipid_class)
                    
                neg_class_add = st.sidebar.multiselect('Add classes to the list of the negative mode classes', neg_class_rest)
                
                neg_class = neg_class + neg_class_add
                
            # picking up the relevant classes from each mode and merging the two datasets 
        
            pos_df = pos_df.loc[pos_df['CLASS'].isin(pos_class)]
            
            pos_intsta_df.dropna(how='all', inplace=True) # removes empty rows
            
            pos_intsta_df['CLASS'] = pos_intsta_df['SPECIES'].apply(lambda x: x.split(' ')[0])
            
            pos_intsta_df = pos_intsta_df.loc[pos_intsta_df['CLASS'].isin(pos_class)]
            
            for column in pos_df:
                
                if 'PRECURSORINTENSITY' in column:
                    
                    pos_cols.append(column)
            
            pos_df = pos_df[pos_cols]
            
            pos_intsta_df = pos_intsta_df[pos_cols]
            
            neg_df = neg_df.loc[neg_df['CLASS'].isin(neg_class)]
            
            neg_intsta_df.dropna(how='all', inplace=True) # removes empty rows
            
            neg_intsta_df['CLASS'] = neg_intsta_df['SPECIES'].apply(lambda x: x.split(' ')[0])
            
            neg_intsta_df = neg_intsta_df.loc[neg_intsta_df['CLASS'].isin(neg_class)]
            
            for column in neg_df:
                
                if 'PRECURSORINTENSITY' in column:
                    
                    neg_cols.append(column)
            
            neg_df = neg_df[neg_cols]
            
            neg_intsta_df = neg_intsta_df[neg_cols]
            
            # if the two datasets have different number of columns (i.e samples)
            
            if len(pos_cols) != len(neg_cols):
                
                st.sidebar.error("""
                                 
                                 The two datasets must have equal number of relevant columns,
                                 otherwise merging is not possible.
                                 The datasets you uploaded do not satisfy the above condition!  
                                 
                                 """)
                                 
                return None, None
            
            else:
            
                unmatch_cols_dict = {'positive' : [], 'negative' : []} # includes the columns with un-matching names between the two datasets 
                
                counter = 0
            
                for i in range(len(pos_cols)):
                    
                    if pos_cols[i] != neg_cols[i]:
                        
                        counter = counter + 1
                        
                        unmatch_cols_dict['positive'].append(pos_cols[i])
                        
                        unmatch_cols_dict['negative'].append(neg_cols[i])
                        
                # if there are no columns with un-matching names 
                
                if counter == 0:
                    
                    pos_df['ERROR'] = pos_df['ERROR'].astype(str)
                    
                    neg_df['ERROR'] = neg_df['ERROR'].astype(str)
                    
                    pos_intsta_df['ERROR'] = pos_intsta_df['ERROR'].astype(str)
                    
                    neg_intsta_df['ERROR'] = neg_intsta_df['ERROR'].astype(str)
            
                    df = pd.concat([pos_df, neg_df])
                    
                    intsta_df = pd.concat([pos_intsta_df, neg_intsta_df])
            
                    return df, intsta_df
            
                # if there are columns with un-matching names 
                
                elif counter > 0:
                
                    st.sidebar.error("""
                                 
                                 The corresponding columns in both datasets must have the EXACT same name, otherwise merging is not possible.
                                 Look at the table below to find out which column names are not matching and correct the column names accordingly
                                 (use the "view full screen" option if necessary).
                                 
                                 """)
                                 
                    unmatch_cols_df = pd.DataFrame.from_dict(unmatch_cols_dict)
                
                    st.sidebar.write(unmatch_cols_df)
                
                    return None, None
        
        else:
            
            return None, None
        
        
        
        
        
# function that builds the side bar for LipidXplorer
def build_sidebar(df):
        
        st.sidebar.subheader("Define Experiment")
            
        n_cond = st.sidebar.number_input('Enter the number of conditions',
                                         min_value = 1, max_value= 20, value = 1, step = 1)
            
        cond_lst, rep_lst = build_cond_rep_lst(n_cond)
        
        # finding the total number of the samples 
        
        counter = 0 
        
        for column in df.columns:
            
            if 'PRECURSORINTENSITY' in column:
                
                counter = counter + 1
                
        # if inputs are not valid
                
        if counter != sum(rep_lst):
            
            st.sidebar.error('The inputs are incomplete and/or inaccurate!')
            
            return False, None, None, None, None, None 
            
        # if inputs are valid 
            
        else:
        
            st.sidebar.subheader('Group Samples')
        
            group_df = group_samples(df, cond_lst, rep_lst)
                
            st.sidebar.subheader("Apply Filters")
            
            st.sidebar.markdown("""
                        
                        A missing value does not necessarily mean zero abundance. It often means "no information available". 
                        If a datapoint has too many missing values, it can not be useful in the process of statistical inference. 
                        LipidCruncher allows the user to apply filters to the data to remove datapoints with too many missing values.
                        
                        For example, if the user sets the minimum required non-zero values at 3, AT LEAST ONE of the conditions 
                        has to have 3 or more non-zero values for the lipid species to pass though the filter. 
                        
                        """)
        
            filtered_conds, passing_abundance_grade = apply_additional_filter(cond_lst, rep_lst)
            
            st.sidebar.subheader("Confirm Inputs")
            
            name_df = update_sample_name(group_df, rep_lst)
        
            st.sidebar.markdown("""
                                
                                LipidCruncher uses the following protocol for naming the samples: s1, s2, ..., sN.
                                The table below shows the updated sample names (there is usually no difference between
                                the old names and the updated names unless there is a missing sample or the samples were 
                                not originally grouped together properly):
                                
                                """)
        
            st.sidebar.write(name_df)
        
            st.sidebar.write("Now, confirm your inputs:")
        
            st.sidebar.write("There are a total of "+str(sum(rep_lst))+" samples.")
            
            for cond in cond_lst:
            
                build_rep_cond_pair(cond, cond_lst, rep_lst) # function defined below 
            
            confirm_data = st.sidebar.checkbox("Confirm the inputs by checking this box")
            
        return confirm_data, cond_lst, rep_lst, group_df, filtered_conds, passing_abundance_grade
    
    
    
    
    
# function that builds cond_lst and rep_lst objects (conditions list and number of replicates list)
def build_cond_rep_lst(n_cond): 
    
        '''
        For example: if we have two conditions, WT and KO, and each have two corresponding replicates, 
        the cond_lst = [WT, KO] and the rep_lst = [2, 2]
        
        '''
            
        cond_lst = []
            
        rep_lst = []
            
        for i in range(1, n_cond+1):
                
            cond_lst.append(cond_text_input(i)) # function defined below 
                
            rep_lst.append(rep_number_input(i)) # function defined below
                
        return cond_lst, rep_lst
    
    
    
    
    
# function that creates a text input box to enter the label for each condition
def cond_text_input(i):  
        
        cond = st.sidebar.text_input('Create a label for condition #'+str(i)+' (e.g. WT or KO)')
                
        return cond
    
    
    
    
            
    # function that creates a text input box to enter the number of replicates corresponding each condition
def rep_number_input(i):  
        
        rep = st.sidebar.number_input('Enter the number of replicates for condition #'+str(i), 
                                              min_value = 1, max_value = 1000, value = 1, step = 1)
                
        return rep
    
    
    
    
    
# function to group together samples that belong to the same condition
def group_samples(df, cond_lst, rep_lst): 
    
        # give user the necessary info
        
        rep_lst_agg = build_rep_lst_agg(rep_lst)
                             
        for cond in cond_lst:
                
            index = cond_lst.index(cond)
                
            if index == 0:
                             
                st.sidebar.write('- Samples indexed from 0 to ' + str(rep_lst_agg[0]-1) + ' must belong to ' + cond)
                    
            else: 
                    
                st.sidebar.write('- Samples indexed from ' + str(rep_lst_agg[index-1]) + ' to ' + str(rep_lst_agg[index]-1) 
                                 + ' must belong to ' + cond)
                
        group_df = build_group_df(df, cond_lst, rep_lst)
        
        st.sidebar.write(group_df)
        
        st.sidebar.write('Are your samples properly grouped together?')
        
        ans = st.sidebar.radio('', ['Yes', 'No'])
        
        if ans == 'Yes':
            
            st.sidebar.write('Go to the next section!')
            
        else:
            
            ordered_col_lst = []
            
            col_lst = group_df['sample name'].tolist()
            
            for cond in cond_lst:
                
                temp = st.sidebar.multiselect('Pick the samples that belong to condition ' + cond, col_lst)
                
                ordered_col_lst += temp
                
                col_lst = [ele for ele in col_lst if ele not in temp]
            
            if len(ordered_col_lst) == sum(rep_lst):
            
                group_df['sample name'] = ordered_col_lst
            
            st.sidebar.write('Check the updated table below to make sure the samples are properly grouped together:') 
            
            st.sidebar.write(group_df)
            
        return group_df
    
    
    
    
    
# function to build the aggregate list of rep_lst elements
@st.cache
def build_rep_lst_agg(rep_lst):  
            
        '''
        For example, if rep_lst = [2, 4, 5], the rep_lst_agg = [2, 6, 11]
        
        '''
            
        rep_lst_agg = [] # the aggregate list of replicates 
            
        for i in range(len(rep_lst)):
                
            rep_lst_agg.append(sum(rep_lst[0:i+1]))
                
        return rep_lst_agg
    
    
    
    
    
# function to build a df with two columns: sample names and corresponding conditions
def build_group_df(df, cond_lst, rep_lst):  
        
        extensive_cond_lst = [] # list includes the cond correponding to each rep - length equal to the total number of reps 
        
        for cond in cond_lst:
            
            index = cond_lst.index(cond)
            
            extensive_cond_lst += [cond for i in range(rep_lst[index])]
            
        group_dict = {'sample name' : [], 'condition' : []}
         
        for col in df.columns:
             
             if 'PRECURSORINTENSITY:' in col:
                 
                 a = col.split(':')
                 
                 b = a[1].split('.')
     
                 group_dict['sample name'].append(b[0]) 
                 
             elif 'PRECURSORINTENSITY[' in col:
                
                 group_dict['sample name'].append(col[19: -1])
        
        group_dict['condition'] = extensive_cond_lst
            
        group_df = pd.DataFrame.from_dict(group_dict) # df including sample_name and corresponding condition 
        
        return group_df 
    
    
    
    
    
# function to remove datapoints with too many missing values 
def apply_additional_filter(cond_lst, rep_lst):
            
        filtered_conds = st.sidebar.multiselect('Specify which conditions to apply the filter to (remove/add conditions)', cond_lst, cond_lst)
            
        filtered_reps = [] # includes the number of reps corresponding to the selected_conds list 
            
        for cond in filtered_conds:
                
            index = cond_lst.index(cond)
                
            filtered_reps.append(rep_lst[index])
            
        passing_abundance_grade = st.sidebar.number_input("Enter the minimum required number of non-zero values"
                                                          , min_value = 0, max_value = np.min(filtered_reps), value = 0, step = 1)
        
        return filtered_conds, passing_abundance_grade
    
    
    
    
    
# function that updates the sample names, so, they are consistent with the following format: s1, s2, s3, ..., sN. 
@st.cache
def update_sample_name(group_df, replicate_lst):  
    
        """
        For example, if you have 4 samples with the following names: WT1, WT2, BQC1, BQC2 or s1, s2, s3, s5, then the updated names are s1, s2, s3 and s4.  
        
        """
        
        name_dict = {"old name" : [] , "updated name" : [], "condition": []}
        
        name_dict['old name'] = group_df['sample name']
        
        name_dict['updated name'] = ['s'+str(i+1) for i in range(sum(replicate_lst))]
        
        name_dict['condition'] = group_df['condition']
                
        name_df = pd.DataFrame.from_dict(name_dict)
        
        return name_df
    
    
    
    
    
# pairs the replicates with their corresponding condition
def build_rep_cond_pair(cond, cond_lst, rep_lst):  
    
        '''
        For example, if cond_lst = [WT, KO] and rep_lst = [2, 2], the function returns:
            s1-s2 correspond to WT
            s3-s4 correspond to KO
        '''
                
        rep_lst_agg = build_rep_lst_agg(rep_lst)
            
        sample_lst = build_sample_lst(cond, cond_lst, rep_lst_agg)
            
        if len(sample_lst) == 1:
                
            return st.sidebar.write("- "+sample_lst[0]+" corresponds to "+cond)
            
        else:
        
            return st.sidebar.write("- "+sample_lst[0]+"-"+sample_lst[-1]+" correspond to "+cond)
        
        
        
        
        
# function to build the list of samples (not the number of replicates) for each condition
@st.cache
def build_sample_lst(cond, cond_lst, rep_lst_agg):  
        
            '''
            For example, if cond_lst = [WT, KO, HA] and rep_lst = [4, 4, 2], then rep_lst_agg = [4, 8, 10].
            In this case, for cond = WT, the sample_lst = [s1, s2, s3, s4], for cond = KO, the sample_lst = 
            [s5, s6, s7, s8] and for cond = HA, the sample_lst = [s9, s10].
            
            '''
            
            if cond_lst.index(cond) == 0:
                    
                sample_lst = ['s'+str(i+1) for i in range(rep_lst_agg[cond_lst.index(cond)])]
                    
            else:
                    
                sample_lst = ['s'+str(i+1) for i in range(rep_lst_agg[cond_lst.index(cond)-1], \
                                                                          rep_lst_agg[cond_lst.index(cond)])]
                        
            return sample_lst
        
        
        
        
        
# function to clean/apply filters to the data
def apply_filter(df, intsta_df, group_df, rep_lst, cond_lst, filtered_conds, passing_abundance_grade):
    
        # initial cleaning 
        
        counter = 0
        
        for column in df.columns: # enforces the LipidSearch sample naming convention on LipidXplorer samples
            
            if 'PRECURSORINTENSITY' in column:
                
                df.rename(columns={column: 'PRECURSORINTENSITY[s'+str(counter+1)+']'}, inplace=True)
                
                intsta_df.rename(columns={column: 'PRECURSORINTENSITY[s'+str(counter+1)+']'}, inplace=True)
                
                counter = counter + 1
                
        total_reps = sum(rep_lst) # total number of all replicates
        
        X = df[['SPECIES', 'FORMULA', 'MASS', 'ERROR', 'CLASS']+
                ['PRECURSORINTENSITY[s' + str(i+1) + ']' for i in range(total_reps)]]
                
        intsta_df = intsta_df[['SPECIES', 'FORMULA', 'MASS', 'ERROR', 'CLASS']+
                ['PRECURSORINTENSITY[s' + str(i+1) + ']' for i in range(total_reps)]]
                
        # cleaning X
        
        X = X[X.SPECIES != '###'] # removes the rows that start with ### (comments)
        
        X.dropna(how='all', inplace=True) # removes empty rows  
        
        X.drop_duplicates(keep = 'first', inplace = True) # removes duplicate datapoints 
        
        auc = ['PRECURSORINTENSITY[s'+str(i+1)+']' for i in range(sum(rep_lst))]
        
        X[auc]=X[auc].mask(X[auc]<0).fillna(0)
        
        X.reset_index(inplace=True)
        
        #first filter
        
        X['ERROR'] = pd.to_numeric(X['ERROR'], downcast="float")
        
        X = X.loc[abs(X['ERROR']) < 5] # removes the datapoint if abs(error) > 5
        
        #second filter: minimum abundance grade 
        
        X = build_abundance_grade_filter(X, rep_lst, cond_lst, filtered_conds, passing_abundance_grade)
        
        X.reset_index(inplace=True)
        
        X.drop(['level_0'], axis=1, inplace=True) # drops an irrelevant column
        
        X.rename(columns={"index": "old_index"}, inplace=True)
        
        X.set_index('SPECIES', inplace = True)
        
        # cleaning internal standards df
        
        intsta_df.dropna(how='all', inplace=True) # removes empty rows
        
        intsta_df.reset_index(inplace=True) 
        
        intsta_df.rename(columns={"index": "old_index"}, inplace=True)
        
        intsta_df.set_index('SPECIES', inplace = True)
        
        st.write('View the cleaned data in conventionl format:')
                
        st.write(X)
        
        csv_download = convert_df(X)
                
        st.download_button(
            
                    label="Download Data",
                    
                    data=csv_download,
                    
                    file_name='cleaned_data.csv',
                    
                    mime='text/csv')
        
        st.write('-----------------------------------------------------------------------------')
        
        st.write('View the cleaned data in log-transformed format:')
                
        log_X = log_transform_df(X, rep_lst)
            
        st.write(log_X)
        
        csv_download = convert_df(log_X)
                
        st.download_button( 
            
                    label="Download Data",
                    
                    data=csv_download,
                    
                    file_name='log_transformed_cleaned_data.csv',
                    
                    mime='text/csv')
        
        st.write('-----------------------------------------------------------------------------')
        
        st.write('View the internal standards data in conventional format: ')
                
        st.write(intsta_df)
        
        csv_download = convert_df(intsta_df)
                
        st.download_button(
            
                    label="Download Data",
                    
                    data=csv_download,
                    
                    file_name='internal_standards.csv',
                    
                    mime='text/csv')
        
        st.write('-----------------------------------------------------------------------------')
        
        st.write('View the internal standards data in log-transformed format:')
                
        log_intsta_df = log_transform_df(intsta_df, rep_lst)
                
        st.write(log_intsta_df)
        
        csv_download = convert_df(log_intsta_df)
                
        st.download_button(
            
                    label="Download Data",
                    
                    data=csv_download,
                    
                    file_name='log_transformed_internal_standards.csv',
                    
                    mime='text/csv')
                
        return X, intsta_df
    
    
    
    
    
# function that removes the lipid species that have an abundance grade lower than the minimum abundance grade
def build_abundance_grade_filter(X, rep_lst, cond_lst, filtered_conds, passing_abundance_grade):
        
        rep_lst_agg = build_rep_lst_agg(rep_lst)  
        
        for cond in filtered_conds:
            
            sample_lst = build_sample_lst(cond, cond_lst, rep_lst_agg)
            
            # adding the 'abundance grade' column for each filtered condition 
            
            X['abundance_grade_'+cond] = X[['PRECURSORINTENSITY['+sample+']' for sample in sample_lst]]\
                .apply(lambda x: abundance_level_calculator(x, passing_abundance_grade),axis=1)
            
        X['abundance_grade'] = X[['abundance_grade_'+cond for cond in filtered_conds]].apply(lambda x: sum(x), axis=1) # adds the 'total abundance grade' column 
        
        X = X.loc[X['abundance_grade'] > 0 ] # keeps only lipid species with non-zero 'total abundance grade'  
        
        X.drop(['abundance_grade_'+cond for cond in filtered_conds]+['abundance_grade'], axis=1, inplace=True) # drops the abundance grade columns
        
        return X
    
    
    
    
    
@st.cache
def convert_df(dataframe):

        return dataframe.to_csv().encode('utf-8')
    
    
    
    
    
# function to log transform the abundance columns of the df
@st.cache
def log_transform_df(X, rep_lst):  
    
        temp = X.copy()
        
        auc = ['PRECURSORINTENSITY[s'+str(i+1)+']' for i in range(sum(rep_lst))]
        
        # filling zero values with 1's to avoid infinity
        
        temp[auc]=temp[auc].mask(temp[auc]<=0).fillna(1)
        
        temp[auc] = temp[auc].apply(lambda x: np.log10(x), axis=0)
        
        return temp
    
    
    
    
    
#function that computes the total abundance grade for each lipid species
@st.cache
def abundance_level_calculator(intensities, passing_abundance_grade): 
        
        '''
        For example: if the passing abundance grade is 3, and condition A has a at least 3 non-zero values, 
        the function returns 1, otherwise it returns 0.
        '''
        
        counter = 0
                    
        for intensity in intensities:
                        
            if float(intensity) > 0:
                            
                counter = counter+1
                            
        if counter >= passing_abundance_grade:
                        
            return 1  
                    
        else:
                        
            return 0
        
        
        
        
        
# function to plot the CoV of lipid species
def plot_cov(X, cond_lst, rep_lst):  
        
        auc = ['PRECURSORINTENSITY[s'+str(i+1)+']' for i in range(sum(rep_lst))]
        
        cond = st.radio('Which of the following conditions corresponds to BQC samples? ', cond_lst+['None of the above'], len(cond_lst))
            
        if cond != 'None of the above':
                
                index = cond_lst.index(cond)
                
                if rep_lst[index] == 1:
                    
                    st.error('The selected condition must have at least two corresponding replicates!')
                    
                else:
                    
                    X[auc]=X[auc].mask(X[auc]<=0).fillna(1) # turning 0's and negative numbers  to 1's so it is possible to log transform
            
                    rep_lst_agg = build_rep_lst_agg(rep_lst)
                
                    sample_lst = build_sample_lst(cond, cond_lst, rep_lst_agg)
                    
                    X, X_plot, X_cov_df = cov_hover(X, sample_lst)
                    
                    show_cov = st.checkbox("View CoV plot")
                    
                    if show_cov:
                
                        st.bokeh_chart(X_plot)
                    
                        csv_download = convert_df(X_cov_df)
                            
                        st.download_button(
                        
                                label="Download Data",
                                
                                data=csv_download,
                                
                                file_name='cov.csv',
                                
                                mime='text/csv')
                    
                    filter_ans = st.radio('Would you like to filter the data by removing datapoints with high CoV?', ['Yes', 'No'], 1)
                    
                    st.warning("Your choice here would affect the rest of the module.")
                    
                    if filter_ans == 'Yes':
                        
                        thresh = st.number_input('Enter the maximum acceptable CoV in %', min_value = 10, max_value = 100, value = 30, step = 1)
                        
                        X = X.loc[X['cov'] <= thresh] # removes datapoints with CoV > thresh
                    
                        X[auc]=X[auc].mask(X[auc]==1).fillna(0) # turning 1's back to 0's
                    
                        X.reset_index(inplace=True)
                    
                        X.drop(['mean', 'cov'], axis=1, inplace=True)
                        
                        X.set_index('SPECIES', inplace = True)
                    
                        st.write(X)
            
                        csv_download = convert_df(X)
                    
                        st.download_button(
                        
                            label="Download Data",
                        
                            data=csv_download,
                        
                            file_name='filtered_data.csv',
                        
                            mime='text/csv')
                    
                    else:
                        
                        X[auc]=X[auc].mask(X[auc]==1).fillna(0) # turning 1's back to 0's
                
                        X.reset_index(inplace=True)
                
                        X.drop(['mean', 'cov'], axis=1, inplace=True)
                        
        return X
    
    
    
    
    
# one CoV plot with hover tool
def cov_hover(X, sample_lst):  
            
        X['cov'] = X[['PRECURSORINTENSITY['+sample+']' for sample in sample_lst]].apply(lambda x: cov_calculator(x), axis=1)
    
        X['mean'] = X[['PRECURSORINTENSITY['+sample+']' for sample in sample_lst]].apply(lambda x: mean_calculator(x), axis=1)
            
        plot = figure(title='CoV - All lipid Species', x_axis_label='Mean of log10(AUC) Across All BQC Samples', y_axis_label='CoV(%)')
            
        x = X['mean'].values.tolist()
            
        y = X['cov'].values.tolist()
        
        species = X.index.values.tolist()
            
        cov_df = pd.DataFrame({"Mean_AUC": x, "CoV": y, 'Species': species})
        
        src = ColumnDataSource(cov_df)
        
        plot.scatter(x="Mean_AUC", y="CoV", name='cov', source=src)
    
        plot.line(x=[i for i in range(12)], y=30, color='red')
        
        hover = HoverTool(tooltips = [('Mean_AUC', '@Mean_AUC'), ('CoV', "@CoV"), ('Species', "@Species")], names=['cov'])
        
        plot.add_tools(hover)
        
        plot.title.text_font_size = "15pt"
        
        plot.xaxis.axis_label_text_font_size = "15pt"
            
        plot.yaxis.axis_label_text_font_size = "15pt"
            
        plot.xaxis.major_label_text_font_size = "15pt"
            
        plot.yaxis.major_label_text_font_size = "15pt"
    
        return X, plot, cov_df 
    
    
    
    
    
# calculate CoV
@st.cache
def cov_calculator(numbers): 
    
    non_zero_lst = [number for number in numbers if (number>1)]
    
    if len(non_zero_lst) > 0:

        cov = np.std(non_zero_lst)/np.mean(non_zero_lst)*100
        
    else:
        
        cov = None 
        
    return cov
    
    
    
    
    
 # calculate mean
@st.cache
def mean_calculator(numbers): 
    
    non_zero_lst = [number for number in numbers if (number>1)]
    
    if len(non_zero_lst) > 0:

        mean = np.mean(non_zero_lst)
    
        mean = np.log10(mean)
        
    else:
        
        mean = None
        
    return mean




# function to plot histograms 
def plot_hist(X, rep_lst, cond_lst):
        
        show_hist = st.checkbox('View the distributions of the Area Under the Curve (AUC)')
        
        if show_hist:
                        
            temp = X.copy()
                
            full_sample_lst = ['s'+str(i+1) for i in range(sum(rep_lst))]
            
            number_rep = sum(rep_lst)
            
            filter_ans = st.radio(' Would you like to filter by lipid class?', ['Yes', 'No'], 1)
            
            if filter_ans == 'Yes':
                
                lipid_class = st.selectbox('Select a lipid class', X['Class'].unique())
                
                temp = temp[temp['Class'] == lipid_class]
            
            temp = temp[['PRECURSORINTENSITY[s'+str(i+1)+']' for i in range(sum(rep_lst))]] # picks the 'MainArea[s1]', ..., 'MainArea[sN]' columns only
            
            temp = temp.mask(temp<=0).fillna(1)
            
            if number_rep < 50:
            
                arrange_hist(number_rep, temp, full_sample_lst)
                
            else:
                
                st.warning("LipidCruncher has limited available resources. You have more than 50 samples, \
                           the best practice is to plot no more than 50 plots at a time.")
                           
                st.write('Use the double-slider below to set a range of samples you want to investigate. \
                         Make sure to not to go beyond the 50 samples limit.')
                         
                srange = st.slider("", value=[0, number_rep])
                
                diff = srange[1] - srange[0]
                
                if diff <= 50: 
                    
                    show_range_hist = st.button("Show Plots!")
                    
                    if show_range_hist:
                    
                        arrange_hist(diff, temp[['PRECURSORINTENSITY[s'+str(i+1)+']' for i in range(srange[0], srange[1])]], full_sample_lst[srange[0]: srange[1]])
                    
                else: 
                    
                    st.error("You have exceeded the 50 samples limit!")
                        
            return
        
        
        
        
        
# function to arrange histograms in in subplot style, three histograms per row
def arrange_hist(number_rep, temp, full_sample_lst):
        
        for i in range(int(round(number_rep/3, 1))):
                            
                col1, col2, col3 = st.columns(3)
                            
                fig = prep_hist(temp, full_sample_lst, 3*i)
                            
                col1.pyplot(fig)
                            
                fig = prep_hist(temp, full_sample_lst, 3*i+1)
                            
                col2.pyplot(fig)
                            
                fig = prep_hist(temp, full_sample_lst, 3*i+2)
                            
                col3.pyplot(fig)
                            
        if number_rep % 3 == 2:
                            
                col1, col2, col3 = st.columns(3)
                            
                fig = prep_hist(temp, full_sample_lst, -2)
                            
                col1.pyplot(fig)
                            
                fig = prep_hist(temp, full_sample_lst, -1)
                            
                col2.pyplot(fig)
                            
        elif number_rep %3 == 1:
                            
                col1, col2, col3 = st.columns(3)
                            
                fig = prep_hist(temp, full_sample_lst, -1)
                            
                col1.pyplot(fig)
        
        return
    
    
    
    
    
# function that prepares the histogram plot 
def prep_hist(temp, full_sample_lst, index):
                        
        plt.rcParams['font.size'] = '35'
                    
        plt.rcParams['axes.linewidth'] = 3
                
        fig, ax = plt.subplots(figsize=(10, 10))
        
        lst = np.log10(temp['PRECURSORINTENSITY[' + full_sample_lst[index] + ']'].values.tolist())
            
        ax.hist(lst, bins = 75, range=(0, 12))
                    
        ax.set_title('Histogram of AUC - '+ full_sample_lst[index], fontsize=50)
                            
        ax.set_xlabel('log10(AUC)', fontsize=50)
                            
        ax.set_ylabel('Count', fontsize=50)
                        
        ax.set_xlim([-0.5, 12])
                        
        return fig
    
    
    
    
    
# pairwise correlation plots
def plot_corr(X, cond_lst, rep_lst): 
    
        show_corr = st.checkbox("Run pairwise correlation tests")
        
        if show_corr:
    
            cond = st.selectbox('Select a condition', cond_lst)
        
            rep_type = st.selectbox('Select the type of your replicates',
                                    ['biological replicates', 'Technical replicates'])
        
            if rep_type == 'biological replicates':
                
                # set the min and the center of the color bar 
            
                v_min = 0.5
            
                thresh = 0.8
            
            elif rep_type == 'Technical replicates':
            
                v_min = 0.75
            
                thresh = 0.9
        
            index = cond_lst.index(cond)
            
            # if there's more than one replicate
        
            if rep_lst[index] > 1:
        
                temp = X[['PRECURSORINTENSITY[s' + str(i+1) + ']' for i in range(sum(rep_lst))]].copy()
        
                temp = temp.apply(lambda x: np.log10(x))
            
                # re-naming the columns from MainArea[s1] to s1 etc
        
                counter = 0
        
                for column in temp.columns:
            
                    temp.rename(columns={column: 's' + str(counter+1)}, inplace=True)
            
                    counter = counter + 1
            
                # selecting the columns that correspond to the selected cond    
            
                rep_lst_agg = build_rep_lst_agg(rep_lst)
                
                sample_lst = build_sample_lst(cond, cond_lst, rep_lst_agg)
                
                temp = temp[[sample for sample in sample_lst]]
        
                fig = plt.figure(figsize=(20, 16))
            
                mask = np.triu(np.ones_like(temp.corr(), dtype=np.bool))
    
                heatmap = sns.heatmap(temp.corr(), mask=mask, vmin=v_min, vmax=1, center = thresh, annot=False, cmap='RdBu', square=False, cbar=True)
            
                heatmap.set_title('Triangle Correlation Heatmap - ' + cond, fontdict={'fontsize':30});
            
                st.pyplot(fig)
                
                st.write('---------------------------------------------------------------------------------------------------')
                
                st.write('Find the exact correlation coefficients in the table below:')
        
                st.write(temp.corr())
            
                csv_download = convert_df(temp.corr())
                            
                st.download_button(
                    
                                label="Download Data",
                                
                                data=csv_download,
                                
                                file_name='Correlation_Matrix_'+str(cond)+'.csv',
                                
                                mime='text/csv')
            
            else:
            
                st.error('The selected condition must have at least two corresponding replicates!')
            
        return
    
    
    
    
    
# plotting PCA
def plot_pca(dataframe, rep_lst, cond_lst):
    
            temp = dataframe.copy()
            
            replicate_lst = rep_lst.copy()
            
            condition_lst = cond_lst.copy()
    
            temp = temp[['PRECURSORINTENSITY[s'+str(i+1)+']' for i in range(sum(replicate_lst))] + ['CLASS']]
            
            full_sample_lst = ['s'+str(i+1) for i in range(sum(replicate_lst))]
            
            show_pca = st.checkbox("Run Principal Component Analysis (PCA)")
            
            if show_pca:
                
                remove_ans = st.radio("Would you like to remove any samples from the analysis?", ['Yes', 'No'], 1)
                
                if remove_ans == 'Yes':
                    
                    r_sample = st.multiselect('Pick the sample(s) that you want to remove from the analysis', full_sample_lst)
                    
                    if (len(full_sample_lst) - len(r_sample)) >= 3:
                    
                        replicate_lst, condition_lst, pca_name_cond_df, temp = update_df_cond_rep_sample_lst(replicate_lst, condition_lst, r_sample, temp)
                        
                        st.write('The table below shows the updated samples names:')
                        
                        st.warning('Samples removed in the PCA section, only affect the PCA analysis not the rest of the module.')
                        
                        st.write(pca_name_cond_df)
                        
                    else:
                        
                        st.error('At least three samples are required for a meanigful analysis!')
                        
                        return None
                        
                filter_ans = st.radio("Would you like to filter by lipid class?", ['Yes', 'No'], 1)
                
                if filter_ans == 'Yes':
                    
                    lipid_class = st.selectbox('Pick a lipid class', temp['CLASS'].value_counts().index.tolist())
                    
                    temp = temp[temp['CLASS'] == lipid_class]
            
                temp.drop(['CLASS'], axis=1, inplace=True)
               
                scores, explained_variance = run_pca(temp)
                
                PC_lst = ['PC'+str(i+1)+' ('+str("{:.0f}".format(explained_variance[i]*100))+'%)' for i in range(3)] # e.g. [PC1 (80%), PC2 (15%), PC3(5%)]
                
                sel_PCx = st.selectbox('Pick the PC shown on the x-axis', PC_lst, 0)
                
                sel_PCy = st.selectbox('Pick the PC shown on the y-axixs', PC_lst, 1)
                
                x_index = int(sel_PCx.split(' ')[0][2]) # 1, 2 or 3
                
                y_index = int(sel_PCy.split(' ')[0][2]) # 1, 2 or 3
                
                PCx = scores[:, x_index - 1].tolist()
                
                PCy = scores[:, y_index - 1].tolist()
            
                rep_lst_agg = build_rep_lst_agg(replicate_lst)
            
                color_lst =[ 'red', 'blue', 'green', 'magenta', 'cyan', 'orange', 'black', 'pink', 'brown', 'yellow', 'purple', \
                            'gray', 'olive', 'chocolate', 'silver', 'darkred', 'khaki', 'skyblue', 'navy', 'orchid']
            
            
                plot = figure(title='PCA', x_axis_label=sel_PCx, y_axis_label=sel_PCy)
            
                x = []
            
                y = []
            
                legend = []
            
                color = []
            
                for cond in condition_lst:
                
                    if condition_lst.index(cond) == 0:
                    
                        x = x + PCx[0: rep_lst_agg[condition_lst.index(cond)]]
                    
                        y = y + PCy[0: rep_lst_agg[condition_lst.index(cond)]]
                    
                        legend = legend + [cond for i in range(rep_lst_agg[condition_lst.index(cond)])]
                    
                        color = color + [color_lst[condition_lst.index(cond)] for i in range(rep_lst_agg[condition_lst.index(cond)])]
                        
                    else:
                    
                        x = x + PCx[rep_lst_agg[condition_lst.index(cond)-1]: rep_lst_agg[condition_lst.index(cond)]]
                    
                        y = y + PCy[rep_lst_agg[condition_lst.index(cond)-1]: rep_lst_agg[condition_lst.index(cond)]]
                    
                        legend = legend + [cond for i in range(rep_lst_agg[condition_lst.index(cond)-1], rep_lst_agg[condition_lst.index(cond)])]
                    
                        color = color + [color_lst[condition_lst.index(cond)] for i in range(rep_lst_agg[condition_lst.index(cond)-1], rep_lst_agg[condition_lst.index(cond)])]
                
                full_sample_lst = ['s'+str(i+1) for i in range(sum(replicate_lst))]
                
                pca_df = pd.DataFrame({"PC1": x, "PC2": y, "sample": full_sample_lst, "legend": legend, "color": color})
            
                src = ColumnDataSource(pca_df)
            
                plot.scatter(x="PC1", y="PC2", legend_group='legend', color='color', source=src)
            
                hover = HoverTool(tooltips = [('PC1', '@PC1'), ('PC2', '@PC2'), ('sample', '@sample')])
            
                plot.add_tools(hover)
            
                plot.title.text_font_size = "15pt"
            
                plot.xaxis.axis_label_text_font_size = "15pt"
            
                plot.yaxis.axis_label_text_font_size = "15pt"
            
                plot.xaxis.major_label_text_font_size = "15pt"
            
                plot.yaxis.major_label_text_font_size = "15pt"
            
                if show_pca:
                
                    st.bokeh_chart(plot)
                    
                    csv_download = convert_df(pca_df[['PC1', 'PC2', 'sample', 'legend']])
                            
                    st.download_button(
                        
                                label="Download Data",
                                
                                data=csv_download,
                                
                                file_name='PCA.csv',
                                
                                mime='text/csv')
                            
            return
        
        
        
        
        
# PCA math 
@st.cache
def run_pca(dataframe): # PCA 
    
            # turning 0's and negative numbers to 1's for log-transformation 
            
            dataframe = dataframe.mask(dataframe<=0).fillna(1)
            
            dataframe = dataframe.apply(lambda x: np.log10(x))
            
            # transposing the dataset 
            
            Y = dataframe.T
                    
            Y = scale(Y)
            
            pca = decomposition.PCA(n_components=3)
            
            pca.fit(Y)
            
            scores = pca.transform(Y)
            
            explained_variance = pca.explained_variance_ratio_
            
            return scores, explained_variance
        
        
        
        
        
# updates the datafraame and rep, cond, full_sample lists after removing a sample
def update_df_cond_rep_sample_lst(replicate_lst, condition_lst, r_sample, dataframe): 
        
        '''
        For example, if we have two conditions, A and B, with 3 replicates each, cond lst = [A, B]
        and repl_lst = [3, 3]. The dataframe will have 6 columns.
        If we remove let's say s5, the following function will update the lists: 
        cond_lst = [A, B], rep_lst = [3, 2], and the dataframe will have 5 columns. 
        
        '''
        
        # dropping bad samples and updating the column names 
        
        dataframe.drop(['PRECURSORINTENSITY[' + sample + ']' for sample in r_sample], axis=1, inplace=True) # updating dataframe
            
        # updating rep_lst and cond_lst
        
        r_sample = [int(item[1:]) for item in r_sample] # separates string and number (e.g. extracts 11 from s11)
        
        rep_lst_agg = build_rep_lst_agg(replicate_lst)
        
        used_sample = [] # useful list to avoid overwriting a sample multiple times 
        
        r_sample_cond_index = [] # lst including the index of the condition that corresponds to the removed sample 
        
        for sample_number in r_sample:
            
            for agg_rep in rep_lst_agg:
                
                if (sample_number <= agg_rep) & (sample_number not in used_sample):
                    
                    used_sample.append(sample_number)
                    
                    r_sample_cond_index.append(rep_lst_agg.index(agg_rep))
                    
        for index in r_sample_cond_index:
            
            replicate_lst[index] = replicate_lst[index] - 1 # updating rep_lst
            
        condition_lst = [cond for cond in condition_lst if replicate_lst[condition_lst.index(cond)] != 0] # removing cond with zero corresponding reps 
            
        replicate_lst = [rep for rep in replicate_lst if rep !=0] # remoing 0's from rep_lst
        
        group_df = build_group_df(dataframe, condition_lst, replicate_lst)
        
        name_df = update_sample_name(group_df, replicate_lst)
        
        # updating column names 
        
        total_reps = sum(replicate_lst) # total number of all replicates
        
        for (sample_1, sample_2) in zip(name_df['old name'], name_df['updated name']):
            
            dataframe.rename(columns={'PRECURSORINTENSITY[' + sample_1 + ']' : 'Area[' + sample_2 + ']'}, inplace = True)
            
        for i in range(total_reps):
            
            dataframe.rename(columns={'Area[s' + str(i+1) + ']' : 'PRECURSORINTENSITY[s' + str(i+1) + ']'}, inplace = True)
        
        return replicate_lst, condition_lst, name_df, dataframe
    
    
    
    
    
# function to to impute missing values and remove bad samples 
def remove_bad_sample(X, rep_lst, cond_lst): 
            
            # removing bad samples
            
            full_sample_lst = ['s'+str(i+1) for i in range(sum(rep_lst))]
                
            remove_ans = st.radio("Would you like to remove any samples from the analysis? ", ['Yes', 'No'], 1)
                    
            if remove_ans == 'Yes':
                        
                r_sample = st.multiselect('Pick the sample(s) that you want to remove from the analysis ', full_sample_lst)
                
                st.warning("Your choice here would be applied to the rest of the module.")
                        
                rep_lst, cond_lst, name_cond_df, X = update_df_cond_rep_sample_lst(rep_lst, cond_lst, r_sample, X)
                
                st.write('The table below shows the updated samples names:')
                
                st.write(name_cond_df)
            
            view_clean_data = st.checkbox('View the final cleaned data')
            
            if view_clean_data:
                
                st.write('View the cleaned data in the conventional format:')
                
                st.write(X)
                
                csv_download = convert_df(X)
                            
                st.download_button(
                    
                                label="Download Data",
                                
                                data=csv_download,
                                
                                file_name='final_cleaned_data.csv',
                                
                                mime='text/csv')
            
                st.write('------------------------------------------------------------------------------------------------')
            
                st.write('View the cleaned data in the log-transformed format:')
                
                log_X = log_transform_df(X, rep_lst)
            
                st.write(log_X)
                
                csv_download = convert_df(log_X)
                            
                st.download_button(
                    
                                label="Download Data",
                                
                                data=csv_download,
                                
                                file_name='final_log_transformed_cleaned_data.csv',
                                
                                mime='text/csv')
            
            return X, rep_lst, cond_lst
        
        
        
        
        
# function to normalize AUC data
def normalize_auc(X, intsta_df, rep_lst):
        
        # extracting the list of classes for which no internal standard is found 
        
        temp = X.copy()
        
        temp.reset_index(inplace=True)
        
        intsta_df.reset_index(inplace=True)
        
        all_class_lst = temp['CLASS'].unique().tolist()
        
        intsta_class_lst = intsta_df['CLASS'].unique().tolist()
        
        no_intsta_class_lst = [x for x in all_class_lst if x not in intsta_class_lst]
        
        st.write("An internal standard is found for the follwing classes:")
        
        for lipid_class in intsta_class_lst:
            
            st.write(lipid_class)
        
        if len(no_intsta_class_lst) > 0:
        
            st.write("However, no internal standards is found for the following classes:")
        
            for lipid_class in no_intsta_class_lst:
            
                st.write(lipid_class)
                
            st.write("For each lipid class that does not have a matching internal standard\
                     , pick one internal standard from the existing list of internal standards \
                     that most closely matches the lipid class:")
            
        # letting the user pick IS for classes with no IS
        
        add_intsta_class_lst = pick_intsta(no_intsta_class_lst, intsta_class_lst)
        
        # the amount of each IS in micro mole 
        
        st.write("Now, enter the concentration of each internal standard:")
        
        intsta_mic_mol_lst = build_intsta_mic_mol_lst(intsta_class_lst)
            
        for i in range(sum(rep_lst)):
                
            temp['normalized_AUC[s'+str(i+1)+']'] = 0
                
        # normalizing the data 
            
        for lipid_class in all_class_lst:
                
            if lipid_class in intsta_class_lst:
                    
                class_index = intsta_class_lst.index(lipid_class)
                    
                mic_mol = intsta_mic_mol_lst[class_index]
                
                species_index_lst = temp[temp['CLASS'] == lipid_class].index.tolist()
                    
                intsta_row = find_intsta_index(lipid_class, intsta_class_lst, no_intsta_class_lst, add_intsta_class_lst, intsta_df)
                    
                for i in range(sum(rep_lst)):
                        
                    nums1 = temp['PRECURSORINTENSITY[s'+str(i+1)+']'].iloc[species_index_lst]
                        
                    num2 = intsta_df['PRECURSORINTENSITY[s'+str(i+1)+']'].iloc[intsta_row]
                        
                    norm_nums = compute_normalized_auc(nums1, num2, mic_mol)
                        
                    temp['normalized_AUC[s'+str(i+1)+']'].iloc[species_index_lst] = norm_nums
                        
            elif lipid_class in no_intsta_class_lst:
                    
                class_index = no_intsta_class_lst.index(lipid_class)
            
                alt_lipid_class = add_intsta_class_lst[class_index]
                    
                alt_class_index = intsta_class_lst.index(alt_lipid_class)
                    
                mic_mol = intsta_mic_mol_lst[alt_class_index]
                    
                species_index_lst = temp[temp['CLASS'] == lipid_class].index.tolist()
                    
                intsta_row = find_intsta_index(lipid_class, intsta_class_lst, no_intsta_class_lst, add_intsta_class_lst, intsta_df)
                    
                for i in range(sum(rep_lst)):
                        
                    nums1 = temp['PRECURSORINTENSITY[s'+str(i+1)+']'].iloc[species_index_lst]
                        
                    num2 = intsta_df['PRECURSORINTENSITY[s'+str(i+1)+']'].iloc[intsta_row]
                        
                    norm_nums = compute_normalized_auc(nums1, num2, mic_mol)
                    
                    temp['normalized_AUC[s'+str(i+1)+']'].iloc[species_index_lst] = norm_nums
                        
        temp.drop(['PRECURSORINTENSITY[s' + str(i+1) + ']' for i in range(sum(rep_lst))], axis=1, inplace=True)
        
        full_sample_lst = ['s'+str(i+1) for i in range(sum(rep_lst))]
        
        temp = temp[['old_index','SPECIES', 'FORMULA', 'MASS', 'CLASS']+
                ['normalized_AUC[' + sample + ']' for sample in full_sample_lst]]
        
        # keeping the naming convention consistent 
        
        for sample in full_sample_lst:
            
            temp.rename(columns={'normalized_AUC[' + sample + ']': 'PRECURSORINTENSITY[' + sample + ']'}, inplace=True)
            
            intsta_df.rename(columns={'normalized_AUC[' + sample + ']': 'PRECURSORINTENSITY[' + sample + ']'}, inplace=True)
        
        confirm = st.checkbox('Confirm the inputs & view the normalized data')
        
        if confirm:
                        
            st.write(temp)
            
            csv_download = convert_df(temp)
                            
            st.download_button(
                
                                label="Download Data",
                                
                                data=csv_download,
                                
                                file_name='normalized_data.csv',
                                
                                mime='text/csv',
                                
                                key='normalized_data')
            
        norm_ans = st.radio('Would you like to use the normalized dataset for the rest of the analysis?', ['Yes', 'No'], 0)
        
        if norm_ans == 'Yes':
            
            dataframe = temp
        
        else:
            
            dataframe = X
        
        return dataframe
    
    
    
    
    
# function that receives the user input on which internal standards to use for classes with no internal standards 
def pick_intsta(no_intsta_class_lst, intsta_class_lst): 
            
        add_intsta_lst = []
            
        for lipid_class in no_intsta_class_lst:
                
            add_intsta_lst.append(build_intsta_selectbox(lipid_class, intsta_class_lst)) # function defined below 
                
        return add_intsta_lst
    
    
    
    
    
# function that creates a select box including the existing list of internal standards 
def build_intsta_selectbox(lipid_class, intsta_class_lst):  
        
        added_intsta_class = st.selectbox('Pick a internal standard for ' + lipid_class + ' species', intsta_class_lst)
                
        return added_intsta_class
    
    
    
    
    
# function for letting the user input the amount of IS 
def build_intsta_mic_mol_lst(intsta_class_lst):
        
        intsta_mic_mol_lst = []
        
        for lipid_class in intsta_class_lst:
            
            intsta_mic_mol_lst.append(build_mic_mol_input_box(lipid_class))
        
        return intsta_mic_mol_lst
    
    
    
    
    
# function for letting the user input the amount of IS
def build_mic_mol_input_box(lipid_class):

        mic_mol = st.number_input('Enter the concentration of ' + lipid_class + ' internal standard species in micromole', 
                                  min_value = 0, max_value = 100000, value = 1, step = 1)
        
        return mic_mol
    
    
    
    
    
# function for finding the index of IS in intsta_df
def find_intsta_index(lipid_class, intsta_class_lst, no_intsta_class_lst, add_intsta_class_lst, intsta_df):
        
        if lipid_class in intsta_class_lst:
            
            intsta_row = intsta_df[intsta_df['CLASS'] == lipid_class].index.tolist()[0]
            
        elif lipid_class in no_intsta_class_lst:
            
            class_index = no_intsta_class_lst.index(lipid_class)
            
            alt_lipid_class = add_intsta_class_lst[class_index]
            
            intsta_row = intsta_df[intsta_df['CLASS'] == alt_lipid_class].index.tolist()[0]
        
        return intsta_row
    
    
    
    
    
# function useful for normalization of data 
def compute_normalized_auc(nums1, num2, mic_mol):
    
        nums1 = [num/num2*mic_mol for num in nums1]
        
        return nums1
    
    
    
    
    
# function that creates all volcano plots
def volcano_plot(X, cond_lst, rep_lst):
    
        temp = X.copy()
    
        # different options for imputing 
        
        impute_ans = st.radio('Pick one option',\
                              ['No imputing', 'Replace missing values in each sample by the minimum detected AUC in that sample'], 0)
            
        if impute_ans == 'Replace missing values in each sample by the minimum detected AUC in that sample':
        
            temp = impute_missing_value(temp)
        
        show_vol = st.checkbox("View Volcano Plots")
        
        if show_vol:
            
            rep_lst_agg = build_rep_lst_agg(rep_lst)
            
            cond_1 = st.selectbox('Pick the first condition', cond_lst, 0)
            
            cond_2 = st.selectbox('Pick the second condition', cond_lst, 1)
            
            for rep in [rep_lst[cond_lst.index(cond_1)], rep_lst[cond_lst.index(cond_2)]]:
                
                if rep < 2:
                    
                    st.error('At least two samples per condition are required for a meanigful analysis!')
                    
                    return None
                
            sample_lst_1 = build_sample_lst(cond_1, cond_lst, rep_lst_agg)
            
            sample_lst_2 = build_sample_lst(cond_2, cond_lst, rep_lst_agg)
            
            temp['p_val_' + cond_1 + '_' + cond_2] = \
                    temp[['PRECURSORINTENSITY[' + sample + ']'  for sample in sample_lst_1 + sample_lst_2]]\
                    .apply(lambda x: p_val_calculator(x[['PRECURSORINTENSITY[' + sample + ']'  for sample in sample_lst_1]], \
                                                      x[['PRECURSORINTENSITY[' + sample + ']'  for sample in sample_lst_2]]), axis=1)
                        
            temp['fc_' + cond_1 + '_' + cond_2] = \
                    temp[['PRECURSORINTENSITY[' + sample + ']'  for sample in sample_lst_1 + sample_lst_2]]\
                    .apply(lambda x: fc_calculator(x[['PRECURSORINTENSITY[' + sample + ']'  for sample in sample_lst_1]], \
                                                   x[['PRECURSORINTENSITY[' + sample + ']'  for sample in sample_lst_2]]), axis=1)
                
            lipid_class_lst = temp['CLASS'].value_counts().index.tolist()
                    
            selected_class = st.multiselect('Add or remove classes (up to 20 classes):', lipid_class_lst, lipid_class_lst[:1])
            
            if len(selected_class) > 20:
                        
                st.error('You can only compare up to 20 lipid classes at a time!')
                        
                return None
                        
            plot, vol_df = volcano_hover(temp, selected_class, cond_1, cond_2)
            
            st.bokeh_chart(plot) 
            
            csv_download = convert_df(vol_df)
                            
            st.download_button(
                
                                label="Download Data",
                                
                                data=csv_download,
                                
                                file_name='volcano_plot.csv',
                                
                                mime='text/csv')
                
        return
    
    
    
    
    
 # function that prepares a single volcano plot with hover tool 
def volcano_hover(dataframe, selected_class, cond_1, cond_2):
        
        unique_color_lst =[ 'red', 'blue', 'green', 'magenta', 'cyan', 'orange', 'black', 'pink', 'brown', 'yellow', 'purple', \
                            'gray', 'olive', 'chocolate', 'silver', 'darkred', 'khaki', 'skyblue', 'navy', 'orchid']
        
        fc = []
        
        pval = []
        
        class_lst = []
        
        color_lst = []
        
        species = []
        
        index_lst = []
        
        for lipid_class in selected_class:
            
            fc = fc + dataframe[dataframe['CLASS'] == lipid_class]['fc_' + cond_1 + '_' + cond_2].values.tolist()
            
            pval = pval + dataframe[dataframe['CLASS'] == lipid_class]['p_val_' + cond_1 + '_' + cond_2].values.tolist()
            
            color_lst = color_lst + [unique_color_lst[selected_class.index(lipid_class)] for i in range(len(dataframe[dataframe['CLASS'] == lipid_class]))]
                
            class_lst = class_lst + dataframe[dataframe['CLASS'] == lipid_class]['CLASS'].values.tolist()
                
            index_lst = index_lst + dataframe[dataframe['CLASS'] == lipid_class]['old_index'].values.tolist()
        
            species = species + dataframe[dataframe['CLASS'] == lipid_class]['SPECIES'].values.tolist()
            
        plot = figure(title='Volcano Plot', x_axis_label='Fold Change (' + cond_1 + '/' + cond_2 + ')', y_axis_label='q-value')
            
        vol_df = pd.DataFrame({"FC": fc, "qvalue": -np.log10(pval), "Species": species, "Class": class_lst, \
                               "Color": color_lst, "Index": index_lst})
        
        src = ColumnDataSource(vol_df)
        
        plot.scatter(x="FC", y="qvalue", legend_group='Class', color='Color', name='volcano', size = 4, source=src)
        
        plot.line(x=[i for i in range(-10, 11)], y=-np.log10(0.05), line_dash = 'dashed', color='black')
        
        plot.line(x=-1, y=[i for i in range(0, 9)], line_dash = 'dashed', color='black')
        
        plot.line(x=1, y=[i for i in range(0, 9)], line_dash = 'dashed', color='black')
            
        hover = HoverTool(tooltips = [('FC', '@FC'), ('p-value', '@qvalue'), ('Species', '@Species'), ('index', '@Index')], names=['volcano'])
            
        plot.add_tools(hover)
        
        plot.title.text_font_size = "15pt"
        
        plot.xaxis.axis_label_text_font_size = "15pt"
            
        plot.yaxis.axis_label_text_font_size = "15pt"
            
        plot.xaxis.major_label_text_font_size = "15pt"
            
        plot.yaxis.major_label_text_font_size = "15pt"
        
        return plot, vol_df
    
    
    
    
    
    # calculate fold change 
@st.cache
def fc_calculator(num_1, num_2): 
        
        non_zero_num_1 = [num for num in num_1 if num > 0]
        
        non_zero_num_2 = [num for num in num_2 if num > 0]
        
        if len(non_zero_num_1) > 0 and len(non_zero_num_2) > 0:
    
            fc = np.log2(np.mean(non_zero_num_1)/np.mean(non_zero_num_2))
            
        else:
            
            fc = None
        
        return fc
    
    
    
    
    # calculate p-value by running a T-test 
@st.cache
def p_val_calculator(num_1, num_2):
        
        non_zero_num_1 = [num for num in num_1 if num > 0]
        
        non_zero_num_2 = [num for num in num_2 if num > 0]
        
        if len(non_zero_num_1) > 0 and len(non_zero_num_2) > 0:
        
            t_value, p_value = stats.ttest_ind(non_zero_num_1,non_zero_num_2)
            
        else: 
            
            p_value = None
        
        return p_value
    
    
    
    
    
 # function for imputing missing values 
def impute_missing_value(X):
    
    full_sample_lst = ['s'+str(i+1) for i in range(sum(rep_lst))]
                
    for sample in full_sample_lst:
                
                lst = [ele for ele in X['PRECURSORINTENSITY[' + sample + ']'].values if ele > 0]
                
                impute_value = min(lst)
                
                X['PRECURSORINTENSITY[' + sample + ']'] = X['PRECURSORINTENSITY[' + sample + ']'].apply(lambda x: impute_value if x<=0 else x)
        
    return X

##########################################################################################################################################
# the main code of the app 

st.header("LipidXplorer Module")

st.markdown("""
            
            The following module allows the user to run lipidomics analysis on LipidXplorer datasets. 
            Start by uploading your dataset and completing the next steps on the side bar.
            
            """)
            
st.info("""
        
        A standard LipidXplorer dataset must have the following columns:
            
        SPECIES: the class that the lipid species belong to and its number of carbon atoms and double bonds

        CLASS: the class that the lipid species belong to

        MASS: the calculated mass of the lipid species

        FORMULA: the chemical formula of the lipid species

        ERROR: the difference between theoritical and calculated mass

        PRECURSORINTENSITY:s1.mzML, ..., PRECURSORINTENSITY:sN.mzML: Area Under the Curve (AUC) representing the relative 
        abundance of the lipid species in samples s1 to sN where N stands for the total number of the sampels
        
        """)

st.sidebar.subheader('Select Mode')
        
mode = st.sidebar.radio('', ['Only positive or negative mode', 'Merge positive and negative modes'])
        
st.sidebar.subheader("Upload Data")
        
if mode == 'Only positive or negative mode':
            
            lipid_xplorer = st.sidebar.file_uploader(label='Upload your LipidXplorer dataset', type=['csv'])
            
            if lipid_xplorer is not None:
            
                df, intsta_df = build_single_lipidxplorer_df(lipid_xplorer)
                
            else: 
                
                df = None
            
elif mode == 'Merge positive and negative modes':
            
            lipid_xplorer_pos = st.sidebar.file_uploader(label='Upload your LipidXplorer dataset in POSITIVE mode', type=['csv'])
            
            lipid_xplorer_neg = st.sidebar.file_uploader(label='Upload your LipidXplorer dataset in NEGATIVE mode', type=['csv'])
            
            if lipid_xplorer_pos is not None:
                
                pos_df, pos_intsta_df = build_single_lipidxplorer_df(lipid_xplorer_pos)
                
            if lipid_xplorer_neg is not None:
                
                neg_df, neg_intsta_df = build_single_lipidxplorer_df(lipid_xplorer_neg)
        
            if (lipid_xplorer_pos is not None) and (lipid_xplorer_neg is not None):
            
                df, intsta_df = build_merged_lipidxplorer_df(pos_df, pos_intsta_df, neg_df, neg_intsta_df)
                
            else:
                
                df = None

if df is not None:
            
            confirm_data, cond_lst, rep_lst, group_df, filtered_conds, passing_abundance_grade = build_sidebar(df)
            
            if confirm_data:
                    
                st.subheader("1) Data Cleaning, Exploration, Quality Check & Anomaly Detection")
                        
                st.markdown("""
                                    The Data Cleaning, Exploration, Quality Check & anomaly detection submodule allows the user to clean, filter
                                    and visualize the data, and runs anomaly detection tests on it.
                    
                                    """)
                    
                st.subheader("1.1) Clean, Filter & Explore Data")
                
                expand_raw_data = st.expander("Raw Data")
            
                with expand_raw_data:
                            
                    st.write('View the raw data:')
                
                    st.write(df)
                    
                expand_clean_data = st.expander('Cleaned Data')
                        
                with expand_clean_data:
                            
                    st.markdown("""
                            
                        The data cleaning process is a six steps process: 
                                
                        1) LipidCruncher removes the empty rwos and deletes the duplicated datapoints.  
                                
                        2) LipidCruncher deletes the datapoints that do not pass through the applied filter: 
                           either their absolute error value is larger than 5 or they have too many missing values.  
                            
                        3) LipidCruncher only keeps the relevant columns: "SPECIES", "FORMULA", "MASS", "CLASS",
                            "PRECURSORINTENSITY[s1]", ..., "PRECURSORINTENSITY[sN]". Where N is the total number of the samples.
                            
                        4) LipidCruncher adds a column named "old_index" to the cleaned dataset 
                            which refers to the index of the lipid species in the raw dataset.
                                
                        5) LipidCruncher extracts the internal standards lipid species and puts them in a separate dataset.
                            An internal standard is a lipid species that is added in a known constant concentration to the samples.
                            As the rate of spraying the lipids in shotgun method is not constant, calibration is required. 
                            Internal standard lipids can be used for calibration purposes.
                            
                        """)
                        
                    X, intsta_f = apply_filter(df, intsta_df, group_df, rep_lst, cond_lst, filtered_conds, passing_abundance_grade) # cleaned data
                    
                expand_filtered_data = st.expander("Filter Data Using BQC Samples")
                    
                with expand_filtered_data:
                    
                    st.markdown(""" 
                                
                        Filtering mass spectrometry data is challenging. However, the most reliable method to filter the data 
                        is using "Batch Quality Control (BQC)" samples. Here is how BQC samples are created: take equal amount 
                        of aliquot from each study sample and pool. Then, create multiple BQC samples by transferring the pooled 
                        mixture to new tubes. 
                        
                        BQC samples that are created as instructed are practically technical replicates. Technical replicates are 
                        expected to be approximately identical. Therefore, the vast majority of the lipid species must have a very
                        low CoV (i.e. CoV < 30%). The coefficient of variation (CoV) is defined as the ratio of the standard deviation to the mean.
                        It shows the extent of variability in relation to the mean of the population and is often expressed as a percentage.
                        The red line in the CoV plot is the line of "30% CoV".
                        
                        Here is how LipidCruncher filter data using BQC samples: lipid species with CoV lower than a set threshold (e.g. 30%) are kept
                        and the rest are removed. The readings with low CoV are highly reliable and the ones with a high CoV are less reliable. 
                        It seems reasonble to provide the user with the option to filter the data using BQC samples. 
                    
                        """)
                        
                    st.info("Creating BQC samples in any lipidomics experiment is a great practice.")
                    
                    X = plot_cov(X, cond_lst, rep_lst)
                    
                expand_hist = st.expander("View Distributions of AUC: Scan Data & Detect Atypical Patterns")
                
                with expand_hist:
                    
                    st.markdown("""

                        In a standard LipidSearch dataset, columns "MainArea[s1]" to "MainArea[sN]" correspond to Area 
                        Under the Curve (AUC) representing the relative abundance of the lipid species 
                        in samples s1 to sN. 
                        
                        To plot the histograms of AUC, LipidCruncher turns all 0 values (i.e. missing values) to 1 
                        and then log-transforms the whole dataset. This allows the user to visualize what portion of the 
                        values are missing without affecting the distribution (as 1 is orders of magnitude smaller than 
                        the minimum detection threshold by mass spectrometry).
                        
                        Visualize the distribution of AUC's in any of the replicates and look for atypical patterns (i.g. too many missing values):

                        """)
                
                    plot_hist(X, rep_lst, cond_lst) # histograms of AUC
                    
                st.subheader("1.2) Detect Anomalies")
            
                expand_corr = st.expander('Pairwise Correlation Analysis') 
                
                with expand_corr:
                    
                    st.markdown("""
                        
                        Typically, the AUC's of any sample is highly linearly correlated to those of its biological replicate
                        (i.e. correlation coefficient > 0.8). This linear correlation is expected to be even stronger for technical replicates 
                        (i.e. correlation coefficient > 0.9).
                        A sample that has a weak correlation with its biological replicates is an outlier. That sample might be an outlier because 
                        of the natural biological variance or an error during sample preparation.
                        
                        Run a correlation test to inspect the degree of correlation between any two biological or technical replicates:
                        
                        """)
                        
                    st.info("LipidCruncher removes the missing values before preforming the correlation test.")
            
                    plot_corr(X, cond_lst, rep_lst) # pairwise correlation plots 
                    
                expand_pca = st.expander('Principal Component Analysis (PCA)')
            
                with expand_pca:
                
                    st.markdown("""
                        
                        Principal Component Analysis, or PCA, is a dimensionality-reduction method that is often used to reduce the 
                        dimensionality of large data sets, by transforming a large set of variables into a smaller one that still contains 
                        most of the information in the large set.
                        
                        Typically, biological replicates are expected to cluster together with the exception of rare outliers that fall 
                        further away from the rest of the replicates. 
                        
                        A plot of the top two principal components (PC1 and PC2) against each other is the best way to inspect the clustering 
                        of different samples. This is because the PC's are ordered based on how much of the variability in the data 
                        they can explain (i.e. PC1 is the PC with the highest variance explained ratio, PC2 is the PC with the second highest 
                        variance expalined ratio and so on).
                        
                        Run PCA to inspect the clustering of different samples:
                        
                        """)
                        
                    st.info("LipidCruncher does NOT remove missing values before performng PCA analysis.")

                    plot_pca(X, rep_lst, cond_lst) # PCA analysis
                    
                st.subheader("2) Data Analysis & Hypothesis Testing")
                
                st.markdown("""
                                    The Data Analysis & Hypothesis Testing submodule allows the user to run statistical analysis on 
                                    the data and test their hypothesis. 
                    
                                    """)
                                    
                st.subheader("2.1) Remove Anomalies")
            
                expand_cleaned_data = st.expander("Remove Anomalies (and/or Blank Samples)")
                
                with expand_cleaned_data:

                    X, rep_lst, cond_lst = remove_bad_sample(X, rep_lst, cond_lst) # cleaned data
                    
                st.subheader('2.2) Normalize Data')
                
                expand_norm_data = st.expander('Normalized Data')
                
                with expand_norm_data:
                    
                    st.info('An internal standard (IS) is substance that is similar to the analyte that is added in a constant amount to the samples.\
                            In order to calibrate the variation in the data, LipidCruncher normalizes the data using using internal standards.')
                            
                    st.write('The following formula is used for data normalization:')
                            
                    latext = r'''
                                    
                            $$ 
                            Concentration (analyte) = \frac{AUC(analyte)}{AUC(IS)} \times Concentration(IS) 
                            $$  
                            
                            '''
                    st.write(latext)
                                        
                    X = normalize_auc(X, intsta_df, rep_lst)
                    
                st.subheader("2.3) Analyze Data")
                    
                expand_vol_plot = st.expander("Volcano Plots")
                
                with expand_vol_plot:
                    
                    st.markdown("""
                                
                                In statistics, a volcano plot is a type of scatter-plot that is used to quickly identify changes in \
                                large data sets composed of replicate data. It plots significance versus fold-change on the y and x axes, respectively. \
                                A volcano plot combines a measure of statistical significance from a statistical test (e.g., a p value from a T-test) \
                                with the magnitude of the change, enabling quick visual identification of those data-points that display large magnitude\
                                changes that are also statistically significant (datapoints at the top left and top right quadrant).
                                
                                Below, q-value (i.e. -log10(p-value)) of each lipid species is plotted versus the fold change of that species. \
                                The p-value is computed from a two-sample T-test and the fold change is computed from the following formula:
                                
                                """)
                        
                    latext = r'''
                            
                    $$ 
                    Fold Change = log2(\frac{Mean AUC(Condition 1)}{Mean AUC(Condition 2)})
                    $$  
                    
                    '''
                    
                    st.write(latext)
                    
                    st.info("The outcome of a T-test is most accurate when both samples that are being compared have the same size (i.e. \
                            equal number of datapoints). Therefore, to create a volcano plot, the best practice is to impute missing values.")
                        
                    st.warning("Your choice here only affects the volcano plot, not the rst of the submodule.")
                    
                    X.reset_index(inplace = True)
                
                    volcano_plot(X, cond_lst, rep_lst)