import csv
import streamlit as st
import pandas as pd
import numpy as np
from scipy import stats
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, HoverTool, Whisker
from bokeh.transform import dodge
import base64

from bokeh.models import BasicTickFormatter, Range1d


def app():
    ##########################################################################################################################################
    def build_lipidsearch_df(lipid_search): # function that creates a pandas dataframe out of the uploaded LipidSearch dataset
    
        with open(lipid_search.name, "wb") as f:
            
            f.write(lipid_search.getbuffer()) # saves the file in the app directory
            
        with open(lipid_search.name, newline='') as f: # counts the number of lines at the beginning of the file that need to be removed 
            
            reader = csv.reader(f)
            
            counter = 0 # counts the number of lines before the line that starts with 'Rej
            
            counter_rej = 0 # 0 means column 'Rej' was not found, 1 means column 'Rej was found 
            
            for row in reader:
                
                if len(row) > 0 and 'Rej' in row[0][0:3]: # All the lines before the line that starts with 'Rej  should be removed
                    
                    counter_rej = counter_rej + 1
                    
                    break
                
                else:
                    
                    counter = counter+1
                        
        if counter_rej == 0:
            
            st.sidebar.error('This is not a valid LipidSearch dataset!')
            
            return None
        
        else:
        
            with open(lipid_search.name,'r') as f: # creates a new file without the first few extra lines 
            
                with open('trimmed_'+lipid_search.name,'w') as f1:
                
                    for i in range(counter):
                    
                        next(f) # skip first line
                    
                    for line in f:
                    
                        f1.write(line)
                    
            
            if ".txt" in lipid_search.name:
                
                df = pd.read_csv('trimmed_'+lipid_search.name, delimiter = "\t", error_bad_lines=False)
                
                if ('LipidMolec' and 'Class' and 'Calc Mass' and 'BaseRt' and 'MainArea[s1]' and 'MainGrade[s1]') in df.columns.values.tolist():
                
                    return df
                
                else:
                    
                    return None
                
            elif ".csv" in lipid_search.name:
                
                df = pd.read_csv('trimmed_'+lipid_search.name, error_bad_lines=False)
                
                if ('LipidMolec' and 'Class' and 'Calc Mass' and 'BaseRt' and 'MainArea[s1]' and 'MainGrade[s1]') in df.columns.values.tolist():
                    
                    return df
                
                else:
                    
                    return None
    
    
    
    
        
    # function that builds the side bar for LipidSearch
    def build_sidebar_lipid_search(df):   
            
        st.sidebar.subheader("Define Experiment")
            
        n_cond = st.sidebar.number_input('Enter the number of conditions',
                                         min_value = 1, max_value= 20, value = 1, step = 1)
            
        cond_lst, rep_lst = build_cond_rep_lst(n_cond)
        
        # finding the total number of the samples 
        
        counter = -1 # initiates at -1 instead of zero because of the extra control sample  
        
        for column in df.columns.values.tolist():
            
            if 'MainArea[' in column:
                
                counter = counter + 1
                
        # if inputs are not valid 
        
        if counter != sum(rep_lst):
            
            st.sidebar.error('The inputs are incomplete and/or inaccurate!')
            
            # assign arbitrary values to vaiables to avoid an error 
            
            passing_letter_grade = 0
            
            pval_tresh = 0
            
            passing_pval_grade = 0 
            
            filter_mode = ''
            
            filtered_conds = ''
            
            passing_abundance_grade = ''
            
            name_df = ''
            
            missing_ans = ''
            
            # assing 'False' to confirm_data, so, the algorithm stops  
            
            confirm_data = False
            
        # if the inputs are valid
        
        else:
            
            st.sidebar.subheader('Group Samples')
        
            group_df = group_samples(df, cond_lst, rep_lst)
                
            st.sidebar.subheader("Apply Built-in Filters")
            
            passing_letter_grade = st.sidebar.number_input('Enter the minimum required MainGrade score (i.e. how many A and/or B grades?)'
                                                           , min_value = 0, max_value = sum(rep_lst), value = 1, step = 1)
        
            pval_tresh = 0.001
            
            passing_pval_grade = st.sidebar.number_input(
                    'Enter the minimum required PValue score '
                    '(i.e. how many replicates need to have a p-value smaller than the threshold?) '
                                                         , min_value = 0, max_value = sum(rep_lst), value = 1, step = 1)
            
            st.sidebar.subheader('Apply Additional Filters')
            
            filter_mode, filtered_conds, passing_abundance_grade = input_filter_mode_cond_grade(cond_lst, rep_lst)
            
            st.sidebar.subheader("Confirm Inputs")

        
            st.sidebar.markdown('''
                            
                                LipidSearch names the samples as following: s1, s2, ..., sN. Where N is the total number of the samples. 
                                However, there are cases where one or more samples are missing. An example is: s1, s2, s4, ..., s10.
                                In that case, LipidCruncher updates the sample names as following: s1, s2, s3, ..., s9. 
                            
                                ''')
                                
            missing_ans = st.sidebar.radio('Do you have any missing samples?', ['Yes', 'No'], 1)
            
            name_df = update_sample_name(group_df)
            
            if missing_ans == 'Yes':
                
                st.sidebar.write('The table below shows the updated sample names.')
        
                st.sidebar.write(name_df)
        
            st.sidebar.write("Now, confirm your inputs:")
            
            st.sidebar.write("There are a total of "+str(sum(rep_lst))+" samples.")
            
            for cond in cond_lst:
            
                build_rep_cond_pair(cond, cond_lst, rep_lst) # function defined below 
            
            confirm_data = st.sidebar.checkbox("Confirm the inputs by checking this box")
            
        return confirm_data, missing_ans, name_df, passing_letter_grade, pval_tresh, passing_pval_grade, filter_mode, \
               filtered_conds, passing_abundance_grade, cond_lst, rep_lst
    
    
    
    
        
    # function that builds cond_lst and rep_lst objects (conditions list and number of replicates list)
    def build_cond_rep_lst(n_cond): 
    
        '''
        For example: if we have two conditions, WT and KO, and each have two corresponding replicates, 
        the cond_lst = [WT, KO] and the rep_lst = [2, 2]
        
        Used for both LipoidSearch and LipidXplorer 
        
        '''
            
        cond_lst = []
            
        rep_lst = []
            
        for i in range(1, n_cond+1):
                
            cond_lst.append(cond_text_input(i)) # function defined below 
                
            rep_lst.append(rep_number_input(i)) # function defined below
                
        return cond_lst, rep_lst
    
    
    
    
            
    # function that creates a text input box to enter the label for each condition
    def cond_text_input(i):  
    
        # Used for both LipidSearch and LipidXplorer 
        
        cond = st.sidebar.text_input('Create a label for condition #'+str(i)+' (e.g. WT or KO)')
                
        return cond
    
    
    
    
            
    # function that creates a text input to enter the number of replicates corresponding each condition
    def rep_number_input(i):  
            
        # Used for both LipidSearch and LipidXplorer
        
        rep = st.sidebar.number_input('Enter the number of replicates for condition #'+str(i), 
                                              min_value = 1, max_value = 1000, value = 1, step = 1)
                
        return rep  
    
    
    
    
              
    # function to build the aggregate list of rep_lst elements
    def build_rep_lst_agg(rep_lst):  
            
        '''
        For example, if rep_lst = [2, 4, 5], the rep_lst_agg = [2, 6, 11]
        
        Used for both LipidSearch and LipidXplorer
        '''
            
        rep_lst_agg = [] # the aggregate list of replicates 
            
        for i in range(len(rep_lst)):
                
            rep_lst_agg.append(sum(rep_lst[0:i+1]))
                
        return rep_lst_agg
    
    
    
    
    
    # function to build the list of replicates (not the number of replicates) for each condition
    def build_sample_lst(cond, cond_lst, rep_lst_agg):  
    
        '''
        For example, if cond_lst = [WT, KO, HA] and rep_lst = [4, 4, 2], then rep_lst_agg = [4, 8, 10].
        In this case, for cond = WT, the sample_lst = [s1, s2, s3, s4], for cond = KO, the sample_lst = 
        [s5, s6, s7, s8] and for cond = HA, the sample_lst = [s9, s10].
        
        Used for both LipidSearch and LipidXplorer
        '''
        
        if cond_lst.index(cond) == 0:
                
            sample_lst = ['s'+str(i+1) for i in range(rep_lst_agg[cond_lst.index(cond)])]
                
        else:
                
            sample_lst = ['s'+str(i+1) for i in range(rep_lst_agg[cond_lst.index(cond)-1], \
                                                                      rep_lst_agg[cond_lst.index(cond)])]
                    
        return sample_lst
    
    
    
    
    
    # pairs the replicates with their corresponding condition
    def build_rep_cond_pair(cond, cond_lst, rep_lst):  
    
        '''
        For example, if cond_lst = [WT, KO] and rep_lst = [2, 2], the function returns:
            s1-s2 correspond to WT
            s3-s4 correspond to KO
            
        Used for both LipidSearch and LipidXplorer
        '''
                
        rep_lst_agg = build_rep_lst_agg(rep_lst) # function defined above 
            
        sample_lst = build_sample_lst(cond, cond_lst, rep_lst_agg) # function defined above 
            
        if len(sample_lst) == 1:
                
            return st.sidebar.write("- "+sample_lst[0]+" corresponds to "+cond)
            
        else:
        
            return st.sidebar.write("- "+sample_lst[0]+"-"+sample_lst[-1]+" correspond to "+cond)
    
    
    
    
    
    # function that computes the total main grade for each lipid species
    def letter_grade_calculator(letters, passing_letter_grade):  
        
                    counter = 0
                    
                    for letter in letters:
                        
                        if letter == 'A' or letter == 'B':
                            
                            counter = counter+1
                            
                    if counter >= passing_letter_grade: 
                        
                        return 1  
                    
                    else:
                        
                        return 0
        
    
    
    
    
    # function that finds out if the lipid species has enough p-values below the trshold
    def pval_grade_calculator(pvals, passing_pval_grade, pval_tresh):
        
                    counter = 0
                    
                    for pval in pvals:
                        
                        if pval < pval_tresh:
                            
                            counter = counter+1
                            
                    if counter >= passing_pval_grade: 
                        
                        return 1  
                    
                    else:
                        
                        return 0
                    
    
    
    
    
    #function that computes the total abundance grade for each lipid species
    def abundance_level_calculator(intensities, passing_abundance_grade): 
        
        '''
        For example: if the passing abundance grade is 3, and condition A has a at least 3 non-zero values, 
        the function returns 1, otherwise it returns 0.
        '''
        
        counter = 0
                    
        for intensity in intensities:
            
            #st.write(float(intensity))
                        
            if float(intensity) > 0:
                            
                counter = counter+1
                            
        if counter >= passing_abundance_grade:
                        
            return 1  
                    
        else:
                        
            return 0
    
    
    
    
    
    # function that removes the lipid species that have an abundance grade lower than the minimum abundance grade
    def build_filter(X, rep_lst, cond_lst, filter_mode, filtered_conds, passing_abundance_grade):
        
        rep_lst_agg = build_rep_lst_agg(rep_lst)  
        
        for cond in filtered_conds:
            
            sample_lst = build_sample_lst(cond, cond_lst, rep_lst_agg)
            
            # adding the 'abundance grade' column for each filtered condition 
            
            X['abundance_grade_'+cond] = X[['MainArea['+sample+']' for sample in sample_lst]]\
                .apply(lambda x: abundance_level_calculator(x, passing_abundance_grade),axis=1)
            
        X['abundance_grade'] = X[['abundance_grade_'+cond for cond in filtered_conds]].apply(lambda x: sum(x), axis=1) # adds the 'total abundance grade' column 
        
        if filter_mode == 'At least one of the follwong conditions':
        
            X = X.loc[X['abundance_grade'] > 0 ] # keeps only lipid species with non-zero 'total abundance grade'  
                
        elif filter_mode == 'Each of the following conditions':
                
            X = X.loc[X['abundance_grade'] == len(filtered_conds) ] 
            # keeps only lipid species with 'total bundance grade' equal to the number of selected conds for filtering
        
        X.drop(['abundance_grade_'+cond for cond in filtered_conds]+['abundance_grade'], axis=1, inplace=True) # drops the abundance grade column
        
        return X
    
    
    
    
    
    # function to receive the following inputs from the user: filter mode, filtered conditions and the minimum renumber of non-zero values required 
    def input_filter_mode_cond_grade(cond_lst, rep_lst):
        
        filter_mode = st.sidebar.radio('Apply filter to ', 
                                       ['Each of the following conditions', 
                                        'At least one of the follwong conditions'])
            
        filtered_conds = st.sidebar.multiselect('Specify which conditions to apply the filter to (remove/add conditions)', cond_lst, cond_lst)
            
        if len(filtered_conds) == 0:
                
            st.sidebar.error('You must choose at least one condition. Otherwise, no filter is applied!')
                
            filtered_conds = [cond_lst[0]]
                
            passing_abundance_grade = 0
                
        else:
            
            filtered_reps = [] # includes the number of reps corresponding to the selected_conds list 
            
            for cond in filtered_conds:
                
                index = cond_lst.index(cond)
                
                filtered_reps.append(rep_lst[index])
            
            passing_abundance_grade = st.sidebar.number_input("Enter the minimum required number of non-zero values "
                                                              , min_value = 0, max_value = np.min(filtered_reps), value = 0, step = 1)
            
            if filter_mode == 'Each of the following conditions':
            
                st.sidebar.write('Each of the selected conditions must have at least ' 
                                 + str(passing_abundance_grade) + ' non-zero value(s).')
                
            elif filter_mode == 'At least one of the follwong conditions':
                
                st.sidebar.write('At least one of the selected conditions must have at least ' 
                                 + str(passing_abundance_grade) + ' non-zero value(s).')
                
            st.sidebar.write('Selected conditions:')
            
            for cond in filtered_conds:
                
                st.sidebar.write('- ' + cond)
        
        return filter_mode, filtered_conds, passing_abundance_grade 
    
    
    
    
    
    # function to log transform the abundance columns of the df
    def log_transform_df(X, rep_lst, full_sample_lst):  
    
        temp = X.copy()
        
        auc = ['MainArea[' + sample + ']' for sample in full_sample_lst]
         
        temp[auc]=temp[auc].mask(temp[auc]<=0).fillna(1)
        
        temp[auc] = temp[auc].apply(lambda x: np.log10(x), axis=0)
        
        return temp
    
    
    
    
    
    # function to apply filters to the LipidSearch data
    def apply_filter_lipid_search(df, rep_lst, cond_lst, missing_ans, name_df, filter_mode, filtered_conds, passing_abundance_grade):  
        
            #first filter
            df = df.loc[df['Rej'] == 0] # removes the datapoint if 'Rej' = 1
        
            # second filter 
            total_reps = sum(rep_lst) # total number of all replicates 
            
            # if there are missing samples (i.g. s1, s2, s4, ..., sN where s3 is missing)
            
            if missing_ans == 'Yes':
                
                # updating the column names 
                
                for (sample_1, sample_2) in zip(name_df['old name'], name_df['updated name']):
                    
                        df.rename(columns={'MainArea['+sample_1+']': 'MainArea['+sample_2+']'}, inplace=True)
                        
                        df.rename(columns={'MainGrade['+sample_1+']': 'MainGrade['+sample_2+']'}, inplace=True)
                        
                        df.rename(columns={'APValue['+sample_1+']': 'APValue['+sample_2+']'}, inplace=True)
                
            df['total_letter_grade'] = df[['MainGrade[s'+str(i+1)+']' for i in range(total_reps)]]\
            .apply(lambda x: letter_grade_calculator(x, passing_letter_grade),axis=1) # adds the 'total letter grade' column
                    
            df = df.loc[df['total_letter_grade'] == 1] # removes the lipid species that do not pass the main grade filter 
                
            # third filter 
            df['total_pval_grade'] = df[['APValue[s'+str(i+1)+']' for i in range(total_reps)]]\
            .apply(lambda x: pval_grade_calculator(x, passing_pval_grade, pval_tresh),axis=1) # adds the 'total p-val grade' column
                    
            df = df.loc[df['total_pval_grade'] == 1] # removes the lipid species that do not pass the p-val grade filter
                 
            # extracting the columns required for analysis  
            X = df[['LipidMolec', 'Class', 'Calc Mass', 'BaseRt'] + ['MainArea[s'+str(i+1)+']' for i in range(total_reps)]]
            
            #second filter: minimum abundance grade 
        
            X = build_filter(X, rep_lst, cond_lst, filter_mode, filtered_conds, passing_abundance_grade)
            
            X.reset_index(inplace=True)
            
            X.rename(columns={'index': 'old_index'}, inplace=True)
            
            # removing bad samples 
            
            rep_lst, cond_lst, full_sample_lst, X, r_sample = remove_bad_samples(cond_lst, rep_lst, X)
            
            # X_impute is a dataframe which includes the info on which datapoint is imputed
            
            X_impute = add_impute_column(X, full_sample_lst, rep_lst, cond_lst)
            
            # different options for imputing 
            
            X, X_impute, impute_ans = impute_missing_value(X, X_impute, full_sample_lst)
            
            view_clean_data = st.checkbox('View the cleaned data')
            
            if view_clean_data:
                
                st.write('View the cleaned data in the conventional format:')
                
                st.write(X)
                
                csv_downloader(X, 'cleaned_data')
            
                st.write('------------------------------------------------------------------------------------------------')
            
                st.write('View the cleaned data in the log-transformed format:')
                
                log_X = log_transform_df(X, rep_lst, full_sample_lst)
            
                st.write(log_X)
            
                csv_downloader(log_X, 'log_transformed_cleaned_data')
            
            return X, X_impute, rep_lst, cond_lst, full_sample_lst, impute_ans
        
        
    
        
    
    # function that adds info on which data points need imputation 
    def add_impute_column(X, full_sample_lst, rep_lst, cond_lst):
        
        X_impute = X.copy()
            
        for sample in full_sample_lst:
                
                X_impute[sample + '_impute'] = X_impute['MainArea[' + sample + ']'].apply(lambda x: 0 if x<=0 else 1)
                
        rep_lst_agg = build_rep_lst_agg(rep_lst)
                
        for cond in cond_lst:
                
            sample_lst = update_sample_lst(cond, cond_lst, rep_lst_agg, full_sample_lst)
                
            index = cond_lst.index(cond)
                
            X_impute[cond + '_impute'] = X_impute[[sample + '_impute' for sample in sample_lst]].\
                apply(lambda x: 0 if sum(x) != rep_lst[index] else 1, axis = 1)
                    
        for cond_1 in cond_lst:
                
            for cond_2 in cond_lst:
                    
                X_impute[cond_1 + '_' + cond_2 + '_impute'] = X_impute[[cond_1 + '_impute', cond_2 + '_impute']].\
                    apply(lambda x: 1 if sum(x) == 2 else 0, axis = 1)
                
        X_impute.drop([sample + '_impute' for sample in full_sample_lst], axis=1, inplace = True)
        
        return X_impute
    
    
    
    
    
    # function for imputing missing values 
    def impute_missing_value(X, X_impute, full_sample_lst):
        
        impute_ans = st.radio('Select a method for imputing missing values',\
                              ['No imputing', 'Replace missing values in each sample by the minimum detected AUC in that sample'], 0)
                
        #auc = ['MainArea[' + sample + ']' for sample in full_sample_lst]
                
        if impute_ans == 'Replace missing values in each sample by the minimum detected AUC in that sample':
                
            for sample in full_sample_lst:
                
                lst = [ele for ele in X['MainArea[' + sample + ']'].values if ele > 0]
                
                impute_value = min(lst)
                
                X['MainArea[' + sample + ']'] = X['MainArea[' + sample + ']'].apply(lambda x: impute_value if x<=0 else x)
                    
                X_impute['MainArea[' + sample + ']'] = X_impute['MainArea[' + sample + ']'].apply(lambda x: impute_value if x<=0 else x)
        
        return X, X_impute, impute_ans 
    
    
    
    
    
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
                
                if ('FORMULA' and 'MASS' and 'ERROR') in df.columns.values.tolist():
                
                    return df, intsta_df
                
                else: 
                    
                   st.sidebar.error('This is not a valid LipidXplorer dataset!')
                
                   return None, None 
            
            else:
                
                st.sidebar.error('This is not a valid LipidXplorer dataset!')
                
                return None, None
            
            
            
            
            
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

        pos_cols = ['SPECIES', 'MOLSPECIES','FORMULA', 'MASS', 'ERROR', 'CLASS']
        
        neg_cols = ['SPECIES', 'MOLSPECIES','FORMULA', 'MASS', 'ERROR', 'CLASS']
        
        # if both datasets are uploaded 
        
        if (pos_df is not None) and (neg_df is not None):
            
            # dealing with classes that are detected but not included in pos_class + neg_class
            
            all_class = pos_class + neg_class
            
            pos_df_unique_class = pos_df['CLASS'].unique()
            
            pos_class_rest = [x for x in pos_df_unique_class if (x not in all_class) and (str(x) != 'nan')]
            
            if len(pos_class_rest) > 0:
            
                st.sidebar.warning('The following additional class(es) are detected in the positive mode. \
                                    Use the drop down menu to add additional class(es) to the list of the positive mode classes.')
            
                for lipid_class in pos_class_rest:
                    
                    st.sidebar.write(lipid_class)
                    
                pos_class_add = st.sidebar.multiselect('Add classes to the list of positive mode classes', pos_class_rest)
                
                pos_class = pos_class + pos_class_add
                
            neg_df_unique_class = neg_df['CLASS'].unique()
            
            neg_class_rest = [x for x in neg_df_unique_class if (x not in all_class) and (x not in pos_class) and (str(x) != 'nan')]
            
            if len(neg_class_rest) > 0:
            
                st.sidebar.warning('The following additional class(es) are detected in the negative mode. \
                                    Use the drop down menu to add additional class(es) to the list of the negative mode classes.')
            
                for lipid_class in neg_class_rest:
                    
                    st.sidebar.write(lipid_class)
                    
                neg_class_add = st.sidebar.multiselect('Add classes to the list of the negative mode classes', neg_class_rest)
                
                neg_class = neg_class + neg_class_add
                
            # picking up the relevant classes from each mode and merging the two datasets 
        
            pos_df = pos_df.loc[pos_df['CLASS'].isin(pos_class)]
            
            #pos_intsta_df.dropna(how='all', inplace=True) # removes empty rows
            
            #pos_intsta_df['CLASS'] = pos_intsta_df['SPECIES'].apply(lambda x: x.split(' ')[0])
            
            pos_intsta_df = pos_intsta_df.loc[pos_intsta_df['CLASS'].isin(pos_class)]
            
            for column in pos_df:
                
                if 'PRECURSORINTENSITY' in column:
                    
                    pos_cols.append(column)
            
            pos_df = pos_df[pos_cols]
            
            pos_intsta_df = pos_intsta_df[pos_cols]
            
            neg_df = neg_df.loc[neg_df['CLASS'].isin(neg_class)]
            
            #neg_intsta_df.dropna(how='all', inplace=True) # removes empty rows
            
            #neg_intsta_df['CLASS'] = neg_intsta_df['SPECIES'].apply(lambda x: x.split(' ')[0])
            
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
    def build_sidebar_lipid_xplorer(df):
        
        st.sidebar.subheader("Define Experiment")
            
        n_cond = st.sidebar.number_input('Enter the number of conditions',
                                         min_value = 1, max_value= 20, value = 1, step = 1)
            
        cond_lst, rep_lst = build_cond_rep_lst(n_cond)
        
        # finding the total number of the samples 
        
        counter = 0 
        
        for column in df.columns.values.tolist():
            
            if 'PRECURSORINTENSITY' in column:
                
                counter = counter + 1
                
        # if inputs are not valid
                
        if counter != sum(rep_lst):
            
            st.sidebar.error('The inputs are incomplete and/or inaccurate!')
            
            group_df = []
            
            filter_mode = ''
            
            filtered_conds = []
            
            passing_abundance_grade = 0
            
            confirm_data = False
            
        # if inputs are valid 
            
        else:
        
            st.sidebar.subheader('Group Samples')
        
            group_df = group_samples(df, cond_lst, rep_lst)
                
            st.sidebar.subheader("Apply Filters")
        
            filter_mode, filtered_conds, passing_abundance_grade = input_filter_mode_cond_grade(cond_lst, rep_lst)
            
            st.sidebar.subheader("Confirm Inputs")
        
            st.sidebar.markdown('''
                            
                                LipidSearch names the samples as following: s1, s2, ..., sN. Where N is the total number of the samples. 
                                However, LipidXplorer allows the user to name the samples.
                                To remain consistent, LipidCruncher unites the naming convention by following the LipidSearch protocol. 
                                Check the following table to see how the sample names have been updated:
                            
                                ''')
                            
            name_df = update_sample_name(group_df)
        
            st.sidebar.write(name_df)
        
            st.sidebar.write("Now, confirm your inputs:")
        
            st.sidebar.write("There are a total of "+str(sum(rep_lst))+" samples.")
            
            for cond in cond_lst:
            
                build_rep_cond_pair(cond, cond_lst, rep_lst) # function defined below 
            
            confirm_data = st.sidebar.checkbox("Confirm the inputs by checking this box")
            
        return confirm_data, cond_lst, rep_lst, group_df, filter_mode, filtered_conds, passing_abundance_grade
    
    
    
    
    
    # function to build a df with two columns: sample names and corresponding conditions
    def build_group_df(df, cond_lst, rep_lst):   
        
        extensive_cond_lst = [] # list includes the cond correponding to each rep - length equal to the total number of reps 
        
        for cond in cond_lst:
            
            index = cond_lst.index(cond)
            
            extensive_cond_lst += [cond for i in range(rep_lst[index])]

        group_dict = {'sample name' : [], 'condition' : []}
        
        column_lst = df.columns.values.tolist()
        
        counter = 0
    
        for column in column_lst:
    
            if (dataset_type == 'LipidXplorer') and ('PRECURSORINTENSITY' in column):
    
                a = column.split(':')
                
                b = a[1].split('.')
    
                group_dict['sample name'].append(b[0])  
                
                counter = counter+1
                    
            elif (dataset_type == 'LipidSearch') and ('MainArea[' in column) and (column != 'MainArea[c]'):
            
                group_dict['sample name'].append(column[9 : -1])
        
        # puts the columns in a sorted order  
        
        if dataset_type == 'LipidSearch':
            
            group_dict['sample name'] = [int(item[1:]) for item in group_dict['sample name']]
            
            group_dict['sample name'].sort()
            
            group_dict['sample name'] = ['s'+str(item) for item in group_dict['sample name']]
            
        for cond in extensive_cond_lst:
            
            group_dict['condition'].append(cond)
            
        group_df = pd.DataFrame.from_dict(group_dict) # df including sample_name and corresponding condition 
        
        return group_df 
    
    
    
    
    
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
            
            # ask user to put in the correct info 
            
            st.sidebar.write("""
                         
                         Samples must be properly grouped together. 
                         
                         Consider the following example:  
                             
                         *Example*: assume there are two conditions, A and B, and they have 3 and 5 
                         corresponding replicates, respectively. In this case, samples indexed from 0 to 2 must correspond to A and 
                         samples indexed from 3 to 7 must correspond to B. However, it turns out the samples are not properly grouped together:
                         samples indexed 0, 2 and 7 belong to A and samples indexed 1, 3, 4, 5, 6 belong to B. 
                         To fix this issue, you have to enter the following in the text box below: 0, 2, 7, 1, 3-6
                         
                         Use the above table and example to properly group your samples together. 
                         
                         """)
            
            group_sample_txt = st.sidebar.text_input('Enter the correct order of the samples:')
            
            # check the validity of the user input
            
            txt_validity_test = []
            
            test = group_sample_txt.split(',')
            
            for item in test:
                
                if '-' not in item:
                    
                    txt_validity_test.append(item)
                    
                else:
                    
                    sub_test = item.split('-')
                    
                    txt_validity_test.append(sub_test[0])
                    
                    txt_validity_test.append(sub_test[1])
                    
            
            counter = 0
            
            for item in txt_validity_test:
                
                if item.replace(" ", "") in [str(i) for i in range(sum(rep_lst))]:
                    
                    counter = counter + 1
                    
                    
            if counter == len(txt_validity_test):
                
                # change the order of the samples based on the user input
                
                group_sample_txt = group_sample_txt.split(',')
                
                group_sample_lst = []
            
                for item in group_sample_txt:
                
                    if '-' not in item:
                    
                        group_sample_lst.append(group_df['sample name'].iloc[int(item)])
                    
                    else:
                    
                        subgroup = item.split('-')
                    
                        group_sample_lst += [group_df['sample name'].iloc[int(ele)] for ele in range(int(subgroup[0]), int(subgroup[1])+1)]
                        
                st.sidebar.write('Check the updated table below to make sure the samples are properly grouped together.')        
            
                group_df.drop('sample name', axis = 1, inplace = True)
            
                group_df['sample name'] = group_sample_lst
            
                group_df = group_df[['sample name', 'condition']]
            
                st.sidebar.write(group_df)
                
            elif group_sample_txt == '':
                
                st.sidebar.write('')
                
            else:
                
                st.sidebar.error('Modify your input to appropriately change the order of the samples!')
            
        return group_df
    
    
    
    
    
    # function that updates the lipidXplorer sample names, so, they are consistent with LipidSearch
    def update_sample_name(group_df):  
    
        """
        For example, if you have 4 samples with the following names: WT1, WT2, BQC1, BQC2, the updated names are s1, s2, s3 and s4.  
        
        """
        
        name_dict = {"old name" : [] , "updated name" : []}
        
        counter = 0 
        
        for name in group_df['sample name'].values.tolist():
            
            name_dict['old name'].append(name)
            
            name_dict['updated name'].append('s'+str(counter+1))
            
            counter = counter + 1
                
        name_df = pd.DataFrame.from_dict(name_dict)
        
        return name_df
    
    
    
    
    
     # function to clean/apply filters to the LipidXplorer data
    def apply_filter_lipid_xplorer(df, intsta_df, group_df, rep_lst, cond_lst, filter_mode, filtered_conds, passing_abundance_grade):
    
        # initial cleaning 
        
        total_reps = sum(rep_lst) # total number of all replicates
        
        X = df[['SPECIES', 'MOLSPECIES', 'FORMULA', 'MASS', 'ERROR', 'CLASS']+
                ['PRECURSORINTENSITY:' + group_df['sample name'].iloc[i] + '.mzML' for i in range(total_reps)]]
                
        intsta_df = intsta_df[['SPECIES', 'MOLSPECIES', 'FORMULA', 'MASS', 'ERROR', 'CLASS']+
                ['PRECURSORINTENSITY:' + group_df['sample name'].iloc[i] + '.mzML' for i in range(total_reps)]]
        
        X.rename(columns={"CLASS": "Class"}, inplace=True)
        
        counter = 0
        
        for column in X.columns: # enforces the LipidSearch sample naming convention on LipidXplorer samples
            
            if 'PRECURSORINTENSITY' in column:
                
                X.rename(columns={column: 'MainArea[s'+str(counter+1)+']'}, inplace=True)
                
                intsta_df.rename(columns={column: 'MainArea[s'+str(counter+1)+']'}, inplace=True)
                
                counter = counter + 1
                
        # cleaning X
        
        X = X[X.SPECIES != '###'] # removes the rows that start with ### (comments)
        
        X.dropna(how='all', inplace=True) # removes empty rows  
        
        X.drop_duplicates(keep = 'first', inplace = True) # removes duplicate datapoints 
        
        X.reset_index(inplace=True)
        
        #first filter
        
        X['ERROR'] = pd.to_numeric(X['ERROR'], downcast="float")
        
        X = X.loc[abs(X['ERROR']) < 5] # removes the datapoint if abs(error) > 5
        
        #second filter: minimum abundance grade 
        
        X = build_filter(X, rep_lst, cond_lst, filter_mode, filtered_conds, passing_abundance_grade)
        
        X.reset_index(inplace=True)
        
        X.drop(['level_0', 'ERROR', 'MOLSPECIES'], axis=1, inplace=True) # drops an irrelevant column
        
        X.rename(columns={"index": "old_index"}, inplace=True)
        
        # cleaning internal standards df
        
        intsta_df.dropna(how='all', inplace=True) # removes empty rows
        
        intsta_df.reset_index(inplace=True)
        
        intsta_df.drop(['MOLSPECIES', 'ERROR'], axis=1, inplace=True) # drops an irrelevant column 
        
        intsta_df.rename(columns={"index": "old_index"}, inplace=True)
        
        # removing bad samples 
            
        rep_lst, cond_lst, full_sample_lst, X, r_sample = remove_bad_samples(cond_lst, rep_lst, X)
        
        intsta_df = update_intsta_df(intsta_df, r_sample)
            
        # imputing missing values  
        
        X_impute = add_impute_column(X, full_sample_lst, rep_lst, cond_lst)
            
        X, X_impute, impute_ans = impute_missing_value(X, X_impute, full_sample_lst)
        
        view_clean_data = st.checkbox('View the cleaned data')
        
        if view_clean_data:
        
            st.write('View the cleaned data in conventionl format:')
                
            st.write(X)
                
            csv_downloader(X, 'cleaned_data')
        
            st.write('-----------------------------------------------------------------------------')
        
            st.write('View the cleaned data in log-transformed format:')
                
            log_X = log_transform_df(X, rep_lst, full_sample_lst)
            
            st.write(log_X)
            
            csv_downloader(log_X, 'log_transformed_cleaned_data')
        
            st.write('-----------------------------------------------------------------------------')
        
            st.write('View the internal standards data in conventional format: ')
                
            st.write(intsta_df)
                
            csv_downloader(intsta_df, 'internal_standards')
        
            st.write('-----------------------------------------------------------------------------')
        
            st.write('View the internal standards data in log-transformed format:')
                
            log_intsta_df = log_transform_df(intsta_df, rep_lst, full_sample_lst)
                
            st.write(log_intsta_df)
                
            csv_downloader(log_intsta_df, 'log_transformed_internal_standards')
                
        return X, X_impute, intsta_df, rep_lst, cond_lst, full_sample_lst, impute_ans
    
    
    
    
    
    # function to extract the internal standards dataframe from the uploaded dataset 
    def extract_internal_standards(df):  
    
        mol_lst = df['MOLSPECIES'].values.tolist()
        
        intsta_lst = []
        
        for ele in mol_lst:
            
            if 'splash' in str(ele):
                
                index = mol_lst.index(ele)
                
                intsta_lst.append(index+1)
                
        intsta_df = df.iloc[intsta_lst, :]
        
        df.drop(intsta_lst,0,inplace=True) # removes the internal standard rows from the dataset
        
        intsta_df.dropna(how='all', inplace=True) # removes empty rows
            
        intsta_df['CLASS'] = intsta_df['SPECIES'].apply(lambda x: x.split(' ')[0])
            
        return df, intsta_df
    
    
    
    
    # function to download data
    def csv_downloader(data, name): 
        
        csvfile = data.to_csv()
        
        b64 = base64.b64encode(csvfile.encode()).decode()
        
        filename = name + ".csv"
        
        href = f'<a href="data:file/csv;base64,{b64}" download="{filename}">Download Data</a>'
        
        st.markdown(href,unsafe_allow_html=True)
        
        return
    
    
    
    
    
    # calculate fold change 
    def fc_calculator(num_1, num_2): 
        
        non_zero_num_1 = [num for num in num_1 if num > 0]
        
        non_zero_num_2 = [num for num in num_2 if num > 0]
        
        if len(non_zero_num_1) > 0 and len(non_zero_num_2) > 0:
    
            fc = np.log2(np.mean(non_zero_num_1)/np.mean(non_zero_num_2))
            
        else:
            
            fc = None
        
        return fc
    
    
    
    
    # calculate p-value by running a T-test 
    def p_val_calculator(num_1, num_2):
        
        non_zero_num_1 = [num for num in num_1 if num > 0]
        
        non_zero_num_2 = [num for num in num_2 if num > 0]
        
        if len(non_zero_num_1) > 0 and len(non_zero_num_2) > 0:
        
            t_value, p_value = stats.ttest_ind(non_zero_num_1,non_zero_num_2)
            
        else: 
            
            p_value = None
        
        return p_value
    
    
    
    
    # function that creates the inputs of volcano plots and builds them 
    def volcano_plot(dataset_type, temp, cond_lst, rep_lst, full_sample_lst, impute_ans):
        
        show_vol = st.checkbox("View Volcano Plots")
        
        if show_vol:
            
            if dataset_type == 'LipidXplorer':
                
                for sample in full_sample_lst:
                    
                    temp.rename(columns={'normalized_AUC['+sample+']': 'MainArea['+sample+']'}, inplace=True)
            
            rep_lst_agg = build_rep_lst_agg(rep_lst)
            
            cond_1 = st.selectbox('Pick the first condition', cond_lst, 0)
            
            cond_2 = st.selectbox('Pick the second condition', cond_lst, 1)
            
            for rep in [rep_lst[cond_lst.index(cond_1)], rep_lst[cond_lst.index(cond_2)]]:
                
                if rep < 2:
                    
                    st.error('At least two samples per condition are required for a meanigful analysis!')
                    
                    return None
                
            sample_lst_1 = update_sample_lst(cond_1, cond_lst, rep_lst_agg, full_sample_lst) # updated sample lst corresponding to cond_1
            
            sample_lst_2 = update_sample_lst(cond_2, cond_lst, rep_lst_agg, full_sample_lst) # updated sample lst corresponding to cond_2
            
            temp['p_val_' + cond_1 + '_' + cond_2] = \
                    temp[['MainArea[' + sample + ']'  for sample in sample_lst_1 + sample_lst_2]]\
                    .apply(lambda x: p_val_calculator(x[['MainArea[' + sample + ']'  for sample in sample_lst_1]], \
                                                      x[['MainArea[' + sample + ']'  for sample in sample_lst_2]]), axis=1)
                        
            temp['fc_' + cond_1 + '_' + cond_2] = \
                    temp[['MainArea[' + sample + ']'  for sample in sample_lst_1 + sample_lst_2]]\
                    .apply(lambda x: fc_calculator(x[['MainArea[' + sample + ']'  for sample in sample_lst_1]], \
                                                   x[['MainArea[' + sample + ']'  for sample in sample_lst_2]]), axis=1)
                
            lipid_class_lst = X['Class'].value_counts().index.tolist()
                    
            selected_class = st.multiselect('Add or remove classes (up to 20 classes):', lipid_class_lst, lipid_class_lst[:1])
            
            if len(selected_class) > 20:
                        
                st.error('You can only compare up to 20 lipid classes at a time!')
                        
                return None
                        
            plot, vol_df = volcano_hover(temp, selected_class, cond_1, cond_2, impute_ans)
            
            st.bokeh_chart(plot)  
            
            csv_downloader(vol_df, 'volcano plot')
                
        return
    
    
    
    
    # function useful for creating appropriate legend for volcano plots 
    def create_class_impute_cond_col(x):
        
        if x[-1] == 0:
        
            return x[0] + ' - includes imputation' 
        
        else:
            
            return x[0] + ' - no imputation '
    
    
    
    
    
    # function that prepares a volcano plot with hover tool 
    def volcano_hover(dataframe, selected_class, cond_1, cond_2, impute_ans):
        
        unique_color_lst =[ 'red', 'blue', 'green', 'magenta', 'cyan', 'orange', 'black', 'pink', 'brown', 'yellow', 'purple', \
                            'gray', 'olive', 'chocolate', 'silver', 'darkred', 'khaki', 'skyblue', 'navy', 'orchid']
            
        dataframe['Class_impute_' + cond_1 + cond_2] = dataframe[['Class', cond_1 + '_' + cond_2 + '_impute']]\
            .apply(lambda x: create_class_impute_cond_col(x), axis = 1)
        
        fc = []
        
        pval = []
        
        class_lst = []
        
        color_lst = []
        
        marker_lst = []
        
        species = []
        
        index_lst = []
        
        for lipid_class in selected_class:
            
            fc = fc + dataframe[dataframe['Class'] == lipid_class]['fc_' + cond_1 + '_' + cond_2].values.tolist()
            
            pval = pval + dataframe[dataframe['Class'] == lipid_class]['p_val_' + cond_1 + '_' + cond_2].values.tolist()
            
            color_lst = color_lst + [unique_color_lst[selected_class.index(lipid_class)] for i in range(len(dataframe[dataframe['Class'] == lipid_class]))]
            
            if impute_ans == 'No imputing':
                
                marker_lst = marker_lst + ['circle' for i in range(len(dataframe[dataframe['Class'] == lipid_class]))]
                
                class_lst = class_lst + dataframe[dataframe['Class'] == lipid_class]['Class'].values.tolist()
                
            else:
            
                marker_lst = marker_lst + \
                            ['circle' if ele == 1 else 'triangle' for ele in dataframe[dataframe['Class'] == lipid_class][cond_1 + '_' + cond_2 + '_impute']]
                            
                class_lst = class_lst + dataframe[dataframe['Class'] == lipid_class]['Class_impute_' + cond_1 + cond_2].values.tolist()
                
            index_lst = index_lst + dataframe[dataframe['Class'] == lipid_class].index.values.tolist()
        
            if dataset_type == 'LipidXplorer':
    
                species = species + dataframe[dataframe['Class'] == lipid_class]['SPECIES'].values.tolist()
        
            else:
        
                species = species + dataframe[dataframe['Class'] == lipid_class]['LipidMolec'].values.tolist()
            
        plot = figure(title='Volcano Plot', x_axis_label='Fold Change (' + cond_1 + '/' + cond_2 + ')', y_axis_label='q-value')
            
        vol_df = pd.DataFrame({"FC": fc, "qvalue": -np.log10(pval), "Species": species, "Class": class_lst, \
                               "Color": color_lst, "Marker": marker_lst, "Index": index_lst})
        
        src = ColumnDataSource(vol_df)
        
        plot.scatter(x="FC", y="qvalue", legend_group='Class', color='Color', marker = 'Marker', name='volcano', size = 4, source=src)
        
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
        
        
        
        
        
    # updates the datafraame and rep, cond, full_sample lists after removing a sample
    def update_df_cond_rep_sample_lst(cond_lst, rep_lst, full_sample_lst, r_sample, df): 
        
        '''
        For example, if we have two conditions, A and B, with 3 replicates each, cond lst = [A, B]
        and repl_lst = [3, 3]. The dataframe will have 6 columns and full_sample_lst = [s1, s2, s3, s4, s5, s6].
        If we remove let's say s5, the following function will update the lists: 
        cond_lst = [A, B], rep_lst = [3, 2], full_sample_lst = [s1, s2, s3, s4, s6]. The dataframe will have 5 columns. 
        
        '''
        
        df.drop(['MainArea[' + sample + ']' for sample in r_sample], axis=1, inplace=True) # updating df
        
        for sample in r_sample: # r_sample is the list of removed samples 
        
            full_sample_lst.remove(sample) # updating full_sample_lst
        
        r_sample = [int(item[1:]) for item in r_sample] # separates string and number (e.g. extracts 11 from s11)
        
        rep_lst_agg = build_rep_lst_agg(rep_lst)
        
        used_sample = [] # useful list to avoid overwriting a sample multiple times 
        
        r_sample_cond_index = [] # lst including the index of the condition that corresponds to the removed sample 
        
        for sample_number in r_sample:
            
            for agg_rep in rep_lst_agg:
                
                if (sample_number <= agg_rep) & (sample_number not in used_sample):
                    
                    used_sample.append(sample_number)
                    
                    r_sample_cond_index.append(rep_lst_agg.index(agg_rep))
                    
        for index in r_sample_cond_index:
            
            rep_lst[index] = rep_lst[index] - 1 # updating rep_lst
            
        cond_lst = [cond for cond in cond_lst if rep_lst[cond_lst.index(cond)] != 0] # removing cond with zero corresponding reps 
            
        rep_lst = [rep for rep in rep_lst if rep !=0] # remoing 0's from rep_lst
        
        return rep_lst, cond_lst, full_sample_lst, df
    
    
    
    
    # function that updates sample_lst corresponding to each condition after the full_sample_lst is updated 
    def update_sample_lst(cond, cond_lst, rep_lst_agg, full_sample_lst):
        
        index = cond_lst.index(cond)
        
        if index == 0:
            
            sample_lst = full_sample_lst[0 : rep_lst_agg[0]]
                
        else:
                    
            sample_lst = full_sample_lst[rep_lst_agg[index-1] : rep_lst_agg[index]]
        
        return sample_lst
    
    
    
    
    # function to normalize AUC data
    def normalize_auc_lipid_xplorer(X, X_impute, intsta_df, rep_lst, full_sample_lst):
        
        # extracting the list of classes for which no internal standard is found 
        
        temp = X_impute.copy()
        
        all_class_lst = temp['Class'].unique().tolist()
        
        intsta_class_lst = intsta_df['CLASS'].unique().tolist()
        
        no_intsta_class_lst = [x for x in all_class_lst if x not in intsta_class_lst]
        
        if len(no_intsta_class_lst) > 0:
        
            st.warning("No internal standards is found for the following classes:")
        
            for lipid_class in no_intsta_class_lst:
            
                st.write(lipid_class)
            
        # letting the user pick IS for classes with no IS
        
        add_intsta_class_lst = build_add_intsta_lst(no_intsta_class_lst, intsta_class_lst)
        
        # the amount of each IS in micro mole 
        
        intsta_mic_mol_lst = build_intsta_mic_mol_lst(intsta_class_lst)
            
        for i in range(sum(rep_lst)):
                
            temp['normalized_AUC[s'+str(i+1)+']'] = 0
                
        # normalizing the data 
            
        for lipid_class in all_class_lst:
                
            if lipid_class in intsta_class_lst:
                    
                class_index = intsta_class_lst.index(lipid_class)
                    
                mic_mol = intsta_mic_mol_lst[class_index]
                
                species_index_lst = temp[temp['Class'] == lipid_class].index.tolist()
                    
                intsta_row = find_intsta_index(lipid_class, intsta_class_lst, no_intsta_class_lst, add_intsta_class_lst, intsta_df)
                    
                for i in range(sum(rep_lst)):
                        
                    nums1 = temp['MainArea[s'+str(i+1)+']'].iloc[species_index_lst]
                        
                    num2 = intsta_df['MainArea[s'+str(i+1)+']'].iloc[intsta_row]
                        
                    norm_nums = normalize_auc(nums1, num2, mic_mol)
                        
                    temp['normalized_AUC[s'+str(i+1)+']'].iloc[species_index_lst] = norm_nums
                        
            elif lipid_class in no_intsta_class_lst:
                    
                class_index = no_intsta_class_lst.index(lipid_class)
            
                alt_lipid_class = add_intsta_class_lst[class_index]
                    
                alt_class_index = intsta_class_lst.index(alt_lipid_class)
                    
                mic_mol = intsta_mic_mol_lst[alt_class_index]
                    
                species_index_lst = temp[temp['Class'] == lipid_class].index.tolist()
                    
                intsta_row = find_intsta_index(lipid_class, intsta_class_lst, no_intsta_class_lst, add_intsta_class_lst, intsta_df)
                    
                for i in range(sum(rep_lst)):
                        
                    nums1 = temp['MainArea[s'+str(i+1)+']'].iloc[species_index_lst]
                        
                    num2 = intsta_df['MainArea[s'+str(i+1)+']'].iloc[intsta_row]
                        
                    norm_nums = normalize_auc(nums1, num2, mic_mol)
                    
                    temp['normalized_AUC[s'+str(i+1)+']'].iloc[species_index_lst] = norm_nums
                        
        temp.drop(['MainArea[s' + str(i+1) + ']' for i in range(sum(rep_lst))], axis=1, inplace=True)
        
        temp_impute = temp
        
        temp = temp[['old_index','SPECIES', 'FORMULA', 'MASS', 'Class']+
                ['normalized_AUC[' + sample + ']' for sample in full_sample_lst]]
        
        confirm = st.checkbox('Confirm the inputs & view the normalized data')
        
        if confirm:
                        
            st.write(temp)
        
        return temp, temp_impute
    
    
    
    
    # function useful for normalization of data 
    def normalize_auc(nums1, num2, mic_mol):
        
        norm_nums = []
    
        for num in nums1:
        
            norm_nums.append(num/num2*mic_mol)
        
        return norm_nums
    
    
    
    
    # function for letting the user input the amount of IS 
    def build_intsta_mic_mol_lst(intsta_class_lst):
        
        intsta_mic_mol_lst = []
        
        for lipid_class in intsta_class_lst:
            
            intsta_mic_mol_lst.append(build_mic_mol_input_box(lipid_class))
        
        return intsta_mic_mol_lst
    
    
    
    
    # function for letting the user input the amount of IS
    def build_mic_mol_input_box(lipid_class):
        
        mic_mol = st.number_input('Enter the amount of ' + lipid_class + ' internal standard species in micro mole', 
                                  min_value = 0, max_value = 100000, value = 1, step = 1)
        
        return mic_mol
    
    
    
    
    
    # function that receives the user input on which internal standards to use for classes with no internal standards 
    def build_add_intsta_lst(no_intsta_class_lst, intsta_class_lst): 
            
        add_intsta_lst = []
            
        for lipid_class in no_intsta_class_lst:
                
            add_intsta_lst.append(build_intsta_selectbox(lipid_class, intsta_class_lst)) # function defined below 
                
        return add_intsta_lst
    
    
    
    
            
    # function that creates a select box including the existing list of internal standards 
    def build_intsta_selectbox(lipid_class, intsta_class_lst):  
        
        added_intsta_class = st.selectbox('Pick a internal standard for ' + lipid_class + ' species', intsta_class_lst)
                
        return added_intsta_class
    
    
    
    
    # function for finding the index of IS in intsta_df
    def find_intsta_index(lipid_class, intsta_class_lst, no_intsta_class_lst, add_intsta_class_lst, intsta_df):
        
        if lipid_class in intsta_class_lst:
            
            intsta_row = intsta_df[intsta_df['CLASS'] == lipid_class].index.tolist()[0]
            
        elif lipid_class in no_intsta_class_lst:
            
            class_index = no_intsta_class_lst.index(lipid_class)
            
            alt_lipid_class = add_intsta_class_lst[class_index]
            
            intsta_row = intsta_df[intsta_df['CLASS'] == alt_lipid_class].index.tolist()[0]
        
        return intsta_row
    
    
    
    
    
    # function to remove bad samples and update dataframes and lists after removing bad samples 
    def remove_bad_samples(cond_lst, rep_lst, dataframe):
        
        full_sample_lst = ['s'+str(i+1) for i in range(sum(rep_lst))]
            
        remove_ans = st.radio("Would you like to remove any samples from the analysis?", ['Yes', 'No'], 1)
                
        if remove_ans == 'Yes':
                    
            r_sample = st.multiselect('Pick the sample(s) that you want to remove from the analysis', full_sample_lst)
                    
            rep_lst, cond_lst, full_sample_lst, dataframe = update_df_cond_rep_sample_lst(cond_lst, rep_lst, full_sample_lst, r_sample, dataframe)
            
        else:
            
            r_sample = []
        
        return rep_lst, cond_lst, full_sample_lst, dataframe, r_sample
    
    
    
    
    
    # function to update intsta_df after removing bad samples 
    def update_intsta_df(dataframe, r_sample):
            
        dataframe.drop(['MainArea[' + sample + ']' for sample in r_sample], axis=1, inplace=True) # updating dataframe
        
        return dataframe
    
    
    
    
    # function to calculate the abundance of sfa, mufa and pufa for each lipid species 
    def calculate_SFA_MUFA_PUFA(dataframe, lipid_class, cond, sample_lst):
        
        dataframe = dataframe[dataframe['Class'] == lipid_class]
        
        dataframe[cond + '_mean_AUC'] = dataframe[['MainArea[' + sample + ']' for sample in sample_lst]]\
            .apply(lambda x: np.mean(x), axis = 1)
            
        dataframe[cond + '_var_AUC'] = dataframe[['MainArea[' + sample + ']' for sample in sample_lst]]\
            .apply(lambda x: np.var(x), axis = 1)
            
        dataframe = dataframe[['LipidMolec', cond + '_mean_AUC', cond + '_var_AUC']]
        
        dataframe['FA_ratio'] = dataframe['LipidMolec'].apply(lambda x: calculate_FA_ratio(x))
        
        dataframe['SFA_ratio'] = dataframe['FA_ratio'].apply(lambda x: x[0])
        
        dataframe['MUFA_ratio'] = dataframe['FA_ratio'].apply(lambda x: x[1])
        
        dataframe['PUFA_ratio'] = dataframe['FA_ratio'].apply(lambda x: x[2])
        
        dataframe['SFA_AUC'] = dataframe[[cond + '_mean_AUC', 'SFA_ratio']].apply(lambda x: x[0]*x[1], axis = 1)
        
        dataframe['MUFA_AUC'] = dataframe[[cond + '_mean_AUC', 'MUFA_ratio']].apply(lambda x: x[0]*x[1], axis = 1)
        
        dataframe['PUFA_AUC'] = dataframe[[cond + '_mean_AUC', 'PUFA_ratio']].apply(lambda x: x[0]*x[1], axis = 1)
        
        dataframe['SFA_var'] = dataframe[[cond + '_var_AUC', 'SFA_ratio']].apply(lambda x: x[0]*x[1]**2, axis = 1)
        
        dataframe['MUFA_var'] = dataframe[[cond + '_var_AUC', 'MUFA_ratio']].apply(lambda x: x[0]*x[1]**2, axis = 1)
        
        dataframe['PUFA_var'] = dataframe[[cond + '_var_AUC', 'PUFA_ratio']].apply(lambda x: x[0]*x[1]**2, axis = 1)
        
        sfa = dataframe['SFA_AUC'].sum()
        
        sfa_var = dataframe['SFA_var'].sum()
        
        mufa = dataframe['MUFA_AUC'].sum()
        
        mufa_var = dataframe['MUFA_var'].sum()
        
        pufa = dataframe['PUFA_AUC'].sum()
        
        pufa_var = dataframe['PUFA_var'].sum()
        
        return sfa, mufa, pufa, sfa_var, mufa_var, pufa_var
    
    
    
    
    # function to calculate the ratio of sfa, mufa and pufa for each lipid species 
    def calculate_FA_ratio(mol_structure):
        
        a = mol_structure.split('(')
        
        b = a[1][:-1]
        
        c = b.split('/')
        
        sfa_ratio = 0
        
        mufa_ratio = 0
        
        pufa_ratio = 0
        
        for item in c:
            
            d = item.split(':')[-1]
            
            if d.isnumeric():
                
                d = int(d)
            
                if d == 0:
                
                    sfa_ratio += 1 
                
                elif d == 1:
                
                    mufa_ratio += 1
                
                else:
                
                    pufa_ratio += 1
                    
            else:
                
                if '0' in d:
                    
                    sfa_ratio += 1
                    
                elif '1' in d:
                    
                    mufa_ratio += 1
                    
                else:
                    
                    pufa_ratio += 1
                    
        total = sfa_ratio + mufa_ratio + pufa_ratio
                
        sfa_ratio = sfa_ratio / total
            
        mufa_ratio = mufa_ratio / total
            
        pufa_ratio = pufa_ratio / total
        
        return sfa_ratio, mufa_ratio, pufa_ratio
    
    
    
    
    # function to perform saturation level analysis 
    def saturation_level_plot(dataset_type, temp, full_sample_lst, cond_lst, rep_lst):
            
        if dataset_type == 'LipidXplorer':
                
            for sample in full_sample_lst:
                    
                temp.rename(columns={'normalized_AUC['+sample+']': 'MainArea['+sample+']'}, inplace=True)
                    
        rep_lst_agg = build_rep_lst_agg(rep_lst)
        
        show_typical_bar = st.checkbox("View Typical Bar Plots")
        
        if show_typical_bar:
            
            for lipid_class in temp['Class'].unique():
                
                sfa_lst = []
                
                mufa_lst = []
                
                pufa_lst = []
                
                sfa_var_lst = []
                
                mufa_var_lst = []
                
                pufa_var_lst = []
                
                for cond in cond_lst:
                    
                    sample_lst = update_sample_lst(cond, cond_lst, rep_lst_agg, full_sample_lst)
                    
                    sfa, mufa, pufa, sfa_var, mufa_var, pufa_var = calculate_SFA_MUFA_PUFA(temp, lipid_class, cond, sample_lst)
                    
                    sfa_lst.append(sfa)
                    
                    mufa_lst.append(mufa)
                    
                    pufa_lst.append(pufa)
                    
                    sfa_var_lst.append(np.sqrt(sfa_var))
                    
                    mufa_var_lst.append(np.sqrt(mufa_var))
                    
                    pufa_var_lst.append(np.sqrt(pufa_var))
                    
                data = {'conditions' : cond_lst,
                        'SFA'   : sfa_lst,
                        'SFA_upper' : [x+e for x,e in zip(sfa_lst, sfa_var_lst)],
                        'SFA_lower' : [x-e for x,e in zip(sfa_lst, sfa_var_lst)],
                        'MUFA'   : mufa_lst,
                        'MUFA_upper' : [x+e for x,e in zip(mufa_lst, mufa_var_lst)],
                        'MUFA_lower' : [x-e for x,e in zip(mufa_lst, mufa_var_lst)],
                        'PUFA'   : pufa_lst,
                        'PUFA_upper' : [x+e for x,e in zip(pufa_lst, pufa_var_lst)],
                        'PUFA_lower' : [x-e for x,e in zip(pufa_lst, pufa_var_lst)]}
                
                max_height = max([max(sfa_lst), max(mufa_lst), max(pufa_lst)])
                    
                source = ColumnDataSource(data=data)
                    
                p = figure(x_range = cond_lst, y_range = (0, max_height * 2), height=250, title="Saturation Level Plot - " + lipid_class,
                           x_axis_label= 'Conditions', y_axis_label= 'Total AUC',toolbar_location='right')

                p.vbar(x=dodge('conditions', -0.25, range=p.x_range), top='SFA', width=0.2, source=source,
                           color="#c9d9d3", legend_label="SFA")
                
                p.add_layout(Whisker(source=source, base=dodge('conditions', -0.25, range=p.x_range), \
                                     upper="SFA_upper", lower="SFA_lower", level='overlay'))

                p.vbar(x=dodge('conditions',  0.0,  range=p.x_range), top='MUFA', width=0.2, source=source,
                           color="#718dbf", legend_label="MUFA")
                
                p.add_layout(Whisker(source=source, base=dodge('conditions', 0.0, range=p.x_range), \
                                     upper="MUFA_upper", lower="MUFA_lower", level='overlay'))

                p.vbar(x=dodge('conditions',  0.25, range=p.x_range), top='PUFA', width=0.2, source=source,
                           color="#e84d60", legend_label="PUFA")
                
                p.add_layout(Whisker(source=source, base=dodge('conditions', 0.25, range=p.x_range), \
                                     upper="PUFA_upper", lower="PUFA_lower", level='overlay'))
                    
                p.x_range.range_padding = 0.1

                p.xgrid.grid_line_color = None

                p.legend.location = "top_center"

                p.legend.orientation = "horizontal"
                
                p.legend.label_text_font_size = "10pt"
                
                p.title.text_font_size = "15pt"
        
                p.xaxis.axis_label_text_font_size = "15pt"
            
                p.yaxis.axis_label_text_font_size = "15pt"
            
                p.xaxis.major_label_text_font_size = "15pt"
            
                p.yaxis.major_label_text_font_size = "15pt"
                
                p.yaxis.formatter = BasicTickFormatter(precision=1)
                    
                st.bokeh_chart(p)
                
                sat_df = pd.DataFrame({"Conditions": cond_lst, "SFA": sfa_lst, "SFA_STDV": np.sqrt(sfa_var_lst), "MUFA": mufa_lst, \
                               "MUFA_STDV": np.sqrt(mufa_var_lst), "PUFA": pufa_lst, "PUFA_STDV": np.sqrt(pufa_lst)})
                
                csv_downloader(sat_df, 'Saturation Level plot')
                
        show_stacked_bar_plots = st.checkbox("View Stacked Bar Plots")
        
        if show_stacked_bar_plots:
            
            for lipid_class in temp['Class'].unique():
                
                sfa_lst = []
                
                mufa_lst = []
                
                pufa_lst = []
                
                for cond in cond_lst:
                    
                    sample_lst = update_sample_lst(cond, cond_lst, rep_lst_agg, full_sample_lst)
                    
                    sfa, mufa, pufa, sfa_var, mufa_var, pufa_var = calculate_SFA_MUFA_PUFA(temp, lipid_class, cond, sample_lst)
                    
                    sfa_lst.append(sfa)
                    
                    mufa_lst.append(mufa)
                    
                    pufa_lst.append(pufa)
                    
                data = {'conditions' : cond_lst,
                        'SFA'   : sfa_lst,
                        'MUFA'   : mufa_lst,
                        'PUFA'   : pufa_lst}
                
                tot_lst = sfa_lst + mufa_lst + pufa_lst
                
                max_height = max(tot_lst)
                
                colors = ["#c9d9d3", "#718dbf", "#e84d60"]
                
                p = figure(x_range = cond_lst, y_range = (0, max_height * 2), height=250, title="Saturation Level Plot - " + lipid_class,
                           x_axis_label= 'Conditions', y_axis_label= 'Total AUC',toolbar_location='right')
                
                p.vbar_stack(['SFA', 'MUFA', 'PUFA'], x='conditions', color = colors, width=0.2, source=data,
                             legend_label=['SFA', 'MUFA', 'PUFA'])
                
                p.x_range.range_padding = 0.1

                p.xgrid.grid_line_color = None

                p.legend.location = "top_center"

                p.legend.orientation = "horizontal"
                
                p.legend.label_text_font_size = "10pt"
                
                p.title.text_font_size = "15pt"
        
                p.xaxis.axis_label_text_font_size = "15pt"
            
                p.yaxis.axis_label_text_font_size = "15pt"
            
                p.xaxis.major_label_text_font_size = "15pt"
            
                p.yaxis.major_label_text_font_size = "15pt"
                
                p.yaxis.formatter = BasicTickFormatter(precision=1)
                
                st.bokeh_chart(p)
                
                sat_df = pd.DataFrame({"Conditions": cond_lst, "SFA": sfa_lst, "MUFA": mufa_lst, \
                                       "PUFA": pufa_lst})
                
                csv_downloader(sat_df, 'Saturation Level plot')
                
        show_stacked_percentage_bar_plots = st.checkbox("View Percentage-Based Stacked Bar Plots")
        
        if show_stacked_percentage_bar_plots:
            
            for lipid_class in temp['Class'].unique():
                
                sfa_lst = []
                
                mufa_lst = []
                
                pufa_lst = []
                
                for cond in cond_lst:
                    
                    sample_lst = update_sample_lst(cond, cond_lst, rep_lst_agg, full_sample_lst)
                    
                    sfa, mufa, pufa, sfa_var, mufa_var, pufa_var = calculate_SFA_MUFA_PUFA(temp, lipid_class, cond, sample_lst)
                    
                    sfa_lst.append(sfa)
                    
                    mufa_lst.append(mufa)
                    
                    pufa_lst.append(pufa)
                    
                percent_sfa_lst = [sfa_lst[i]/(sfa_lst[i]+mufa_lst[i]+pufa_lst[i])*100 for i in range(len(sfa_lst))]
                
                percent_mufa_lst = [mufa_lst[i]/(sfa_lst[i]+mufa_lst[i]+pufa_lst[i])*100 for i in range(len(sfa_lst))]
                
                percent_pufa_lst = [pufa_lst[i]/(sfa_lst[i]+mufa_lst[i]+pufa_lst[i])*100 for i in range(len(sfa_lst))]
                    
                data = {'conditions' : cond_lst,
                        'SFA'   : percent_sfa_lst,
                        'MUFA'   : percent_mufa_lst,
                        'PUFA'   : percent_pufa_lst}
                
                colors = ["#c9d9d3", "#718dbf", "#e84d60"]
                
                p = figure(x_range = cond_lst, y_range = (0, 100), height=250, title="Saturation Level Plot - " + lipid_class,
                           x_axis_label= 'Conditions', y_axis_label= 'AUC',toolbar_location='right')
                
                p.vbar_stack(['SFA', 'MUFA', 'PUFA'], x='conditions', color = colors, width=0.2, source=data,
                             legend_label=['SFA', 'MUFA', 'PUFA'])
                
                p.x_range.range_padding = 0.1

                p.xgrid.grid_line_color = None

                p.legend.location = "center_right"

                p.legend.orientation = "vertical"
                
                p.legend.label_text_font_size = "6.5pt"
                
                p.title.text_font_size = "15pt"
        
                p.xaxis.axis_label_text_font_size = "15pt"
            
                p.yaxis.axis_label_text_font_size = "15pt"
            
                p.xaxis.major_label_text_font_size = "15pt"
            
                p.yaxis.major_label_text_font_size = "15pt"
                
                p.yaxis.formatter = BasicTickFormatter(precision=1)
                
                st.bokeh_chart(p)
                
                sat_df = pd.DataFrame({"Conditions": cond_lst, "SFA": percent_sfa_lst, "MUFA": percent_mufa_lst, \
                                       "PUFA": percent_pufa_lst})
                
                csv_downloader(sat_df, 'Saturation Level plot')
            
        return 
    
    ##########################################################################################################################################
    # the main code of the app 
        
    st.header("Data Analysis & Hypothesis Testing Module")
        
    dataset_type = st.radio('Select the type of your dataset:', ['LipidSearch', 'LipidXplorer'])
    
    if dataset_type == 'LipidSearch':
    
        st.sidebar.subheader("Upload Data")
            
        lipid_search = st.sidebar.file_uploader(label='Upload your LipidSearch dataset', type=['csv', 'txt'])
            
        if lipid_search is not None:
                
                df = build_lipidsearch_df(lipid_search)
                
                if df is not None:
                    
                    # building the side bar 
                
                    confirm_data, missing_ans, name_df, passing_letter_grade, pval_tresh, passing_pval_grade, filter_mode,\
                        filtered_conds, passing_abundance_grade, cond_lst, rep_lst = build_sidebar_lipid_search(df)
            
                    # if user confirms the inputs:
            
                    if confirm_data:
                    
                        st.subheader("Further Process Cleaned Data")
                    
                        expand_cleaned_data = st.beta_expander("Impute Missing Values & Remove Low Quality Samples")
                        
                        with expand_cleaned_data:
                    
                            X, X_impute, rep_lst, cond_lst, full_sample_lst, impute_ans = apply_filter_lipid_search\
                                (df, rep_lst, cond_lst, missing_ans, name_df, filter_mode, filtered_conds, passing_abundance_grade) # cleaned data 
                    
                        st.subheader("Analyze Data")
                        
                        expand_vol_plot = st.beta_expander("Volcano Plots")
                        
                        with expand_vol_plot:
                            
                            st.info('In statistics, a volcano plot is a type of scatter-plot that is used to quickly identify changes in \
                                    large data sets composed of replicate data. It plots significance versus fold-change on the y and x axes, respectively. \
                                    A volcano plot combines a measure of statistical significance from a statistical test (e.g., a p value from a T-test) \
                                    with the magnitude of the change, enabling quick visual identification of those data-points that display large magnitude\
                                    changes that are also statistically significant.')
                                    
                            st.markdown('Below, q-value (i.e. -log10(p-value)) of each lipid species is plotted versus the fold change of that species. \
                                The p-value is computed from a two-sample T-test and the fold change is computed from the following formula:')
                                
                            latext = r'''
                                    
                            $$ 
                            Fold Change = log2(\frac{Mean AUC(Condition 1)}{Mean AUC(Condition 2)})
                            $$  
                            
                            '''
                            st.write(latext)
                            
                            temp = X_impute.copy()
                        
                            volcano_plot(dataset_type, temp, cond_lst, rep_lst, full_sample_lst, impute_ans)
                            
                        expand_sat_plot = st.beta_expander("Saturation Level Plots")
                        
                        with expand_sat_plot:
                            
                            st.info('Saturation level plots show the saturation profile of each lipid class.')
                                    
                            st.info('First, for each lipid species, the ratio of Saturated Fatty Acids (SFA), Mono Unsaturated \
                                      Fatty Acids (MUFA) and Poly Unsaturated Fatty Acids (PUFA) is calculated as following:')
                                    
                            st.write('SFA ratio = total number of saturated fatty acids / total number of fatty acids')
                            
                            st.write('MUFA ratio = total number of mono unsaturated fatty acids / total number of fatty acids')
                            
                            st.write('PUFA ratio = total number of poly unsaturated fatty acids / total number of fatty acids')
                            
                            st.info('Then, for each lipid species, the abundance of SFA, MUFA and PUFA is calculated as following:')
                            
                            st.write('AUC(SFA) = (AUC averaged over replicates).(SFA ratio)')
                            
                            st.write('AUC(MUFA) = (AUC averaged over replicates).(MUFA ratio)')
                            
                            st.write('AUC(PUFA) = (AUC averaged over replicates).(PUFA ratio)')
                            
                            st.info('Finally, total AUC(SFA), AUC(MUFA) and AUC(PUFA) for each lipid class is calculated by taking the sum of \
                                    AUC(SFA), AUC(MUFA) and AUC(PUFA) over all lipid species that belong to that class. ')
                            
                            temp = X.copy()
                            
                            saturation_level_plot(dataset_type, temp, full_sample_lst, cond_lst, rep_lst)
                    
    elif dataset_type == 'LipidXplorer':
        
        st.sidebar.subheader('Select Mode')
        
        mode = st.sidebar.radio('', ['Only positive or negative mode', 'Merge positive and negative modes'])
        
        st.sidebar.subheader("Upload Data")
        
        if mode == 'Only positive or negative mode':
            
            lipid_xplorer = st.sidebar.file_uploader(label='Upload your LipidXplorer dataset', type=['csv'])
            
            lipid_xplorer_pos = None
            
            lipid_xplorer_neg = None
            
            if lipid_xplorer is not None:
            
                df, intsta_df = build_single_lipidxplorer_df(lipid_xplorer)
                
            else: 
                
                df = None
            
        elif mode == 'Merge positive and negative modes':
            
            lipid_xplorer_pos = st.sidebar.file_uploader(label='Upload your LipidXplorer dataset in POSITIVE mode', type=['csv'])
            
            lipid_xplorer_neg = st.sidebar.file_uploader(label='Upload your LipidXplorer dataset in NEGATIVE mode', type=['csv'])
            
            lipid_xplorer = None
            
            if lipid_xplorer_pos is not None:
                
                pos_df, pos_intsta_df = build_single_lipidxplorer_df(lipid_xplorer_pos)
                
            if lipid_xplorer_neg is not None:
                
                neg_df, neg_intsta_df = build_single_lipidxplorer_df(lipid_xplorer_neg)
        
            if (lipid_xplorer_pos is not None) and (lipid_xplorer_neg is not None):
            
                df, intsta_df = build_merged_lipidxplorer_df(pos_df, pos_intsta_df, neg_df, neg_intsta_df)
                
            else:
                
                df = None

        if df is not None:
            
            confirm_data, cond_lst, rep_lst, group_df, filter_mode, filtered_conds, passing_abundance_grade = build_sidebar_lipid_xplorer(df)
            
            if confirm_data:
                    
                st.subheader("Further Process Cleaned Data")
                            
                expand_clean_data = st.beta_expander('Impute Missing Values & Remove Low Quality Samples')
                        
                with expand_clean_data:
                        
                    X, X_impute, intsta_df, rep_lst, cond_lst, full_sample_lst, impute_ans \
                        = apply_filter_lipid_xplorer(df, intsta_df, group_df, rep_lst, cond_lst, filter_mode, filtered_conds, passing_abundance_grade) # cleaned data
                        
                st.subheader('Normalize Data')
                
                expand_norm_data = st.beta_expander('Normalized Data')
                
                with expand_norm_data:
                    
                    st.info('An internal standard is a known concentration of a substance that is present in every sample that is analyzed.\
                            In order to calibrate the variation in the data, LipidCruncher normalizes the data using using internal standards.')
                            
                    st.write('The following formula is used for data normalization:')
                            
                    latext = r'''
                                    
                            $$ 
                            Concentration (analyte) = \frac{AUC(analyte)}{AUC(IS)} \times Concentration(IS) 
                            $$  
                            
                            '''
                    st.write(latext)
                                        
                    norm_X, norm_X_impute = normalize_auc_lipid_xplorer(X, X_impute, intsta_df, rep_lst, full_sample_lst)
                
                st.subheader('Analyze Data')
                        
                expand_vol_plot = st.beta_expander("Volcano Plots")
                        
                with expand_vol_plot:
                    
                    st.info('In statistics, a volcano plot is a type of scatter-plot that is used to quickly identify changes in \
                            large data sets composed of replicate data. It plots significance versus fold-change on the y and x axes, respectively. \
                            A volcano plot combines a measure of statistical significance from a statistical test (e.g., a p value from a T-test) \
                            with the magnitude of the change, enabling quick visual identification of those data-points that display large magnitude\
                            changes that are also statistically significant.')
                            
                    st.markdown('Below, q-value (i.e. -log10(p-value)) of each lipid species is plotted versus the fold change of that species. \
                                The p-value is computed from a two-sample T-test and the fold change is computed from the following formula:')
                                
                    latext = r'''
                                    
                            $$ 
                            Fold Change = log2(\frac{Mean AUC(Condition 1)}{Mean AUC(Condition 2)})
                            $$  
                            
                            '''
                    st.write(latext)
                    
                    data_ans = st.radio('Select which data you would like to do analysis on', ['non-normalized data', 'normalized data'], 0)
            
                    if data_ans == 'non-normalized data':
                
                        temp = X_impute.copy()
                
                    else:
                
                        temp = norm_X_impute.copy()
                        
                    volcano_plot(dataset_type, temp, cond_lst, rep_lst, full_sample_lst, impute_ans)
                
                
                
                
                
                

                
                
                
    
    
        
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
        
