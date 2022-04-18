import csv
import streamlit as st
import pandas as pd
import numpy as np
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, HoverTool
import base64
import matplotlib.pyplot as plt
import seaborn as sns 
from sklearn.preprocessing import scale # Data scaling
from sklearn import decomposition #PCA


def app():
    ##########################################################################################################################################
    # functions used in the main code of the app 
    
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
    def log_transform_df(X, rep_lst):  
    
        temp = X.copy()
        
        auc = ['MainArea[s'+str(i+1)+']' for i in range(sum(rep_lst))]
        
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
            
            st.write('View the cleaned data in the conventional format:')
                
            st.write(X)
                
            csv_downloader(X, 'cleaned_data')
            
            st.write('------------------------------------------------------------------------------------------------')
            
            st.write('View the cleaned data in the log-transformed format:')
                
            log_X = log_transform_df(X, rep_lst)
            
            st.write(log_X)
            
            csv_downloader(log_X, 'log_transformed_cleaned_data')
            
            return X 
    
    
    
    
    
    # function to plot histograms 
    def plot_hist(X, rep_lst, cond_lst):
        
        show_hist = st.checkbox('View the distributions of the Area Under the Curve (AUC)')
        
        if show_hist:
                        
            rep_lst_agg = build_rep_lst_agg(rep_lst)
                        
            temp = X.copy()
            
            cond_ans = st.radio('Would you like to filter by condition?', ['Yes', 'No'], 1)
            
            if cond_ans == 'Yes':
                                
                cond = st.selectbox('select a condition', cond_lst)  
                
                sample_lst = build_sample_lst(cond, cond_lst, rep_lst_agg)
                        
                index = cond_lst.index(cond)
                        
                number_rep = rep_lst[index]
                
            else:
                
                sample_lst = ['s'+str(i+1) for i in range(sum(rep_lst))]
            
                number_rep = sum(rep_lst)
            
            filter_ans = st.radio(' Would you like to filter by lipid class?', ['Yes', 'No'], 1)
            
            if filter_ans == 'Yes':
                
                lipid_class = st.selectbox('Select a lipid class', X['Class'].value_counts().index.tolist())
                
                temp = temp[temp['Class'] == lipid_class]
            
            temp = temp[['MainArea[s'+str(i+1)+']' for i in range(sum(rep_lst))]] # picks the 'MainArea[s1]', ..., 'MainArea[sN]' columns only
            
            temp = temp.mask(temp <= 0).fillna(1) # turning 0's and negative numbers  to 1's so it is possible to log transform
            
            if number_rep < 50:
            
                arrange_hist(number_rep, temp, sample_lst)
                
            else:
                
                st.warning("LipidCruncher has limited available resources. You have more than 50 samples, \
                           the best is to plot no more than 50 plots at a time.")
                           
                st.write('Use the double-slider below to set a range of samples you want to investigate. \
                         Make sure to not to go beyond the 50 samples limit.')
                         
                srange = st.slider("", value=[0, number_rep])
                
                diff = srange[1] - srange[0]
                
                if diff <= 50: 
                    
                    show_range_hist = st.button("Show Plots!")
                    
                    if show_range_hist:
                    
                        arrange_hist(diff, temp[['MainArea[s'+str(i+1)+']' for i in range(srange[0], srange[1])]], sample_lst[srange[0]: srange[1]])
                    
                else: 
                    
                    st.error("You have exceeded the 50 samples limit!")
                        
            return
                    
    
    
    
    
    
    # function to arrange histograms in in subplot style, three histograms per row
    def arrange_hist(number_rep, temp, sample_lst):
        
        for i in range(int(round(number_rep/3, 1))):
                            
                col1, col2, col3 = st.columns(3)
                            
                fig = prep_hist(temp, sample_lst, 3*i)
                            
                col1.pyplot(fig)
                            
                fig = prep_hist(temp, sample_lst, 3*i+1)
                            
                col2.pyplot(fig)
                            
                fig = prep_hist(temp, sample_lst, 3*i+2)
                            
                col3.pyplot(fig)
                            
        if number_rep % 3 == 2:
                            
                col1, col2, col3 = st.columns(3)
                            
                fig = prep_hist(temp, sample_lst, -2)
                            
                col1.pyplot(fig)
                            
                fig = prep_hist(temp, sample_lst, -1)
                            
                col2.pyplot(fig)
                            
        elif number_rep %3 == 1:
                            
                col1, col2, col3 = st.columns(3)
                            
                fig = prep_hist(temp, sample_lst, -1)
                            
                col1.pyplot(fig)
        
        return
    
    
    
    
    
    # function that prepares the histogram plot 
    def prep_hist(temp, sample_lst, index):
                        
        plt.rcParams['font.size'] = '35'
                    
        plt.rcParams['axes.linewidth'] = 3
                
        fig, ax = plt.subplots(figsize=(10, 10))
        
        lst = np.log10(temp['MainArea[' + sample_lst[index] + ']'].values.tolist())
            
        ax.hist(lst, bins = 75, range=(0, 12))
                    
        ax.set_title('Histogram of AUC - '+ sample_lst[index], fontsize=50)
                            
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
            
                thresh = 0.7
            
            elif rep_type == 'Technical replicates':
            
                v_min = 0.75
            
                thresh = 0.9
        
            index = cond_lst.index(cond)
            
            # if there's more than replicate
        
            if rep_lst[index] > 1:
        
                temp = X[['MainArea[s' + str(i+1) + ']' for i in range(sum(rep_lst))]].copy()
        
                temp = temp.apply(lambda x: np.log10(x))
            
                # re-naming the columns from MainArea[s1] to s1 etc
        
                counter = 0
        
                for column in temp.columns.values.tolist():
            
                    temp.rename(columns={column: 's' + str(counter+1)}, inplace=True)
            
                    counter = counter + 1
            
                # selecting the columns that correspond to the selected cond    
            
                rep_lst_agg = build_rep_lst_agg(rep_lst)
        
                if index == 0:
        
                    temp = temp[['s' + str(i+1) for i in range(rep_lst_agg[index])]]
        
                else:
            
                    temp = temp[['s' + str(i+1) for i in range(rep_lst_agg[index-1], rep_lst_agg[index])]]
        
                fig = plt.figure(figsize=(20, 16))
            
                mask = np.triu(np.ones_like(temp.corr(), dtype=np.bool))
    
                heatmap = sns.heatmap(temp.corr(), mask=mask, vmin=v_min, vmax=1, center = thresh, annot=False, cmap='RdBu', square=False, cbar=True)
            
                heatmap.set_title('Triangle Correlation Heatmap - ' + cond, fontdict={'fontsize':30});
            
                st.pyplot(fig)
                
                st.write('---------------------------------------------------------------------------------------------------')
                
                st.write('Find the exact correlation coefficients in the table below:')
        
                st.write(temp.corr())
                
                csv_downloader(temp.corr(), 'Correlation_Matrix_'+str(cond))
            
            else:
            
                st.error('The selected condition must have at least two corresponding replicates!')
            
        return
    
    
    
    
    
    # creates list objects useful for retention time plotting
    def build_per_class_auc_retention_label_lst(X, cond_lst, rep_lst, class_name):  
            
            rep_lst_agg = build_rep_lst_agg(rep_lst)
            
            auc = []
                
            retention = []
                
            rep_label = []
            
            index_lst = []
            
            for cond in cond_lst:
                
                sample_lst = build_sample_lst(cond, cond_lst, rep_lst_agg)
                
                column_lst = ['MainArea' + '[' + sample + ']' for sample in sample_lst] 
                # column list corresponding to each condition 
                
                for column in column_lst:
                    
                    rep_label = rep_label + [sample_lst[column_lst.index(column)] \
                                             for i in range(len(X[X['Class'] == class_name][column].values.tolist()))]
                    
                    auc = auc + X[X['Class'] == class_name][column].values.tolist() 
                    # list of  all AUC data for each combination of lipid class/condition (i.e. including all relevant replicates)
                    
                    retention = retention + X[X['Class'] == class_name]['BaseRt'].values.tolist()
                    # list of  all retention time data for each combination of lipid class/condition (i.e. including all relevant replicates)
                    
                    index_lst = index_lst + X[X['Class'] == class_name].index.values.tolist()
            
            return auc, retention, rep_label, index_lst 
        
    
    
    
    
    # plots a single retention time plot with hover tool
    def retention_hover(X, rep_lst, cond_lst, class_name):  
            
            retention = X[X['Class'] == class_name]['BaseRt'].values.tolist()
            
            mass = X[X['Class'] == class_name]['Calc Mass'].values.tolist()
            
            species = X[X['Class'] == class_name]['LipidMolec'].values.tolist()
            
            retention_df = pd.DataFrame({"Mass": mass, "Retention": retention, "Species": species})
            
            # bokeh plot with hover tool
            src = ColumnDataSource(retention_df)
                
            plot = figure(title=class_name, x_axis_label='Calculated Mass', y_axis_label='Retention Time (mins)')
            
            plot.scatter(x="Mass", y="Retention", source=src)
                
            # hover tool
            hover = HoverTool(tooltips = [('Mass', '@Mass'), ('Retention_time', '@Retention'), ('Species', '@Species')])
            
            plot.add_tools(hover)
            
            plot.title.text_font_size = '15pt'
            
            plot.xaxis.axis_label_text_font_size = "15pt"
            
            plot.yaxis.axis_label_text_font_size = "15pt"
            
            plot.xaxis.major_label_text_font_size = "15pt"
            
            plot.yaxis.major_label_text_font_size = "15pt"
            
            return plot, retention_df
    
    
    
    
    # comparison mode - retention time plot 
    def retention_multi(X, rep_lst, selected_class): 
        
        unique_color_lst =[ 'red', 'blue', 'green', 'magenta', 'cyan', 'orange', 'black', 'pink', 'brown', 'yellow', 'purple', \
                            'gray', 'olive', 'chocolate', 'silver', 'darkred', 'khaki', 'skyblue', 'navy', 'orchid']
        
        retention = []
        
        mass = []
        
        class_lst = []
        
        color_lst = []
        
        for lipid_class in selected_class:
            
            retention = retention + X[X['Class'] == lipid_class]['BaseRt'].values.tolist()
            
            mass = mass + X[X['Class'] == lipid_class]['Calc Mass'].values.tolist()
            
            class_lst = class_lst + X[X['Class'] == lipid_class]['Class'].values.tolist()
            
            color_lst = color_lst + [unique_color_lst[selected_class.index(lipid_class)] for i in range(len(X[X['Class'] == lipid_class]))]
        
        retention_df = pd.DataFrame({"Mass": mass, "Retention": retention, "Class": class_lst, "Color": color_lst})
            
        # bokeh plot with hover tool
        src = ColumnDataSource(retention_df)
        
        plot = figure(title='Retention Time vs. Mass - Comparison Mode', x_axis_label='Calculated Mass', y_axis_label='Retention Time (mins)')
            
        plot.scatter(x="Mass", y="Retention", legend_group='Class', color='Color', source=src)
            
        plot.title.text_font_size = '15pt'
            
        plot.xaxis.axis_label_text_font_size = "15pt"
            
        plot.yaxis.axis_label_text_font_size = "15pt"
            
        plot.xaxis.major_label_text_font_size = "15pt"
            
        plot.yaxis.major_label_text_font_size = "15pt"
           
        return plot, retention_df     
        
    
    
    
    # all retention time plots
    def plot_retention(X, cond_lst, rep_lst):  
            
            show_retention = st.checkbox("View the retention time plots")
            
            if show_retention:
                
                mode = st.radio('', ['Comparison Mode', 'Individual Mode'])
                
                if mode == 'Individual Mode':
                
                    lipid_class_lst = X['Class'].value_counts().index.tolist()
                
                    for lipid_class in lipid_class_lst:
                    
                        class_name = lipid_class
                    
                        plot, retention_df = retention_hover(X, rep_lst, cond_lst, class_name)
                    
                        st.bokeh_chart(plot)
                    
                        file_name = 'retention_plot_' + class_name
                    
                        csv_downloader(retention_df, file_name)
                    
                        st.write("----------------------------------------------------------------------------------------------------")
                        
                else:
                    
                    lipid_class_lst = X['Class'].value_counts().index.tolist()
                    
                    selected_class = st.multiselect('Add or remove classes (up to 20 classes):', lipid_class_lst, lipid_class_lst[:2])
                    
                    if len(selected_class) > 20:
                        
                        st.error('You can only compare up to 20 lipid classes at a time!')
                        
                        return None
                        
                    else:
                    
                        plot, retention_df = retention_multi(X, rep_lst, selected_class)
                    
                        st.bokeh_chart(plot)
                        
                        csv_downloader(retention_df, 'Retention_Time_Comparison_Mode')
                    
            return
        
    
    
    
    
    # PCA math 
    def run_pca(df): # PCA 
    
            # turning 0's and negative numbers to 1's for log-transformation 
            
            df = df.mask(df<=0).fillna(1)
            
            df = df.apply(lambda x: np.log10(x))
            
            # transposing the dataset 
            
            X = df.T
                    
            X = scale(X)
            
            pca = decomposition.PCA(n_components=3)
            
            pca.fit(X)
            
            scores = pca.transform(X)
            
            explained_variance = pca.explained_variance_ratio_
            
            return scores, explained_variance
            
        
    
    
    
    
    # plotting PCA
    def plot_pca(df, rep_lst, cond_lst):
    
            df = df[['MainArea[s'+str(i+1)+']' for i in range(sum(rep_lst))] + ['Class']]
            
            full_sample_lst = ['s'+str(i+1) for i in range(sum(rep_lst))]
            
            show_pca = st.checkbox("Run Principal Component Analysis (PCA)")
            
            if show_pca:
                
                remove_ans = st.radio("Would you like to remove any samples from the analysis?", ['Yes', 'No'], 1)
                
                if remove_ans == 'Yes':
                    
                    r_sample = st.multiselect('Pick the sample(s) that you want to remove from the analysis', full_sample_lst)
                    
                    if (len(full_sample_lst) - len(r_sample)) >= 3:
                    
                        rep_lst, cond_lst, full_sample_lst, df = update_df_cond_rep_sample_lst(cond_lst, rep_lst, full_sample_lst, r_sample, df)
                        
                    else:
                        
                        st.error('At least three samples are required for a meanigful analysis!')
                        
                        return None
                        
                filter_ans = st.radio("Would you like to filter by lipid class?", ['Yes', 'No'], 1)
                
                if filter_ans == 'Yes':
                    
                    lipid_class = st.selectbox('Pick a lipid class', X['Class'].value_counts().index.tolist())
                    
                    df = df[df['Class'] == lipid_class]
            
                df.drop(['Class'], axis=1, inplace=True)
            
                scores, explained_variance = run_pca(df)
                
                PC_lst = ['PC'+str(i+1)+' ('+str("{:.0f}".format(explained_variance[i]*100))+'%)' for i in range(3)]
                
                sel_PCx = st.selectbox('Pick the PC shown on the x-axis', PC_lst, 0)
                
                sel_PCy = st.selectbox('Pick the PC shown on the y-axixs', PC_lst, 1)
                
                x_index = int(sel_PCx.split(' ')[0][2])
                
                y_index = int(sel_PCy.split(' ')[0][2])
                
                PCx = scores[:, x_index - 1].tolist()
                
                PCy = scores[:, y_index - 1].tolist()
            
                rep_lst_agg = build_rep_lst_agg(rep_lst)
            
                color_lst =[ 'red', 'blue', 'green', 'magenta', 'cyan', 'orange', 'black', 'pink', 'brown', 'yellow', 'purple', \
                            'gray', 'olive', 'chocolate', 'silver', 'darkred', 'khaki', 'skyblue', 'navy', 'orchid']
            
            
                plot = figure(title='PCA', x_axis_label=sel_PCx, y_axis_label=sel_PCy)
            
                x = []
            
                y = []
            
                legend = []
            
                color = []
            
                for cond in cond_lst:
                
                    if cond_lst.index(cond) == 0:
                    
                        x = x + PCx[0: rep_lst_agg[cond_lst.index(cond)]]
                    
                        y = y + PCy[0: rep_lst_agg[cond_lst.index(cond)]]
                    
                        legend = legend + [cond for i in range(rep_lst_agg[cond_lst.index(cond)])]
                    
                        color = color + [color_lst[cond_lst.index(cond)] for i in range(rep_lst_agg[cond_lst.index(cond)])]
                        
                    else:
                    
                        x = x + PCx[rep_lst_agg[cond_lst.index(cond)-1]: rep_lst_agg[cond_lst.index(cond)]]
                    
                        y = y + PCy[rep_lst_agg[cond_lst.index(cond)-1]: rep_lst_agg[cond_lst.index(cond)]]
                    
                        legend = legend + [cond for i in range(rep_lst_agg[cond_lst.index(cond)-1], rep_lst_agg[cond_lst.index(cond)])]
                    
                        color = color + [color_lst[cond_lst.index(cond)] for i in range(rep_lst_agg[cond_lst.index(cond)-1], rep_lst_agg[cond_lst.index(cond)])]
            
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
                
                    csv_downloader(pca_df[['PC1', 'PC2', 'sample', 'legend']], 'PCA')
                            
            return
           
    
    
    
    
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
            
                st.sidebar.markdown('The following additional class(es) are detected in the positive mode. \
                                    Use the drop down menu to add additional classes to the list of the positive mode classes.')
            
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
        
        st.write('View the cleaned data in conventionl format:')
                
        st.write(X)
                
        csv_downloader(X, 'cleaned_data')
        
        st.write('-----------------------------------------------------------------------------')
        
        st.write('View the cleaned data in log-transformed format:')
                
        log_X = log_transform_df(X, rep_lst)
            
        st.write(log_X)
            
        csv_downloader(log_X, 'log_transformed_cleaned_data')
        
        st.write('-----------------------------------------------------------------------------')
        
        st.write('View the internal standards data in conventional format: ')
                
        st.write(intsta_df)
                
        csv_downloader(intsta_df, 'internal_standards')
        
        st.write('-----------------------------------------------------------------------------')
        
        st.write('View the internal standards data in log-transformed format:')
                
        log_intsta_df = log_transform_df(intsta_df, rep_lst)
                
        st.write(log_intsta_df)
                
        csv_downloader(log_intsta_df, 'log_transformed_internal_standards')
                
        return X, intsta_df
    
    
    
    
    
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
            
        return df, intsta_df 
    
    
    
    
    
    # function to plot the CoV of lipid species
    def plot_cov(X, cond_lst, rep_lst, dataset_type):  
        
        auc = ['MainArea[s'+str(i+1)+']' for i in range(sum(rep_lst))]
    
        X[auc]=X[auc].mask(X[auc]<=0).fillna(1) # turning 0's and negative numbers  to 1's so it is possible to log transform
        
        show_cov = st.checkbox("Run Coefficient of Variation Analysis (CoV)")
        
        if show_cov:
            
            cond = st.radio('Which of the following conditions has corresponding technical replicates? ', cond_lst+['None of the above'], len(cond_lst))
            
            if cond != 'None of the above':
                
                index = cond_lst.index(cond)
                
                if rep_lst[index] == 1:
                    
                    st.error('The selected condition must have at least two corresponding replicates!')
                    
                else:
            
                    rep_lst_agg = build_rep_lst_agg(rep_lst)
                
                    sample_lst = build_sample_lst(cond, cond_lst, rep_lst_agg)
                
                    X_plot, X_cov_df = cov_hover(X, sample_lst, dataset_type)
                
                    st.bokeh_chart(X_plot)
                
                    csv_downloader(X_cov_df, 'CoV_All_Lipid_Species')
                
        return
    
    
    
    
    
    # one CoV plot with hover tool
    def cov_hover(df, sample_lst, dataset_type):  
            
        df['cov'] = df[['MainArea['+sample+']' for sample in sample_lst]].apply(lambda x: cov_calculator(x), axis=1)
    
        df['mean'] = df[['MainArea['+sample+']' for sample in sample_lst]].apply(lambda x: mean_calculator(x), axis=1)
            
        plot = figure(title='CoV - All lipid Species', x_axis_label='Mean of log10(AUC) Across All BQC Samples', y_axis_label='CoV(%)')
            
        x = df['mean'].values.tolist()
            
        y = df['cov'].values.tolist()
        
        max_y = int(np.max(y))
    
        if dataset_type == 'LipidXplorer':
    
            species = df['SPECIES'].values.tolist()
        
        else:
        
            species = df['LipidMolec'].values.tolist()
            
        cov_df = pd.DataFrame({"Mean_AUC": x, "CoV": y, 'Species': species})
        
        src = ColumnDataSource(cov_df)
        
        plot.scatter(x="Mean_AUC", y="CoV", name='cov', source=src)
    
        plot.line(x=[i for i in range(9)], y=25, color='red')
        
        plot.line(x=3, y=[i for i in range(max_y)], color='green')
        
        hover = HoverTool(tooltips = [('Mean_AUC', '@Mean_AUC'), ('CoV', "@CoV"), ('Species', "@Species")], names=['cov'])
        
        plot.add_tools(hover)
        
        plot.title.text_font_size = "15pt"
        
        plot.xaxis.axis_label_text_font_size = "15pt"
            
        plot.yaxis.axis_label_text_font_size = "15pt"
            
        plot.xaxis.major_label_text_font_size = "15pt"
            
        plot.yaxis.major_label_text_font_size = "15pt"
    
        return plot, cov_df  
    
    
    
    
    
    # calculate CoV
    def cov_calculator(numbers): 
        
        non_zero_lst = [number for number in numbers if (number>0)]
        
        if len(non_zero_lst) > 0:
    
            cov = np.std(non_zero_lst)/np.mean(non_zero_lst)*100
            
        else:
            
            cov = None 
        
        return cov
    
    
    
    
    
    # calculate mean
    def mean_calculator(numbers): 
        
        non_zero_lst = [number for number in numbers if (number>0)]
        
        if len(non_zero_lst) > 0:
    
            mean = np.mean(non_zero_lst)
        
            mean = np.log10(mean)
            
        else:
            
            mean = None
        
        return mean
    
    
    
    
    # function to download data
    def csv_downloader(data, name): 
        
        csvfile = data.to_csv()
        
        b64 = base64.b64encode(csvfile.encode()).decode()
        
        filename = name + ".csv"
        
        href = f'<a href="data:file/csv;base64,{b64}" download="{filename}">Download Data</a>'
        
        st.markdown(href,unsafe_allow_html=True)
        
        return
    
    
    
    
    
    @st.cache
    def convert_df(dataframe):

        return dataframe.to_csv().encode('utf-8')
    
    ##########################################################################################################################################
    # the main code of the app 
        
    st.header("Data Exploration & Quality Check Module")
        
    st.markdown("""
    
                    The quality control module of LipidCruncher provides the user with a visual scan through the data and runs few sanity tests on it.
                    
                    Start by selecting the type of your dataset and completing the next steps on the side bar. 
    
                    """)
        
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
                    
                        st.subheader("Explore Data")
                    
                        expand_raw_data = st.expander("Raw Data")
            
                        with expand_raw_data:
                            
                            st.write("View the raw data:")
                
                            st.write(df)
                    
                            csv_download = convert_df(df)
                            
                            st.download_button(
                                label="Download data",
                                data=csv_download,
                                file_name='data.csv',
                                mime='text/csv')
                            
                            #csv_downloader(df, 'data')
                    
                        expand_cleaned_data = st.expander("Cleaned Data")
                        
                        with expand_cleaned_data:
                            
                            st.info("""
                            
                            The data cleaning process is a three steps process: 
                                
                            1) LipidCruncher deletes the datapoints that do not pass through the filters: 
                                either their associated "Rej" value is 1, or their MainGrade or APvalue score is lower than the set minimum.
                            
                            2) LipidCruncher only keeps the relevant columns: "LipidMolec", "Class", "Calc Mass", "BaseRt", "MainArea[s1]", 
                                ..., "MainArea[sN]". Where N is the total number of the samples.
                            
                            3) LipidCruncher adds a column named "old_index" to the cleaned dataset 
                                which refers to the index of the lipid species in the raw dataset.
                            
                            """)
                    
                            X = apply_filter_lipid_search(df, rep_lst, cond_lst, missing_ans, name_df, filter_mode, filtered_conds, passing_abundance_grade) # cleaned data 
                            
                        expand_hist = st.expander("Distributions of the Area Under the Curve (AUC)")
                        
                        with expand_hist:
                            
                            st.info("""
        
                                In a standard LipidSearch dataset, columns "MainArea[s1]" to "MainArea[sN]" correspond to Area 
                                Under the Curve (AUC) representing the relative abundance of the lipid species 
                                in samples s1 to sN. 
        
                                """)
                        
                            st.info("""
                            
                                To plot the histograms of AUC, LipidCruncher turns all 0 values (i.e. missing values) to 1 
                                and then log-transforms the whole dataset. This allows the user to visualize what portion of the 
                                values are missing without affecting the distribution (as 1 is orders of magnitude smaller than 
                                the minimum detection threshold by mass spectrometry).
                        
                                """)
                                
                            st.write(""" Visualize the distribution of AUC's in any of the replicates and compare with other replicates: """)
                        
                            plot_hist(X, rep_lst, cond_lst) # histograms of AUC
                    
                        expand_retention = st.expander('Retention Time of Lipid Species')
                        
                        with expand_retention:
                            
                            st.info("""
                                
                                The retention time of a lipid species is a function of its degree of hydrophobicity. 
                                The more hydrophobic the lipid species, the longer the retention time. When retention time is 
                                plotted versus molecular mass, lipid species tend to form separate clusters based upon which lipid clas they belong to. 
                                This plot (comparison mode) is a useful tool for detecting misidentified lipid species. 
                                
                                """)
                            
                            st.write('Inspect the retention time of lipid species within any lipid class and compare with other lipid classes:')    
                    
                            plot_retention(X, cond_lst, rep_lst) # retention time plots
                    
                        st.subheader("Run Quality Checks")
                    
                        expand_corr = st.expander('Pairwise Correlation Analysis') 
                        
                        with expand_corr:
                            
                            st.info("""
                                
                                Typically, the AUC's of lipid species in any sample is highly linearly correlated to those in its biological replicate
                                (i.e. correlation coefficient > 0.7). This linear correlation is expected to be even stronger for technical replicates 
                                (i.e. correlation coefficient > 0.9).
                                A weaker correlation should raise suspicion on the quality of the involving replicates.
                                
                                """)
                                
                            st.info("LipidCruncher removes the missing values before preforming the correlation test.")
                            
                            st.write('Run a correlation test to inspect the degree of correlation between any two biological or technical replicates:')
                    
                            plot_corr(X, cond_lst, rep_lst) # pairwise correlation plots 
                    
                        expand_pca = st.expander('Principal Component Analysis (PCA)')
                    
                        with expand_pca:
                        
                            st.info("""
                                
                                Principal Component Analysis, or PCA, is a dimensionality-reduction method that is often used to reduce the 
                                dimensionality of large data sets, by transforming a large set of variables into a smaller one that still contains 
                                most of the information in the large set.
                                
                                """)
                                
                            st.info("""
                                
                                Typically, biological replicates are expected to cluster together with the exception of rare outliers that fall 
                                further away from the rest of the replicates. However, technical replicates are always expected to tightly 
                                cluster together. 
                                
                                """)
                                
                            st.info(""" 
                                
                                A plot of the top two principal components (PC1 and PC2) against each other is the best way to inspect the clustering 
                                of different samples. This is because the PC's are ordered based on how much of the variability in the data 
                                they can explain (i.e. PC1 is the PC with the highest variance explained ratio, PC2 is the PC with the second highest 
                                variance expalined ratio and so on).
                                 
                                
                                """)
                        
                            st.write('Run PCA to inspect the clustering of different samples:')
                    
                            plot_pca(X, rep_lst, cond_lst) # PCA analysis 
                    
                    
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
                    
                st.subheader("Explore Data")
                    
                expand_raw_data = st.expander("Raw Data")
            
                with expand_raw_data:
                            
                    st.write('View the raw data:')
                
                    st.write(df)
                            
                expand_clean_data = st.expander('Cleaned Data')
                        
                with expand_clean_data:
                            
                    st.info("""
                            
                        The data cleaning process is a six steps process: 
                                
                        1) LipidCruncher removes the empty rwos and deletes the duplicated datapoints.  
                                
                        2) LipidCruncher deletes the datapoints that do not pass through the applied filter.
                            
                        3) LipidCruncher only keeps the relevant columns: "SPECIES", "FORMULA", "MASS", "CLASS",
                            "PRECURSORINTENSITY[s1]", ..., "PRECURSORINTENSITY[sN]". Where N is the total number of the samples.
                                
                        4) To keep the naming convention consistent between LipidSearch and LipidXplorer, 
                            LipidCruncher changes the name of the following columns from
                            "PRECURSORINTENSITY[s1]", ..., "PRECURSORINTENSITY[sN]"
                            to
                            "MainArea[s1]", ..., "MainArea[sN]". 
                            
                        5) LipidCruncher adds a column named "old_index" to the cleaned dataset 
                            which refers to the index of the lipid species in the raw dataset.
                                
                        6) LipidCruncher extracts the internal standards lipid species and puts them in a separate dataset.
                            An internal standard is a lipid that is added in a known constant amount to the samples.
                            As the rate of spraying the lipids in shotgun method is not constant, calibration is required. 
                            Internal standard lipids can be used for calibration purposes.
                            
                        """)
                        
                    X, intsta_f = apply_filter_lipid_xplorer(df, intsta_df, group_df, rep_lst, cond_lst, filter_mode, filtered_conds, passing_abundance_grade) # cleaned data
                    
                expand_hist = st.expander('Distributions of the Area Under the Curve (AUC)')
                    
                with expand_hist:
                            
                        st.info("""
        
                            In a cleaned LipidXplorer dataset, columns "MainArea[s1]" to "MainArea[sN]" correspond to Area 
                            Under the Curve (AUC) representing the relative abundance of the lipid species 
                            in samples s1 to sN. 
    
                            """)
                        
                        st.info("""
                            
                            To plot the histograms of AUC, LipidCruncher turns all 0 values (i.e. missing values) to 1 
                            and then log-transforms the whole dataset. This allows the user to visualize what portion of the 
                            values are missing without affecting the distribution (as 1 is orders of magnitude smaller than 
                            the minimum detection threshold by mass spectrometry).
                        
                            """)
                                
                        st.write(""" Visualize the distribution of AUC's in any of the replicates and compare with other replicates: """)
                    
                        plot_hist(X, rep_lst, cond_lst) # histograms of AU
                    
                st.subheader('Run Quality Checks')
                    
                expand_corr = st.expander('Pairwise Correlation Analysis')
                    
                with expand_corr:
                            
                        st.info("""
                                
                            Typically, the AUC's of lipid species in any sample is highly linearly correlated to those in its biological replicate
                            (i.e. correlation coefficient > 0.7). This linear correlation is expected to be even stronger for technical replicates 
                            (i.e. correlation coefficient > 0.9).
                            A weaker correlation should raise suspicion on the quality of the involving replicates.
                                
                            """)
                                
                        st.info("LipidCruncher removes the missing values before preforming the correlation test.")
                            
                        st.write('Run a correlation test to inspect the degree of correlation between any two biological or technical replicates:')
                    
                        plot_corr(X, cond_lst, rep_lst) # pairwise correlation plots
                    
                expand_pca = st.expander('Principal Component Analysis')
                    
                with expand_pca:
                            
                        st.info("""
                                
                            Principal Component Analysis, or PCA, is a dimensionality-reduction method that is often used to reduce the 
                            dimensionality of large data sets, by transforming a large set of variables into a smaller one that still contains 
                            most of the information in the large set.
                                
                            """)
                                
                        st.info("""
                                
                            Typically, biological replicates are expected to cluster together with the exception of rare outliers that fall 
                            further away from the rest of the replicates. However, technical replicates are always expected to tightly 
                            cluster together. 
                            
                            """)
                                
                        st.info(""" 
                                
                            A plot of the top two principal components (PC1 and PC2) against each other is the best way to inspect the clustering 
                            of different samples. This is because the PC's are ordered based on how much of the variability in the data 
                            they can explain (i.e. PC1 is the PC with the highest variance explained ratio, PC2 is the PC with the second highest variance 
                            expalined ratio and so on).
                                 
                                
                            """)
                        
                        st.write('Run PCA to inspect the clustering of different samples:')
                        
                        plot_pca(X, rep_lst, cond_lst) # PCA analysis
                    
                st.subheader('Evaluate the Quality of Technical Replicates')
                    
                expand_cov = st.expander('Coefficient of Variation Analysis')
                    
                with expand_cov:
                            
                        st.info(""" 
                                
                            The coefficient of variation (CoV) is defined as the ratio of the standard deviation to the mean.
                            It shows the extent of variability in relation to the mean of the population and is often expressed as a percentage.
                                
                            """)
                                
                        st.info(""" 
                                
                            Technical replicates are expected to be approximately identical. Therefore, the vast majority of the lipid species
                            with a abundance higher than the "limit of detection" must have a very low CoV (i.e. CoV < 25%).
                            The red line in the CoV plot is the line of "25% CoV" and the green line is the "limit of detection" line.
                                
                            """)
                                
                        st.write('Run a CoV analysis to evaluate the quality of your technical replicates:')
                    
                        plot_cov(X, cond_lst, rep_lst, dataset_type)
                
                
                
                
                
                

                
                
                
    
    
        
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
        
