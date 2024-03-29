import streamlit as st
import pandas as pd
import numpy as np
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, HoverTool, Whisker, BasicTickFormatter
import matplotlib.pyplot as plt
import seaborn as sns 
from sklearn.preprocessing import scale # Data scaling
from sklearn import decomposition # PCA
from scipy import stats
from bokeh.transform import dodge
import math

def app():

##########################################################################################################################################
# functions used in the main code of the app

    def build_lipidsearch_df(lipidsearch):
        
        with open(lipid_search.name, "wb") as f:
            
            f.write(lipid_search.getbuffer()) # saves the file in the app directory
            
        if ".txt" in lipid_search.name:
                
                df = txt_to_df(lipid_search.name)
                
                #df, intsta_df = extract_internal_standards(df)
                
                cols = df.columns.values.tolist()
                
                if ('LipidMolec' in cols) and ('ClassKey' in cols) and ('FAKey' in cols) and ('CalcMass in cols') and ('BaseRt' in cols):
                
                    return df
                
                else:
                    
                    st.sidebar.error('This is not a valid LipidSearch dataset!')
                    
                    return None
                
        elif ".csv" in lipid_search.name:
                
                df = csv_to_df(lipid_search.name)
                
                #df, intsta_df = extract_internal_standards(df)
                
                cols = df.columns.values.tolist()
                
                if ('LipidMolec' in cols) and ('ClassKey' in cols) and ('FAKey' in cols) and ('CalcMass in cols') and ('BaseRt' in cols):
                    
                    return df
                
                else:
                    
                    st.sidebar.error('This is not a valid LipidSearch dataset!')
                    
                    return None
                
                
                
                

# function to extract the internal standards dataframe from the uploaded dataset 
    def extract_internal_standards(df):
    
        mol_lst = df['LipidMolec'].values.tolist()
        
        intsta_lst = []
        
        for ele in mol_lst:
            
            #ele_split = str(ele).split(' ')
            
            if '(s)' in ele:
                
                index = mol_lst.index(ele)
                
                intsta_lst.append(index)
                
        intsta_df = df.iloc[intsta_lst, :]
        
        df.drop(intsta_lst,0,inplace=True) # removes the internal standard rows from the dataset
            
        return df, intsta_df
            
            
            
            
    def txt_to_df(dataset):
        
        return pd.read_csv(lipid_search.name, delimiter = "\t", error_bad_lines=False)
    
    
    
    
    
    @st.cache
    def csv_to_df(dataset):
        
        return pd.read_csv(lipid_search.name, error_bad_lines=False, encoding = "ISO-8859-1")





# function that builds the side bar
    def build_sidebar(df):   
            
        st.sidebar.subheader("Define Experiment")
            
        n_cond = st.sidebar.number_input('Enter the number of conditions',
                                         min_value = 1, max_value= 20, value = 1, step = 1)
            
        cond_lst, rep_lst = build_cond_rep_lst(n_cond)
        
        if check_input_validity(df, rep_lst) == False:
            
            st.sidebar.error('The inputs are incomplete and/or inaccurate!')
            
            return False, None, None, None, None, None
            
        # if the inputs are valid
        
        else:
            
            st.sidebar.subheader('Group Samples')
        
            group_df = group_samples(df, cond_lst, rep_lst)
            
            st.sidebar.subheader('Apply Filter')
            
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
            
        return confirm_data, name_df, filtered_conds, passing_abundance_grade, cond_lst, rep_lst
    
    
    
    
    
    @st.cache
    def check_input_validity(df, rep_lst):
        
        counter = 0 
        
        for column in df.columns:
            
            if ('MeanArea[s' in column) and ('Org' not in column):
                
                counter = counter + 1
                
        if counter == sum(rep_lst):
            
            return True
        
        else:
        
            return False 
    
    
    
    
    
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
            
        main_area_col_lst = build_main_area_col_lst(df)

        group_dict = {'sample name' : [], 'condition' : []}
        
        group_dict['sample name'] = ['s' + str(ele) for ele in main_area_col_lst]
        
        group_dict['condition'] = extensive_cond_lst
            
        group_df = pd.DataFrame.from_dict(group_dict) # df including sample_name and corresponding condition 
        
        return group_df 
    
    
    
    
    
    @st.cache
    def build_main_area_col_lst(df):
        
        main_area_col_lst = []
         
        for col in df.columns:
             
             if ('MeanArea[s' in col) and ('Org' not in col):
                 
                 main_area_col_lst.append(int(col[10:-1]))
        
        main_area_col_lst.sort()
                 
        return main_area_col_lst
    
    
    
    
    
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
        
        
        
        
        
    @st.cache
    def convert_df(dataframe):

        return dataframe.to_csv().encode('utf-8')
    
    
    
    
    
# function to apply filters to the data
    def clean_data(df, rep_lst, cond_lst, name_df, filtered_conds, passing_abundance_grade): 
    
            temp = df[['LipidMolec', 'ClassKey', 'FAKey', 'CalcMass', 'BaseRt', 'TotalGrade'] + ['MeanArea[' + sample + ']' for sample in name_df['old name']]]
            
            #first filter
            
            temp = temp[temp['TotalGrade'].isin(['A', 'B'])]
             
            # updating the column names
            
            total_reps = sum(rep_lst) # total number of all replicates
            
            for (sample_1, sample_2) in zip(name_df['old name'], name_df['updated name']):
                
                temp.rename(columns={'MainArea[' + sample_1 + ']' : 'Area[' + sample_2 + ']'}, inplace = True)
                
            for i in range(total_reps):
                
                temp.rename(columns={'Area[s' + str(i+1) + ']' : 'MainArea[s' + str(i+1) + ']'}, inplace = True)
                
            # creating new df including total AUC for each lipid species 
            
            X1, X2 = build_agg_df(temp, total_reps)
            
            X = pd.merge(X1, X2, on='LipidMolec')
            
            X = X[['LipidMolec', 'ClassKey', 'FAKey', 'CalcMass', 'BaseRt'] + ['MeanArea[s' + str(i+1) + ']' for i in range(total_reps)]]
            
            # extracting internal standards 
            
            X, intsta_df = extract_internal_standards(X)
            
            #third filter: minimum abundance grade 
        
            X = build_abundance_grade_filter(X, rep_lst, cond_lst, filtered_conds, passing_abundance_grade)
            
            # organizing and printing the dataset
            
            X.reset_index(inplace=True)
            
            X.drop(['index'], axis=1, inplace=True) # drops an irrelevant column
            
            cleaned_data = st.checkbox("View Cleaned Data")
            
            if cleaned_data:
            
                st.write('View the cleaned data in the conventional format:')
                
                st.write(X)
            
                csv_download = convert_df(X)
                            
                st.download_button(
                                label="Download Data",
                                
                                data=csv_download,
                                
                                file_name='cleaned_data.csv',
                                
                                mime='text/csv')
            
                st.write('------------------------------------------------------------------------------------------------')
            
                st.write('View the cleaned data in the log-transformed format:')
                
                log_X = log_transform_df(X, rep_lst)
            
                st.write(log_X)
            
                csv_download = convert_df(log_X)
                            
                st.download_button(
                    
                                label="Download Data",
                                
                                data=csv_download,
                                
                                file_name='log_transformed_cleaned_data.csv',
                                
                                mime='text/csv')
                
                st.write('------------------------------------------------------------------------------------------------')
                
                st.write('View the internal standards dataset in the conventional format:')
                
                st.write(intsta_df)
            
                csv_download = convert_df(intsta_df)
                            
                st.download_button(
                                label="Download Data",
                                
                                data=csv_download,
                                
                                file_name='internal_standards.csv',
                                
                                mime='text/csv')
            
                st.write('------------------------------------------------------------------------------------------------')
            
                st.write('View the internal standards dataset in the log-transformed format:')
                
                log_intsta_df = log_transform_df(intsta_df, rep_lst)
            
                st.write(log_intsta_df)
            
                csv_download = convert_df(log_intsta_df)
                            
                st.download_button(
                    
                                label="Download Data",
                                
                                data=csv_download,
                                
                                file_name='log_transformed_intsta_df.csv',
                                
                                mime='text/csv')
            
            return X 
        
        
        
        
        
    @st.cache
    def build_agg_df(dataframe, total_reps):
        
        X1 = dataframe.groupby(["LipidMolec", 'ClassKey', 'FAKey'])[["MeanArea[s" + str(i+1) + ']' for i in range(total_reps)]].sum()
        
        X1.reset_index(inplace = True)
        
        X2 = dataframe.groupby(["LipidMolec"])[['CalcMass', 'BaseRt']].mean()
        
        X2.reset_index(inplace = True)
        
        return X1, X2 
        
        
        
        
        
# function to log transform the abundance columns of the df
    @st.cache
    def log_transform_df(X, rep_lst):  
        
            temp = X.copy()
            
            auc = ['MeanArea[s'+str(i+1)+']' for i in range(sum(rep_lst))]
            
            # filling zero values with 1's to avoid infinity
            
            temp[auc]=temp[auc].mask(temp[auc]<=1).fillna(1)
            
            temp[auc] = temp[auc].apply(lambda x: np.log10(x), axis=0)
            
            return temp





 # function that removes the lipid species that have an abundance grade lower than the minimum abundance grade
    def build_abundance_grade_filter(X, rep_lst, cond_lst, filtered_conds, passing_abundance_grade):
        
        rep_lst_agg = build_rep_lst_agg(rep_lst)  
        
        for cond in filtered_conds:
            
            sample_lst = build_sample_lst(cond, cond_lst, rep_lst_agg)
            
            # adding the 'abundance grade' column for each filtered condition 
            
            X['abundance_grade_'+cond] = X[['MeanArea['+sample+']' for sample in sample_lst]]\
                .apply(lambda x: abundance_level_calculator(x, passing_abundance_grade),axis=1)
            
        X['abundance_grade'] = X[['abundance_grade_'+cond for cond in filtered_conds]].apply(lambda x: sum(x), axis=1) # adds the 'total abundance grade' column 
        
        X = X.loc[X['abundance_grade'] > 0 ] # keeps only lipid species with non-zero 'total abundance grade'  
        
        X.drop(['abundance_grade_'+cond for cond in filtered_conds]+['abundance_grade'], axis=1, inplace=True) # drops the abundance grade columns
        
        return X
    
    
    
    
    
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
            
            temp = temp[['MeanArea[s'+str(i+1)+']' for i in range(sum(rep_lst))]] # picks the 'MainArea[s1]', ..., 'MainArea[sN]' columns only
            
            temp = temp.mask(temp<=1).fillna(1)
            
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
                    
                        arrange_hist(diff, temp[['MeanArea[s'+str(i+1)+']' for i in range(srange[0], srange[1])]], full_sample_lst[srange[0]: srange[1]])
                    
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
        
        lst = np.log10(temp['MeanArea[' + full_sample_lst[index] + ']'].values.tolist())
            
        ax.hist(lst, bins = 75, range=(0, 12))
                    
        ax.set_title('Histogram of AUC - '+ full_sample_lst[index], fontsize=50)
                            
        ax.set_xlabel('log10(AUC)', fontsize=50)
                            
        ax.set_ylabel('Count', fontsize=50)
                        
        ax.set_xlim([-0.5, 12])
                        
        return fig
    
    
    
    
    
# all retention time plots
    def plot_retention(X, cond_lst, rep_lst):  
    
            temp = X.copy()
            
            show_retention = st.checkbox("View the retention time plots")
            
            if show_retention:
                
                mode = st.radio('', ['Comparison Mode', 'Individual Mode'])
                
                if mode == 'Individual Mode':
                
                    lipid_class_lst = temp['ClassKey'].value_counts().index.tolist()
                
                    for lipid_class in lipid_class_lst:
                    
                        class_name = lipid_class
                    
                        plot, retention_df = retention_hover(temp, rep_lst, cond_lst, class_name)
                    
                        st.bokeh_chart(plot)
                    
                        file_name = 'retention_plot_' + class_name
                        
                        csv_download = convert_df(retention_df)
                            
                        st.download_button(
                            
                                label="Download Data",
                                
                                data=csv_download,
                                
                                file_name=file_name+'.csv',
                                
                                mime='text/csv')
                    
                        st.write("----------------------------------------------------------------------------------------------------")
                        
                else:
                    
                    lipid_class_lst = temp['ClassKey'].value_counts().index.tolist()
                    
                    selected_class = st.multiselect('Add or remove classes (up to 20 classes):', lipid_class_lst, lipid_class_lst[:2])
                    
                    if len(selected_class) > 20:
                        
                        st.error('You can only compare up to 20 lipid classes at a time!')
                        
                        return None
                        
                    else:
                    
                        plot, retention_df = retention_multi(temp, rep_lst, selected_class)
                    
                        st.bokeh_chart(plot)
                        
                        csv_download = convert_df(retention_df)
                            
                        st.download_button(
                            
                                label="Download Data",
                                
                                data=csv_download,
                                
                                file_name='Retention_Time_Comparison_Mode.csv',
                                
                                mime='text/csv')
                    
            return
        
        
        
        
        
    # plots a single retention time plot with hover tool
    def retention_hover(temp, rep_lst, cond_lst, class_name):  
            
            retention = temp[temp['ClassKey'] == class_name]['BaseRt'].values.tolist()
            
            mass = temp[temp['ClassKey'] == class_name]['CalcMass'].values.tolist()
            
            species = temp[temp['ClassKey'] == class_name]['LipidMolec'].values.tolist()
            
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
    def retention_multi(temp, rep_lst, selected_class): 
        
        unique_color_lst =[ 'red', 'blue', 'green', 'magenta', 'cyan', 'orange', 'black', 'pink', 'brown', 'yellow', 'purple', \
                            'gray', 'olive', 'chocolate', 'silver', 'darkred', 'khaki', 'skyblue', 'navy', 'orchid']
        
        retention = []
        
        mass = []
        
        class_lst = []
        
        color_lst = []
        
        for lipid_class in selected_class:
            
            retention = retention + temp[temp['ClassKey'] == lipid_class]['BaseRt'].values.tolist()
            
            mass = mass + temp[temp['ClassKey'] == lipid_class]['CalcMass'].values.tolist()
            
            class_lst = class_lst + temp[temp['ClassKey'] == lipid_class]['ClassKey'].values.tolist()
            
            color_lst = color_lst + [unique_color_lst[selected_class.index(lipid_class)] for i in range(len(temp[temp['ClassKey'] == lipid_class]))]
        
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
        
                temp = X[['MeanArea[s' + str(i+1) + ']' for i in range(sum(rep_lst))]].copy()
        
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
    
            temp = temp[['MeanArea[s'+str(i+1)+']' for i in range(sum(replicate_lst))] + ['ClassKey']]
            
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
                    
                    lipid_class = st.selectbox('Pick a lipid class', temp['ClassKey'].value_counts().index.tolist())
                    
                    temp = temp[temp['ClassKey'] == lipid_class]
            
                temp.drop(['ClassKey'], axis=1, inplace=True)
               
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
            
            dataframe = dataframe.mask(dataframe<=1).fillna(1)
            
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
        
        dataframe.drop(['MeanArea[' + sample + ']' for sample in r_sample], axis=1, inplace=True) # updating dataframe
            
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
            
            dataframe.rename(columns={'MeanArea[' + sample_1 + ']' : 'Area[' + sample_2 + ']'}, inplace = True)
            
        for i in range(total_reps):
            
            dataframe.rename(columns={'Area[s' + str(i+1) + ']' : 'MeanArea[s' + str(i+1) + ']'}, inplace = True)
        
        return replicate_lst, condition_lst, name_df, dataframe
    
    
    
    
    
# function to plot the CoV of lipid species
    def plot_cov(X, cond_lst, rep_lst):  
        
        auc = ['MeanArea[s'+str(i+1)+']' for i in range(sum(rep_lst))]
        
        cond = st.radio('Which of the following conditions corresponds to BQC samples? ', cond_lst+['None of the above'], len(cond_lst))
            
        if cond != 'None of the above':
                
                index = cond_lst.index(cond)
                
                if rep_lst[index] == 1:
                    
                    st.error('The selected condition must have at least two corresponding replicates!')
                    
                else:
                    
                    X[auc]=X[auc].mask(X[auc]<=1).fillna(1) # turning 0's and negative numbers  to 1's so it is possible to log transform
            
                    rep_lst_agg = build_rep_lst_agg(rep_lst)
                
                    sample_lst = build_sample_lst(cond, cond_lst, rep_lst_agg)
                    
                    X, X_plot, X_cov_df = cov_hover(X, sample_lst)
                    
                    st.write('--------------------------------------------------------------------------------------------------------------')
                    
                    show_cov = st.checkbox("View CoV plot")
                    
                    if show_cov:
                
                        st.bokeh_chart(X_plot)
                    
                        csv_download = convert_df(X_cov_df)
                            
                        st.download_button(
                        
                                label="Download Data",
                                
                                data=csv_download,
                                
                                file_name='cov.csv',
                                
                                mime='text/csv')
                        
                    st.write('--------------------------------------------------------------------------------------------------------------')
                    
                    filter_ans = st.radio('Would you like to filter the data by removing datapoints with high CoV?', ['Yes', 'No'], 1)
                    
                    st.warning("Your choice here would affect the rest of the module.")
                    
                    if filter_ans == 'Yes':
                        
                        thresh = st.number_input('Enter the maximum acceptable CoV in %', min_value = 10, max_value = 100, value = 30, step = 1)
                        
                        X = X.loc[X['cov'] <= thresh] # removes datapoints with CoV > thresh
                    
                        X[auc]=X[auc].mask(X[auc]==1).fillna(0) # turning 1's back to 0's
                    
                        X.drop(['mean', 'cov'], axis=1, inplace=True)
                        
                        filtered_data = st.checkbox('View Filtered Data')
                        
                        if filtered_data:
                    
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
            
        X['cov'] = X[['MeanArea['+sample+']' for sample in sample_lst]].apply(lambda x: cov_calculator(x), axis=1)
    
        X['mean'] = X[['MeanArea['+sample+']' for sample in sample_lst]].apply(lambda x: mean_calculator(x), axis=1)
            
        plot = figure(title='CoV - All lipid Species', x_axis_label='Mean of log10(AUC) Across All BQC Samples', y_axis_label='CoV(%)')
            
        x = X['mean'].values.tolist()
            
        y = X['cov'].values.tolist()
        
        species = X['LipidMolec'].values.tolist()
            
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
    
        if len(non_zero_lst) > 1:

            cov = np.std(non_zero_lst)/np.mean(non_zero_lst)*100
        
        else:
        
            cov = None 
        
        return cov
    
    
    
    
    
 # calculate mean
    @st.cache
    def mean_calculator(numbers): 
    
        non_zero_lst = [number for number in numbers if (number>1)]
    
        if len(non_zero_lst) > 1:

            mean = np.mean(non_zero_lst)
    
            mean = np.log10(mean)
        
        else:
        
            mean = None
        
        return mean
    
    
    
    
    
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
        
        
        
        
        
# function for imputing missing values 
    def impute_missing_value(X):
    
        full_sample_lst = ['s'+str(i+1) for i in range(sum(rep_lst))]
                
        for sample in full_sample_lst:
                
                lst = [ele for ele in X['MeanArea[' + sample + ']'].values if ele > 0]
                
                impute_value = min(lst)
                
                X['MeanArea[' + sample + ']'] = X['MeanArea[' + sample + ']'].apply(lambda x: impute_value if x<=0 else x)
        
        return X





    # function that creates all volcano plots
    def volcano_plot(X, cond_lst, rep_lst):
    
        temp = X.copy()
    
        # different options for imputing 
        
        impute_ans = st.radio('Pick one option',\
                              ['No imputing', 'Replace missing values in each sample by the minimum detected AUC in that sample'], 1)
            
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
                    temp[['MeanArea[' + sample + ']'  for sample in sample_lst_1 + sample_lst_2]]\
                    .apply(lambda x: p_val_calculator(x[['MeanArea[' + sample + ']'  for sample in sample_lst_1]], \
                                                      x[['MeanArea[' + sample + ']'  for sample in sample_lst_2]]), axis=1)
                        
            temp['fc_' + cond_1 + '_' + cond_2] = \
                    temp[['MeanArea[' + sample + ']'  for sample in sample_lst_1 + sample_lst_2]]\
                    .apply(lambda x: fc_calculator(x[['MeanArea[' + sample + ']'  for sample in sample_lst_1]], \
                                                   x[['MeanArea[' + sample + ']'  for sample in sample_lst_2]]), axis=1)
                
            lipid_class_lst = temp['ClassKey'].value_counts().index.tolist()
                    
            selected_class = st.multiselect('Add or remove classes:', lipid_class_lst, lipid_class_lst)
            
            if len(selected_class) > 30:
                        
                st.error('You can only compare up to 30 lipid classes at a time!')
                        
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
                            'gray', 'olive', 'chocolate', 'silver', 'darkred', 'khaki', 'skyblue', 'navy', 'orchid', 'deepskyblue', \
                            'tan', 'lime', 'tomato', 'deeppink', 'lavender', 'salmon', 'darkkhaki', 'lavenderblush', 'indigo']
        
        fc = []
        
        pval = []
        
        class_lst = []
        
        color_lst = []
        
        species = []
        
        for lipid_class in selected_class:
            
            fc = fc + dataframe[dataframe['ClassKey'] == lipid_class]['fc_' + cond_1 + '_' + cond_2].values.tolist()
            
            pval = pval + dataframe[dataframe['ClassKey'] == lipid_class]['p_val_' + cond_1 + '_' + cond_2].values.tolist()
            
            color_lst = color_lst + [unique_color_lst[selected_class.index(lipid_class)] for i in range(len(dataframe[dataframe['ClassKey'] == lipid_class]))]
                
            class_lst = class_lst + dataframe[dataframe['ClassKey'] == lipid_class]['ClassKey'].values.tolist()
        
            species = species + dataframe[dataframe['ClassKey'] == lipid_class]['LipidMolec'].values.tolist()
            
        plot = figure(title='Volcano Plot', x_axis_label='Fold Change (' + cond_1 + '/' + cond_2 + ')', y_axis_label='q-value')
            
        vol_df = pd.DataFrame({"FC": fc, "qvalue": -np.log10(pval), "Species": species, "Class": class_lst, \
                               "Color": color_lst})
        
        src = ColumnDataSource(vol_df)
        
        plot.scatter(x="FC", y="qvalue", legend_group='Class', color='Color', name='volcano', size = 4, source=src)
        
        plot.line(x=[i for i in range(-10, 11)], y=-np.log10(0.05), line_dash = 'dashed', color='black')
        
        plot.line(x=-1, y=[i for i in range(0, 9)], line_dash = 'dashed', color='black')
        
        plot.line(x=1, y=[i for i in range(0, 9)], line_dash = 'dashed', color='black')
            
        hover = HoverTool(tooltips = [('FC', '@FC'), ('p-value', '@qvalue'), ('Species', '@Species')], names=['volcano'])
            
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
    
    
    
    
    
# function to perform saturation level analysis 
    def saturation_level_plot(X, cond_lst, rep_lst):
    
        temp = X.copy()
        
        auc = ['MeanArea[s'+str(i+1)+']' for i in range(sum(rep_lst))]
        
        temp[auc]=temp[auc].mask(temp[auc]<=1).fillna(1)
                    
        rep_lst_agg = build_rep_lst_agg(rep_lst)
        
        show_typical_bar = st.checkbox("View Typical Bar Plots")
        
        if show_typical_bar:
            
            saturation_level_plot_bar(temp, rep_lst_agg, cond_lst)
                
        show_stacked_bar_plots = st.checkbox("View Stacked Bar Plots")
        
        if show_stacked_bar_plots:
            
            saturation_level_plot_stack(temp, rep_lst_agg, cond_lst)
                
        show_stacked_percentage_bar_plots = st.checkbox("View Percentage-Based Stacked Bar Plots")
        
        if show_stacked_percentage_bar_plots:
            
            saturation_level_plot_percentage(temp, rep_lst_agg, cond_lst)
            
        return 
    
    
    
    
        
    def saturation_level_plot_bar(dataframe, rep_lst_agg, cond_lst):
    
        for lipid_class in dataframe['ClassKey'].unique():
        
            sfa_lst = []
        
            mufa_lst = []
        
            pufa_lst = []
        
            sfa_var_lst = []
        
            mufa_var_lst = []
        
            pufa_var_lst = []
        
            for cond in cond_lst:
            
                sample_lst = build_sample_lst(cond, cond_lst, rep_lst_agg)
                
                sfa, mufa, pufa, sfa_var, mufa_var, pufa_var = calculate_SFA_MUFA_PUFA(dataframe, lipid_class, cond, sample_lst)
            
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
            
            sat_df = pd.DataFrame({"Conditions": cond_lst, "SFA": sfa_lst, "SFA_STDV": sfa_var_lst, "MUFA": mufa_lst, \
                           "MUFA_STDV": mufa_var_lst, "PUFA": pufa_lst, "PUFA_STDV": pufa_lst})
            
            csv_download = convert_df(sat_df)
                        
            st.download_button(
                
                            label="Download Data",
                            
                            data=csv_download,
                            
                            file_name='saturation_level_plot_bar.csv',
                            
                            mime='text/csv',
                            
                            key=lipid_class)
            
            st.write('------------------------------------------------------------------------------------------------')
        
        return 





    def saturation_level_plot_stack(dataframe, rep_lst_agg, cond_lst):
    
        for lipid_class in dataframe['ClassKey'].unique():
            
            sfa_lst = []
            
            mufa_lst = []
            
            pufa_lst = []
            
            for cond in cond_lst:
                
                sample_lst = build_sample_lst(cond, cond_lst, rep_lst_agg)
                
                sfa, mufa, pufa, sfa_var, mufa_var, pufa_var = calculate_SFA_MUFA_PUFA(dataframe, lipid_class, cond, sample_lst)
                
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
            
            csv_download = convert_df(sat_df)
                        
            st.download_button(
                
                            label="Download Data",
                            
                            data=csv_download,
                            
                            file_name='saturation_level_plot_stacked.csv',
                            
                            mime='text/csv',
                            
                            key=lipid_class)
            
            st.write('------------------------------------------------------------------------------------------------')
        
        return 





    def saturation_level_plot_percentage(dataframe, rep_lst_agg, cond_lst):
    
        for lipid_class in dataframe['ClassKey'].unique():
            
            sfa_lst = []
            
            mufa_lst = []
            
            pufa_lst = []
            
            for cond in cond_lst:
                
                sample_lst = build_sample_lst(cond, cond_lst, rep_lst_agg)
                
                sfa, mufa, pufa, sfa_var, mufa_var, pufa_var = calculate_SFA_MUFA_PUFA(dataframe, lipid_class, cond, sample_lst)
                
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
            
            csv_download = convert_df(sat_df)
                        
            st.download_button(
                
                            label="Download Data",
                            
                            data=csv_download,
                            
                            file_name='saturation_level_plot_percentage.csv',
                            
                            mime='text/csv',
                            
                            key = lipid_class)
            
            st.write('------------------------------------------------------------------------------------------------')
        
        return 
    
    
    
    
    
    # function to calculate the ratio of sfa, mufa and pufa for each lipid species 
    @st.cache
    def calculate_FA_ratio(mol_structure):
            
        a = mol_structure.split('(')
        
        b = a[1][:-1]
        
        c = b.split('_')
        
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
    
    
    
    
    
 # function to calculate the abundance of sfa, mufa and pufa for each lipid species 
    @st.cache
    def calculate_SFA_MUFA_PUFA(dataframe, lipid_class, cond, sample_lst):
        
        dataframe = dataframe[dataframe['ClassKey'] == lipid_class]
        
        dataframe[cond + '_mean_AUC'] = dataframe[['MeanArea[' + sample + ']' for sample in sample_lst]]\
            .apply(lambda x: np.mean(x), axis = 1)
            
        dataframe[cond + '_var_AUC'] = dataframe[['MeanArea[' + sample + ']' for sample in sample_lst]]\
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
    
    
    
    
    
    def plot_pathway(X, cond_lst, rep_lst):
    
        pathway_viz = st.checkbox("Visualize Lipidomic Pathways")
        
        if pathway_viz:
        
            temp = X.copy()
            
            control = st.selectbox('Select the control group', cond_lst)
            
            experimental = st.selectbox('Select the group you would like to compare with the control group (experimental group)', cond_lst)
            
            # ratio_lst: list of abundance ratio between any pair of classes 
            
            class_lst, ratio_lst = calculate_class_fold_change(X, cond_lst, rep_lst, control, experimental)
            
            temp['saturation_comb'] = temp['LipidMolec'].apply(lambda x: calculate_saturated_unsaturated_chain_number(x))
            
            # what percentage of chains are saturated 
            
            sat_ratio_lst = calculate_saturation_ratio(temp, class_lst)
            
            plt.rcParams["figure.figsize"] = [10, 10]
        
            fig, ax = plt.subplots()
        
            ax.set_title('Lipid Pathway Visualization', fontsize=20)
        
            ax.set_xlim([-25, 25])
            ax.set_ylim([-20, 30])
        
            ax.axes.xaxis.set_visible(False)
            ax.axes.yaxis.set_visible(False)
        
            plt.gca().set_aspect('equal', adjustable='box')
            
            ax.add_patch(draw_circle(5, 0, 0, 'b'))
            ax.add_patch(draw_circle(2.5, -7.5 * math.cos(math.pi/6), -7.5 * math.cos(math.pi/3), 'b'))
            ax.add_patch(draw_circle(2.5, 7.5 * math.cos(math.pi/6), -7.5 * math.cos(math.pi/3), 'b'))
            ax.add_patch(draw_circle(2.5, 10 + 2.5 * math.cos(math.pi/4), 15 + 2.5 * math.sin(math.pi/4), 'b'))
            ax.add_patch(draw_circle(2.5, 12.5, 10, 'b'))
            ax.add_patch(draw_circle(2.5, 10 + 2.5 * math.cos(math.pi/4), 5 - 2.5 * math.sin(math.pi/4), 'b'))
        
            ax.plot([0, 0], [0, 20], c='b')
            ax.plot([0, -5], [15, 20], c='b')
            ax.plot([-5, -10], [20, 15], c='b')
            ax.plot([-10, -10], [15, 5], c='b')
            ax.plot([0, 5], [10, 10], c='b')
            ax.plot([5, 10], [10, 15], c='b')
            ax.plot([5, 10], [10, 10], c='b')
            ax.plot([5, 10], [10, 5], c='b')
            ax.plot([5*math.cos(math.pi/6), 10], [-5*math.sin(math.pi/6), 5], c='b')
            
            ax.annotate('G3P', xy=(0, 20), xytext=(0, 23),arrowprops=dict(facecolor='black'), fontsize=15)
            ax.annotate('Fatty Acids', xy=(-5, 20), xytext=(-5, 25),arrowprops=dict(facecolor='black'), fontsize=15)
            
            ax.text(-3.5, 14, 'LPA', fontsize=15)
            ax.text(-2.5, 9.5, 'PA', fontsize=15)
            ax.text(-4, 5.5, 'DAG', fontsize=15)
            ax.text(-4, 0.5, 'TAG', fontsize=15)
            ax.text(-4, -2, 'PC', fontsize=15)
            ax.text(-12, -6.5, 'LPC', fontsize=15)
            ax.text(2.5, -2, 'PE', fontsize=15)
            ax.text(9, -6, 'LPE', fontsize=15)
            ax.text(-14, 15, 'LCBs', fontsize=15)
            ax.text(-13.5, 10, 'Cer', fontsize=15)
            ax.text(-13, 5, 'SM', fontsize=15)
            ax.text(2, 11, 'CDP-DAG', fontsize=15)
            ax.text(10.5, 15.5, 'PI', fontsize=15)
            ax.text(14, 19, 'LPI', fontsize=15)
            ax.text(10.5, 9.5, 'PG', fontsize=15)
            ax.text(15.5, 9.5, 'LPG', fontsize=15)
            ax.text(10.5, 4, 'PS', fontsize=15)
            ax.text(14.5, 1, 'LPS', fontsize=15)
            
            pathway_class_lst = ['TG', 'DG', 'PA', 'LPA', 'LCB', 'Cer', 'SM', 'PE', 'LPE', 'PC', 'LPC', 'PI', 'LPI', 'CDP-DAG', \
                                     'PG', 'LPG', 'PS', 'LPS']
                
            pathway_ratio_lst = [0] * len(pathway_class_lst)
            
            pathway_sat_ratio_lst = [0] * len(pathway_class_lst)
            
            i = 0
            
            for lipid_class in pathway_class_lst:
                
                if lipid_class in class_lst:
                    
                    class_index = class_lst.index(lipid_class)
                    
                    pathway_ratio_lst[i] = ratio_lst[class_index]
                    
                    pathway_sat_ratio_lst[i] = sat_ratio_lst[class_index]
                
                i += 1
                
            pathway_dict = {'class': pathway_class_lst, 'abundance ratio': pathway_ratio_lst, \
                            'saturated fatty acids ratio': pathway_sat_ratio_lst}
                    
            pathway_ratio_lst = [50*ele**2 for ele in pathway_ratio_lst]
            
            color_contour = pathway_sat_ratio_lst
            size = pathway_ratio_lst
            
            points_x = [0, 0, 0, 0, -10, -10, -10, 5*math.cos(math.pi/6), 10*math.cos(math.pi/6), \
                        -5*math.cos(math.pi/6), -10*math.cos(math.pi/6), 10, 10+5*math.cos(math.pi/4), \
                        5, 10, 15, 10, 10+5*math.cos(math.pi/4)]
            points_y = [0, 5, 10, 15, 15, 10, 5, -5*math.sin(math.pi/6), -10*math.sin(math.pi/6), \
                        -5*math.sin(math.pi/6), -10*math.sin(math.pi/6), 15, 15+5*math.sin(math.pi/4), \
                        10, 10, 10, 5, 5-5*math.sin(math.pi/4)]
            
            ax.add_patch(draw_circle(0.5, 0, 0, 'black'))
            ax.add_patch(draw_circle(0.5, 0, 5, 'black'))
            ax.add_patch(draw_circle(0.5, 0, 10, 'black'))
            ax.add_patch(draw_circle(0.5, -10, 10, 'black'))
            ax.add_patch(draw_circle(0.5, -10, 5, 'black'))
            ax.add_patch(draw_circle(0.5, 5*math.cos(math.pi/6), -5*math.sin(math.pi/6), 'black'))
            ax.add_patch(draw_circle(0.5, 10*math.cos(math.pi/6), -10*math.sin(math.pi/6), 'black'))
            ax.add_patch(draw_circle(0.5, -5*math.cos(math.pi/6), -5*math.sin(math.pi/6), 'black'))
            ax.add_patch(draw_circle(0.5, -10*math.cos(math.pi/6), -10*math.sin(math.pi/6), 'black'))
            ax.add_patch(draw_circle(0.5, 10, 15, 'black'))
            ax.add_patch(draw_circle(0.5, 10+5*math.cos(math.pi/4), 15+5*math.sin(math.pi/4), 'black'))
            ax.add_patch(draw_circle(0.5, 10, 10, 'black'))
            ax.add_patch(draw_circle(0.5, 15, 10, 'black'))
            ax.add_patch(draw_circle(0.5, 10, 5, 'black'))
            ax.add_patch(draw_circle(0.5, 10+5*math.cos(math.pi/4), 5-5*math.sin(math.pi/4), 'black'))
        
            points = ax.scatter(points_x, points_y, c=color_contour, s=size, cmap="plasma")
            
            cbar = fig.colorbar(points)
            
            cbar.set_label(label='Unsaturated <--- Saturatted Fatty Acids Ratio ---> Saturated', size=15)
            
            cbar.ax.tick_params(labelsize=15)
            
            st.pyplot(fig)
                
            pathway_df = pd.DataFrame.from_dict(pathway_dict)
            
            pathway_df.set_index('class', inplace=True)
            
            st.write(pathway_df)
            
            csv_download = convert_df(pathway_df)
                        
            st.download_button(
                            label="Download Data",
                            
                            data=csv_download,
                            
                            file_name='pathway_viz.csv',
                            
                            mime='text/csv')
        
        return 





    def draw_circle(radius, x0, y0, color):
    
        return plt.Circle((x0, y0), radius, color = color, fill = False)





    def calculate_class_fold_change(X, cond_lst, rep_lst, control, experimental):
    
        temp = X.copy()
        
        total_reps = sum(rep_lst)
        
        rep_lst_agg = build_rep_lst_agg(rep_lst)
        
        group_temp = temp.groupby('ClassKey')[["MeanArea[s" + str(i+1) + ']' for i in range(total_reps)]].sum()
        
        control_index = cond_lst.index(control)
        
        experimental_index = cond_lst.index(experimental)
        
        class_lst = group_temp.index.tolist()
        
        ratio_lst = [0] * len(class_lst)
        
        i = 0
        
        for lipid_class in group_temp.index:
            
            if control_index == 0:
        
                denom = np.mean(group_temp.iloc[i][: rep_lst_agg[0]])
                
            else:
                
                denom = np.mean(group_temp.iloc[i][rep_lst_agg[control_index - 1] : rep_lst_agg[control_index]])
                
            if experimental_index == 0:
                
                nume = np.mean(group_temp.iloc[i][: rep_lst_agg[0]])
                
            else:
                
                nume = np.mean(group_temp.iloc[i][rep_lst_agg[experimental_index - 1] : rep_lst_agg[experimental_index]])
                
            ratio = nume / denom 
            
            ratio_lst[i] = ratio
            
            i += 1
        
        return class_lst, ratio_lst





    def calculate_saturated_unsaturated_chain_number(mol_structure):
            
        a = mol_structure.split('(')
        
        b = a[1][:-1]
        
        c = b.split('/')
        
        sat_chain = 0
        
        unsat_chain = 0
        
        for item in c:
            
            d = item.split(':')[-1]
            
            if d.isnumeric():
                
                d = int(d)
            
                if d == 0:
                
                    sat_chain += 1 
                
                else:
                
                    unsat_chain += 1
                    
            else:
                
                if '0' in d:
                    
                    sat_chain += 1
                    
                else:
                    
                    unsat_chain += 1
        
        return sat_chain, unsat_chain
    
    
    
    
    
    def calculate_saturation_ratio(dataframe, class_lst):
    
        dataframe['total_sat'] = dataframe['saturation_comb'].apply(lambda x: x[0])
        
        dataframe['total_unsat'] = dataframe['saturation_comb'].apply(lambda x: x[1])
        
        group_dataframe = dataframe.groupby('ClassKey')[['total_sat', 'total_unsat']].sum()
        
        group_dataframe['saturation_ratio'] = group_dataframe[['total_sat', 'total_unsat']].apply(lambda x: x[0]/(x[0]+x[1]), axis = 1)
        
        sat_ratio_lst = group_dataframe['saturation_ratio'].tolist()
        
        return sat_ratio_lst





    def plot_total_abundance(X, cond_lst, rep_lst):
    
        view_abundance_chart = st.checkbox("View Total Abundance Bar Charts")
        
        if view_abundance_chart: 
        
            temp = X.copy()
            
            total_reps = sum(rep_lst)
            
            rep_lst_agg = build_rep_lst_agg(rep_lst)
            
            group_temp = temp.groupby('ClassKey')[["MeanArea[s" + str(i+1) + ']' for i in range(total_reps)]].sum()
            
            class_lst = group_temp.index.tolist()
            
            for cond in cond_lst:
                
                index = cond_lst.index(cond)
                
                if index == 0:
                
                    group_temp["mean_AUC_" + cond] = group_temp[["MeanArea[s" + str(i+1) + ']' for i in range(rep_lst_agg[0])]].apply(lambda x: np.mean(x), axis = 1)
                    
                    group_temp["std_AUC_" + cond] = group_temp[["MeanArea[s" + str(i+1) + ']' for i in range(rep_lst_agg[0])]].apply(lambda x: np.std(x), axis = 1)
                    
                else:
                    
                    group_temp["mean_AUC_" + cond] = group_temp[["MeanArea[s" + str(i+1) + "]" for i in range(rep_lst_agg[index - 1], rep_lst_agg[index])]]\
                                 .apply(lambda x: np.mean(x), axis = 1)
                                 
                    group_temp["std_AUC_" + cond] = group_temp[["MeanArea[s" + str(i+1) + "]" for i in range(rep_lst_agg[index - 1], rep_lst_agg[index])]]\
                                 .apply(lambda x: np.std(x), axis = 1)
            
            for cond in cond_lst:
                
                auc = group_temp['mean_AUC_' + cond].tolist()
                
                std = group_temp['std_AUC_' + cond].tolist()
                
                log_auc = np.log10(auc)
                
                log_std = [(np.log10(auc[i]+std[i]) - np.log10(auc[i]-std[i])) for i in range(len(auc))] 
                
                y_pos = np.arange(len(class_lst))
                
                plt.rcdefaults()
                fig, ax = plt.subplots()
                
                ax.barh(y_pos, log_auc, xerr=log_std, align='center')
                ax.set_yticks(y_pos)
                ax.set_yticklabels(class_lst)
                ax.invert_yaxis()  # labels read top-to-bottom
                ax.set_xlabel('log10(Total AUC) averaged over all replicates')
                ax.set_ylabel('Lipid Class')
                ax.set_title('Total Abundance Bar Chart - ' + cond)
                
                st.pyplot(fig)
                
                abundance_df = pd.DataFrame({"class": class_lst, "total_AUC": log_auc, "total_AUC_STDV": log_std})
                
                csv_download = convert_df(abundance_df)
                            
                st.download_button(
                    
                                label="Download Data",
                                
                                data=csv_download,
                                
                                file_name='abundance_bar_chart.csv',
                                
                                mime='text/csv')
                
                st.write('------------------------------------------------------------------------------------------------')
    
        return

##########################################################################################################################################
# the main code of the app 

    st.header("LipidSearch 5.0 Module")

    st.markdown("""
            
            The following module allows the user to run lipidomics analysis on LipidSearch 5.0 datasets. 
            Start by uploading your dataset and completing the next steps on the side bar.
            
            """)
            
    st.info("""
        
        A standard LipidSearch dataset must have the following columns:

        LipidMolec: the class that the lipid species belong to and its number of carbon atoms and double bonds

        ClassKey: the class that the lipid species belong to
        
        FAKey: the number of carbon atoms and double bonds corresponding to the lipid species 

        CalcMass: the calculated mass of the lipid species

        BaseRt: the retention time of the lipid species in the chromatography column

        MeanArea[s1], ..., MeanArea[sN]: Area Under the Curve (AUC) representing 
        the relative abundance of the lipid species in samples s1 to sN where N stands for the total number of the sampels
        
        """)

    st.sidebar.subheader("Upload Data")
            
    lipid_search = st.sidebar.file_uploader(label='Upload your LipidSearch 5.0 dataset', type=['csv', 'txt'])
            
    if lipid_search is not None:
        
        df = build_lipidsearch_df(lipid_search)
        
        if df is not None:
            
            # building the side bar 
        
            confirm_data, name_df, filtered_conds, passing_abundance_grade, cond_lst, rep_lst = build_sidebar(df)
                
            # if user confirms the inputs:
    
            if confirm_data:
                
                st.subheader("1) Data Cleaning, Exploration, Quality Check & Anomaly Detection")
                        
                st.markdown("""
                                    The Data Cleaning, Exploration, Quality Check & Anomaly Detection submodule allows the user to clean, filter
                                    and visualize the data, and runs anomaly detection tests on it.
                    
                                    """)
            
                st.subheader("1.1) Clean, Filter & Explore Data")
            
                expand_raw_data = st.expander("Raw Data")
    
                with expand_raw_data:
                    
                    raw_data = st.checkbox('View Raw Data')
                    
                    if raw_data:
        
                        st.write(df)
            
                        csv_download = convert_df(df)
                    
                        st.download_button(
                        
                            label="Download Data",
                        
                            data=csv_download,
                        
                            file_name='data.csv',
                        
                            mime='text/csv') 
                        
                expand_cleaned_data = st.expander("Cleaned Data")
                
                with expand_cleaned_data:
                    
                    st.markdown("""
                    
                    The data cleaning process is a four steps process: 
                        
                    1) LipidCruncher deletes the datapoints that do not pass through the filters: 
                        either their associated "Total Grade" value is "C" or "D", or they have too many missing values.
                    
                    2) LipidCruncher only keeps the relevant columns.
                    
                    3) 'MeanArea' columns in the cleaned dataset are the aggregated version of those in the raw dataset.
                        In the raw dataset, each lipid species appears as multiple datapoints. For each lipid species, LipidCruncher takes the "sum" 
                        of AUC over all the corresponding datapoints. Therefore, in the cleaned dataset, each lipid species is represented by a single 
                        datapoint. 
                        
                    4) In the raw dataset, each lipid species appears as multiple datapoints. For each lipid species, LipidCruncher takes the "mean" 
                        of the retention times over all the corresponding datapoints. Therefore, in the cleaned dataset, there is only one retention time 
                        corresponding to each lipid species and that is the mean retention time. 
                    
                    """)
            
                    X = clean_data(df, rep_lst, cond_lst, name_df, filtered_conds, passing_abundance_grade) # cleaned data 
                    
                expand_filtered_data = st.expander("Filter Data Using BQC Samples")
                    
                with expand_filtered_data:
                    
                    st.markdown(""" 
                                
                        Filtering mass spectrometry data is challenging. However, one reliable method to filter the data 
                        is using "Batch Quality Control (BQC)" samples. Here is how BQC samples are created: take equal amount 
                        of aliquot from each study sample and pool. Then, create multiple BQC samples by transferring the pooled 
                        mixture to new tubes. If the pooled mixture is perfectly homogeneous, the BQC samples will be identical. Therefore, 
                        there should not be any difference between the readings among all BQC samples. In reality, the pooled mixture is not 
                        perfectly homogeneous, therefore a low variance in readings is expected. The vast majority of the lipid species 
                        are expected to have a very
                        low CoV (i.e. CoV < 30%). The coefficient of variation (CoV) is defined as the ratio of the standard deviation to the mean.
                        It shows the extent of variability in relation to the mean of the population and is often expressed as a percentage.
                        The red line in the CoV plot is the line of "30% CoV".
                        
                        When a datapoint (i.e. a lipid species) has a low CoV (i.e. consistent readings among all BQC samples), that datapoint is highly reliable.
                        However, when a datapoint has a high CoV (i.e. inconsistent readings among BQC samples), that datapoint is less reliable. 
                        The high CoV is either the result of the inhomogeneity of the pooled mixture (i.e. valid datapoint) or the error of the instrument
                        and/or the lipid identification software (i.e. invalid datapoint). Not doing any filtering comes with the risk of accepting some 
                        invalid datapoints and filtering comes with the risk of rejecting some valid datapoints. It seems reasonable to have the option to run
                        the analysis with and without BQC filtering and compare the results.   
                        
                        Here is how LipidCruncher filters the data using BQC samples: lipid species with CoV lower than a set threshold (e.g. 30%) are kept
                        and the rest are removed. 
                    
                        """)
                        
                    st.info("""
                            
                            Creating BQC samples in any lipidomics experiment is a great practice:
                                
                            1) The vast majority of the lipid species living in BQC samples are expected to have a low CoV. 
                               Confirming the above statement is a great check for the overall quality of the dataset. 
                            
                            2) BQC samples provide a mechanism to separate datapoints with a very high reliablity from the ones that are less reliable. 
                            
                            """)
                    
                    X = plot_cov(X, cond_lst, rep_lst)
                    
                expand_hist = st.expander("View Distributions of AUC: Scan Data & Detect Atypical Patterns")
                
                with expand_hist:
                    
                    st.markdown("""
    
                        In a standard LipidSearch dataset, columns "MeanArea[s1]" to "MeanArea[sN]" correspond to Area 
                        Under the Curve (AUC) representing the relative abundance of the lipid species 
                        in samples s1 to sN. 
                        
                        To plot the histograms of AUC, LipidCruncher turns all 0 values (i.e. missing values) to 1 
                        and then log-transforms the whole dataset. This allows the user to visualize what portion of the 
                        values are missing without affecting the distribution (as 1 is orders of magnitude smaller than 
                        the minimum detection threshold by mass spectrometry).
                        
                        Visualize the distribution of AUC's in any of the replicates and look for atypical patterns (e.g. too many missing values):
    
                        """)
                
                    plot_hist(X, rep_lst, cond_lst) # histograms of AUC
                    
                expand_retention = st.expander('View Retention Time Plots: Check Sanity of Data')
                
                with expand_retention:
                    
                    st.markdown("""
                        
                        The retention time of a lipid species is a function of its degree of hydrophobicity. 
                        The more hydrophobic the lipid species, the longer the retention time. When retention time is 
                        plotted versus molecular mass, lipid species tend to form separate clusters based upon which lipid class they belong to.
                        
                        Inspect the retention time of lipid species within any lipid class and compare with other lipid classes. 
                        Does everything make sense?
                        
                        """) 
            
                    plot_retention(X, cond_lst, rep_lst) # retention time plots
                    
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
                    
                st.subheader("2.2) Analyze Data")
                
                expand_abundance_bar_chart = st.expander("Class Abundance Bar Chart")
                
                with expand_abundance_bar_chart:
                    
                    st.markdown("""
                                
                                For each condition, the following bar charts show the total abundance of each lipid class in log scale.
                                The total abundance of a class is computed by summing the abundances of all the lipid species belonging 
                                to that class.   
                                
                                """)
                                
                    plot_total_abundance(X, cond_lst, rep_lst)
                    
                expand_vol_plot = st.expander("Volcano Plots - Test Hypothesis")
                
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
                        
                    st.warning("Your choice here only affects the volcano plot, not the rest of the submodule.")
                
                    volcano_plot(X, cond_lst, rep_lst)
                    
                expand_sat_plot = st.expander("Saturation Level Plots - Investigate Saturation Level of Different Lipid Classes")
                
                with expand_sat_plot:
                    
                    st.markdown("""
                                
                                Saturation level plots show the saturation profile of each lipid class. 
                                First, for each lipid species, the ratio of Saturated Fatty Acids (SFA), Mono Unsaturated \
                                Fatty Acids (MUFA) and Poly Unsaturated Fatty Acids (PUFA) is calculated as following:
                                    
                                SFA ratio = total number of saturated fatty acids / total number of fatty acids
                                
                                MUFA ratio = total number of mono unsaturated fatty acids / total number of fatty acids
                                
                                PUFA ratio = total number of poly unsaturated fatty acids / total number of fatty acids
                                
                                Then, for each lipid species, the abundance of SFA, MUFA and PUFA is calculated as following:
                                    
                                AUC(SFA) = (AUC averaged over all replicates).(SFA ratio)
                                
                                AUC(MUFA) = (AUC averaged over all replicates).(MUFA ratio)
                                
                                AUC(PUFA) = (AUC averaged over all replicates).(PUFA ratio)
                                
                                Finally, total AUC(SFA), AUC(MUFA) and AUC(PUFA) for each lipid class is calculated by taking the sum of \
                                AUC(SFA), AUC(MUFA) and AUC(PUFA) over all lipid species that belong to that class. 
                                
                                """)
                    
                    saturation_level_plot(X, cond_lst, rep_lst)
                    
                expand_pathway_plot = st.expander("Lipid Pathway Visualization")
                    
                with expand_pathway_plot:
                        
                    st.markdown("""
                                    
                                    The following is a visualization of lipidomic pathways. Each solid circle corresponds to a lipid class. 
                                    The size of each solid circle is set by the ratio of 
                                    the total abundance of the experimental group to that of the control group. Total abundance of a class is 
                                    computed by summing the abundances of all the lipid species belonging to that class averaged over
                                    all the replicates. 
                                    
                                    Size equal one means 
                                    the total abundance of both control and experimental groups is equal. A black circle with no fill and size 
                                    equal one is plotted on top of each solid circle to make the comparison more convenient. 
                                    
                                    The color of each 
                                    solid circle is set by the ratio of the total saturated fatty acid chains to the total fatty acid chains 
                                    (saturated and unsaturated).       
                                    
                                    """)
                        
                    plot_pathway(X, cond_lst, rep_lst)
