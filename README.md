# LipidCruncher
LipidCruncher is a web app that allows the user to perform lipidomics analysis on LipidSearch and LipidXplorer datasets.  
Find the app at: https://share.streamlit.io/fw-lab/lipidcruncher/main/lipidcruncher.py

Mass spectrometry is an analytical technique that is used to measure the mass-to-charge ratio of ions. The results are presented as a mass spectrum, a plot of intensity as a function of the mass-to-charge ratio. In the context of lipidomics, the output of mass spectrometry is the relative abundance of the lipid species that make up the sample under study. This output is in the form of a spectrum in which the peaks represent an identified lipid species and the area underneath each peak reperesnts the relative abundance of the corresponding lipid species. There are two pieces of software that turn this spectrum into a lipidomics dataset: LipidSearch and LipidXplorer. LipidCruncher is a web app that allows the user to perform lipidomics analysis on the LipidSearch and LipidXplorer datasets.

## Try The App in under 5 minutes!
The app is simple and intuitive. In order to try the app follow the steps below:
1) Download the sample dataset: search_test1.csv
2) Open the app and click on the "Data Exploration & Quality Check Module" from the navigation bar 
3) Select "LipidSearch" as the type of your dataset 
4) Upload the dataset
5) Enter the inputs: 
   3 conditions: WT, KO, Control 
   Each condition has 6 corresponding replicates
6) In the "Apply Built-in Filters" section, set the minimum required MainGrade Score at 6
7) Skip to the bottom of the sidebar and check the box to confirm the inputs
8) Explore the app! 
