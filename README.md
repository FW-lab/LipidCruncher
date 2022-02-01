# LipidCruncher  
Mass spectrometry is an analytical technique that is used to measure the mass-to-charge ratio of ions. In the context of lipidomics, the output of mass spectrometry is the relative abundance of the lipid species that make up the sample under study. This output is in the form of a spectrum in which the peaks represent an identified lipid species and the area underneath each peak reperesnts the relative abundance of the corresponding lipid species. There are two pieces of software that turn this spectrum into a lipidomics dataset: LipidSearch and LipidXplorer. LipidCruncher is a web app that allows the user to perform lipidomics analysis on the LipidSearch and LipidXplorer datasets.

Find LipidCruncher at: https://share.streamlit.io/fw-lab/lipidcruncher/main/lipidcruncher.py

Watch a demo on the "Data Exploration & Quality Check Module" of LipidCruncher at: https://www.youtube.com/watch?v=EF406iMQ0qM

## Try The App
The app is simple and intuitive. In order to try the app follow the steps below:
1) Download the sample dataset: search_test1.csv
3) Open the app and click on the "Data Exploration & Quality Check Module" from the navigation bar 
4) Select "LipidSearch" as the type of your dataset 
5) Upload the dataset
6) Define the experiment by enter the inputs on the sidebar:
 
   3 biological conditions: C (control), K (Knock down), O (Over Expression) 
   
   Each condition has 6 corresponding replicates (i.e. 6 samples per condition)
   
6) In the "Apply Built-in Filters" section, set the minimum required MainGrade Score at 6 (i.e. removes the data points that are highly likely to be misidentifications) and skip to the end of the sidebar!
7) Check the box at the end of the sidebar to confirm the inputs
8) Explore the app! 
9) Watch the demo for more info: https://www.youtube.com/watch?v=EF406iMQ0qM
