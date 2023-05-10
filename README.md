# Project Description

This repository contains Group 5's BF528 Project 1 files. As a result of reproducing the study - "Gene expression classification of colon cancer into molecular subtypes: characterization, validation, and prognostic value", this repository mainly consists of R scripts, output files and a final report.  

For the first project, I was reponsible for the data curation part which invloved navigating the Gene Expression Omnibus (GEO) database and downloading publicly available microarray data to share to my group members (programmer, analyst and biologist) in this pipeline.

# Contributors

Manasa Rupuru - manasarapuru (Analyst)

Pooja Savla - poojas4998 (Data Curator)

Pragya Rawat - rpragya17 (Biologist)

Vrinda Jethalia - vrindajethalia799 (Programmer)

# Repository Contents
BF528_Project1_Report_G5.pdf - Final Report 

project-1-programmer.R - Script written for data preprocessing & quality control. It can be executed using the R console or an integrated development environment (IDE) like RStudio. The programmer also added it as .rmd and knitted the report to html (project-1-programmer.rmd, project-1-programmer.html). All of this is in the programmer_scripts folder.  

Proj_1_Analyst.R - Script written for noise filtering, dimensionality reduction, hierarchical clustering and subtype discovery. It can be executed using the R console or an integrated development environment (IDE) like RStudio.

biologist_script.R - Script written to identify top 10 up- and down-regulated probesets with gene symbol, t-statistic, nominal p-value, and adjusted p-value columns;  and identify top 3 enriched gene sets for each geneset type. It can be executed using the R console or an integrated development environment (IDE) like RStudio.
