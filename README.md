---
title: "Shared methylation regions and multispecies piscine epigenetic clock"
author: "Gabriel Ecker-Eckhofen"
date: "June 2024"
output: html_document
---

### Overview

This repository has been made for a masters thesis project: **Shared Methylation Regions (SMRs): a New Method in Comparative Epigenetics Used to Create a Multispecies Piscine Epigenetic Clock**, which has not been made publicly available yet. 

### Introduction
In this project we introduced a new workflow for identifying comparable methylation data across species. This workflow involves 1) expanding CpG methylation sites from methylation data to sequences (in our case 1kb long) using the species genome, 2) aligning these sequences onto a common reference genome for all species and 3) identifying overlapping aligned sequences which we named **shared methylation regions (SMRs)**.

These SMRs now included CpG sites from all species which 4) we then used for age correlation testing. *Note, that this could be done with any variable you can test correlation for*. 5) SMRs were grouped into two classes. One class contained CpGs which were positively correlating with age and the other was containing negatively correlating ones. 6) Now, we selected only the highest correlating CpG for each species and CpG in every SMR group. This leaves one CpG (for each species) in each SMR. 7) If there was no agreement in correlation direction within an SMR (not all species had CpGs correlating in the same direction), the SMR was excluded. If there were SMRs with positively correlating and negatively correlating CpGs, only the negatively ones were retained. 

This finally allowed us to 8) use the retained SMRs as independent variables for creating a methylation-based age prediction model (also known as an "epigenetic clock"). 9) We tested various models such as multivariate linear regression and non-parametric random forest regression. We achieved notable accuracy in age prediction for four species using a single model ("multispecies epigegentic clock"), the results of which will be published in a paper separate from the thesis. 

### using this repo
As of now, crucial data which has been used in this project is not published yet. It is to be expected that, during the course of this year, all data sets will be made publicly available in separate articles. To reproduce results which have been shown in the master thesis, scripts in the folder *001_scripts* can be looked at and some can be used. The scripts are numerated and everything above "003_...R" can be run with the data made avialble in this repository. Please use the main branch as it is up to data and sufficient for all plots and graphs in the manuscript. 


#### Disclaimer 
- Sequencing data is not published yet and therefore, only scripts involved in the downstream processes can be used as of now.
- This repository is under development and is therefore subject to changes.

for any questions, recommendations or remarks, reach out to me via gabriel.eckhofen@imbresea.eu