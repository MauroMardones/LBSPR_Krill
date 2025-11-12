# **Testing Intrinsic productivity in Krill at Subarea 48.1 (Antarctic Peninsula)**  

## **Overview**  

This project applies a length-based approach to assess the productivity of Antarctic krill (*Euphausia superba*) in Subarea 48.1 (Antarctic Peninsula). Specifically, we use the **Length-Based Spawning Potential Ratio (LBSPR)** methodology [@Hordyk2016] to estimate the sustainability of krill exploitation under data-limited conditions. The study is based on **Yield per Recruit (YPR) models**, adapted for krill fisheries using size-structured population data.  

This research is part of an **academic project** within the **Doctoral Program in Antarctic and Subantarctic Sciences** at the **Universidad de Magallanes, Chile**.  

## **Objectives**  

- Evaluate the reproductive potential of krill using **LBSPR**, a data-limited stock assessment method.  
- Analyze the impact of fishing pressure on krill size distribution and productivity.  
- Provide insights into the sustainability of krill exploitation in Subarea 48.1.  

## **Methodology**  
This repository contains the code and processed data used to reproduce the analyses presented in the manuscript.

The complete code to replicate the exercise is provided in:

**/index.Rmd**

This R Markdown file includes all steps required to run the **LBSPR** model — from data loading and formatting to analysis and figure generation.

### Data
The input data are stored in the **`data/`** folder and include:

- **`lengthJOIN.csv`**  
- **`lengthSSIW.csv`** 
- **`lengthEI.csv`** 
- **`lengthBS.csv`** 
- **`lengthExtra.csv`** 

They are already formatted according to the requirements of the **LBSPR** model and are ready for direct use.

### How to Reproduce the Analysis

1. Open `code/index.Rmd` in **RStudio** (or any compatible R environment).  
2. Set the working directory to the root of this repository.  
3. Knit or run the R Markdown file to reproduce all results and figures.

## **Reproducibility & Code**  
An user friendly script are available in the following repository:  
[🔗 LBSPR Krill Analysis](https://mauromardones.github.io/LBSPR_Krill/)  

## **Contact**  
For questions or collaboration, feel free to reach out.