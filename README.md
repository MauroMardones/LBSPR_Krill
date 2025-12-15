# **Testing Intrinsic Productivity in Antarctic Krill at Subarea 48.1 (Antarctic Peninsula)**

## **Overview**

This project applies a length-based approach to assess the intrinsic productivity of Antarctic krill (*Euphausia superba*) in Subarea 48.1 (Antarctic Peninsula). Specifically, we use the **Length-Based Spawning Potential Ratio (LBSPR)** methodology [@Hordyk2016] to evaluate the sustainability of krill exploitation under data-limited conditions. The analysis is grounded in **Yield per Recruit (YPR)** theory and adapted to krill fisheries using size-structured population data.

This research is part of an **academic project** within the **Doctoral Program in Antarctic and Subantarctic Sciences** at the **Universidad de Magallanes, Chile**.

---

## **Objectives**

* Evaluate the reproductive potential of Antarctic krill using the **LBSPR** framework.
* Analyze the effects of fishing pressure on size structure and population productivity.
* Provide quantitative insights into the sustainability of krill exploitation in Subarea 48.1.

---

## **Methodology**

This repository contains the **source code** required to reproduce the analyses presented in the associated manuscript.

The complete workflow to replicate the exercise is provided in:

**`index.Rmd`**

This R Markdown file includes all steps necessary to run the **LBSPR** model, from data loading and formatting to model execution and figure generation.

---

## **Data Availability and Required Input Files**

The **processed input data** required to run the LBSPR analyses are archived in Zenodo and must be downloaded prior to execution of the code.

### **Zenodo repository**

**DOI:** [https://doi.org/10.5281/zenodo.17936869](https://doi.org/10.5281/zenodo.17936869)

**Citation:**
Mardones, M. (2025). *LBSPR-based stock assessment of Antarctic krill (Euphausia superba): code and processed data* (v1.0.0) [Data set]. Zenodo. [https://doi.org/10.5281/zenodo.17936869](https://doi.org/10.5281/zenodo.17936869)

---

## **How to Reproduce the Analysis**

1. Clone or download this GitHub repository.
2. Download the processed data from Zenodo using the DOI above.
3. Copy the following files into the `data/` directory of the repository:

   * `lenghtBS.csv`
   * `lenghtEI.csv`
   * `lenghtJOIN.csv`
   * `lenghtSSIW.csv`
   * `lenghtExtra.csv`

These files contain **processed length–frequency distributions**, formatted according to the input requirements of the **LBSPR** model.

4. Open `index.Rmd` in **RStudio** (or a compatible R environment).
5. Set the working directory to the root of this repository.
6. Knit or run the R Markdown file to reproduce all results and figures.


> **Note:**
> These datasets are **processed and aggregated** and are provided exclusively to enable reproducibility of the analytical workflow. They do not represent original raw observations and do not substitute primary data sources.

---

## **Reproducibility & Code**

A user-friendly rendered version of the analysis and documentation is available at:

🔗 **LBSPR Krill Analysis**
[https://mauromardones.github.io/LBSPR_Krill/](https://mauromardones.github.io/LBSPR_Krill/)

---

## **Contact**

For questions, comments, or collaboration, feel free to reach out.

