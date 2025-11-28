---
title: "Supporting Information 3"
subtitle: "Disparate estimates of intrinsic productivity for Antarctic krill (Euphausia superba) across small spatial scales, under a rapidly changing ocean."
author: "Mardones, M; Jarvis Mason, E.T.;  Santa Cruz, F.; Watters, G.; Cárdenas, C.A"
date:  "28 November, 2025"
bibliography: LBSPR.bib
csl: apa.csl
link-citations: yes
linkcolor: blue
output:
  bookdown::pdf_document2:
    fig_caption: yes
    keep_md: true
    toc: true
    toc_deep: 3
    toc_float:
      collapsed: false
      smooth_scroll: false
    theme: cosmo
    fontsize: 0.9em
    linestretch: 1.7
    html-math-method: katex
    self-contained: true
    code-tools: true
    number_sections: false
editor_options: 
  markdown: 
    wrap: 72
---





```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```
## Loading required package: gridExtra
```

```
## 
## Attaching package: 'gridExtra'
```

```
## The following object is masked from 'package:dplyr':
## 
##     combine
```
# Methodological Approach: Operating Model Construction  

To evaluate the performance of the Length-Based Spawning Potential Ratio (LBSPR) model in estimating fishing mortality (*F*) and Spawning Potential Ratio (*SPR*) in krill (*Euphausia superba*), we developed an Operating Model (OM) that simulates population dynamics under different fishing pressures. The OM provides a controlled environment to generate synthetic length-frequency data, allowing us to assess the accuracy and bias of LBSPR estimates.  



### Population Dynamics Simulation  
The population was modeled using a size-structured framework based on the von Bertalanffy growth function (VBGF):  

\[
L_t = L_{\infty} \left(1 - e^{-k(t - t_0)}\right)
\]

where:  
- \( L_t \) is the length at age \( t \),  
- \( L_{\infty} \) is the asymptotic maximum length,  
- \( k \) is the growth coefficient,  
- \( t_0 \) is the theoretical age at length zero.  

Natural mortality (\( M \)) was assumed constant and taken from the literature. Recruitment followed a Beverton-Holt stock-recruitment relationship, and selectivity was modeled as a logistic function.  

### Fishing Scenarios  
We simulated multiple fishing mortality scenarios (\( F \)) to assess the robustness of LBSPR estimates:  
- Low exploitation: \( F = 0.2 \)  
- Moderate exploitation: \( F = 0.5 \)  
- High exploitation: \( F = 1.0 \)  
- Overexploitation: \( F = 1.5 \)  

For each scenario, a length-based catch-at-size distribution was generated using a deterministic fishing process combined with stochastic recruitment.  

### Generation of Length-Frequency Data  


The OM generated synthetic length-frequency data by sampling individuals from the simulated population at different time steps. The sampling process mimicked empirical fisheries-dependent and fisheries-independent data collection.   

We applied sampling variability to reflect realistic observation errors and deviations due to gear selectivity.  

### Estimation Procedure Using LBSPR  

The synthetic length data were analyzed using the LBSPR framework (R package *LBSPR*). The model estimated:  
1. Fishing mortality (*F*)  
2. Spawning Potential Ratio (*SPR*)  

The estimated values were compared against true values from the OM to assess bias and accuracy.  

### Performance Evaluation Metrics  

To quantify the performance of LBSPR, we used:  
- Bias: \( Bias = \frac{F_{est} - F_{true}}{F_{true}} \)  
- Mean Absolute Error (MAE): \( MAE = \frac{1}{n} \sum |F_{est} - F_{true}| \)  
- Root Mean Squared Error (RMSE): \( RMSE = \sqrt{\frac{1}{n} \sum (F_{est} - F_{true})^2} \)  

Additionally, we plotted estimated vs. true values to visualize systematic deviations in LBSPR estimates.  

### Interpretation of Results  
- If estimates align with true values, LBSPR is accurate and unbiased.  
- If estimates consistently overestimate or underestimate \( F \) and \( SPR \), LBSPR exhibits systematic bias.  
- Performance under different exploitation levels highlights the robustness of LBSPR in krill stock assessments.  


## MO 




![](MO_LBSPR_files/figure-latex/unnamed-chunk-3-1.pdf)<!-- --> 

  


![](MO_LBSPR_files/figure-latex/unnamed-chunk-4-1.pdf)<!-- --> 







Ahora con L inf


```
## A blank LB_pars object created
```

```
## Default values have been set for some parameters
```

```
## File not found. A blank LB_lengths object created
```

```
## Fitting model
```

```
## Year:
```

```
## 1
```

```
## A blank LB_pars object created
```

```
## Default values have been set for some parameters
```

```
## File not found. A blank LB_lengths object created
```

```
## Fitting model
```

```
## Year:
```

```
## 1
```

```
## A blank LB_pars object created
```

```
## Default values have been set for some parameters
```

```
## File not found. A blank LB_lengths object created
```

```
## Fitting model
```

```
## Year:
```

```
## 1
```

```
## A blank LB_pars object created
```

```
## Default values have been set for some parameters
```

```
## File not found. A blank LB_lengths object created
```

```
## Fitting model
```

```
## Year:
```

```
## 1
```

```
## A blank LB_pars object created
```

```
## Default values have been set for some parameters
```

```
## File not found. A blank LB_lengths object created
```

```
## Fitting model
```

```
## Year:
```

```
## 1
```

```
## A blank LB_pars object created
```

```
## Default values have been set for some parameters
```

```
## File not found. A blank LB_lengths object created
```

```
## Fitting model
```

```
## Year:
```

```
## 1
```

```
## A blank LB_pars object created
```

```
## Default values have been set for some parameters
```

```
## File not found. A blank LB_lengths object created
```

```
## Fitting model
```

```
## Year:
```

```
## 1
```

```
## A blank LB_pars object created
```

```
## Default values have been set for some parameters
```

```
## File not found. A blank LB_lengths object created
```

```
## Fitting model
```

```
## Year:
```

```
## 1
```

```
## A blank LB_pars object created
```

```
## Default values have been set for some parameters
```

```
## File not found. A blank LB_lengths object created
```

```
## Fitting model
```

```
## Year:
```

```
## 1
```

```
## A blank LB_pars object created
```

```
## Default values have been set for some parameters
```

```
## File not found. A blank LB_lengths object created
```

```
## Fitting model
```

```
## Year:
```

```
## 1
```

```
## A blank LB_pars object created
```

```
## Default values have been set for some parameters
```

```
## File not found. A blank LB_lengths object created
```

```
## Fitting model
```

```
## Year:
```

```
## 1
```

```
## A blank LB_pars object created
```

```
## Default values have been set for some parameters
```

```
## File not found. A blank LB_lengths object created
```

```
## Fitting model
```

```
## Year:
```

```
## 1
```

![](MO_LBSPR_files/figure-latex/unnamed-chunk-6-1.pdf)<!-- --> ![](MO_LBSPR_files/figure-latex/unnamed-chunk-6-2.pdf)<!-- --> ![](MO_LBSPR_files/figure-latex/unnamed-chunk-6-3.pdf)<!-- --> 

