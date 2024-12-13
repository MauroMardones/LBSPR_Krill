---
title: "Supporting Information 1"
subtitle: "Spatial and temporal variability in the intrinsic productivity of Antarctic krill (Euphausia superba) along the Western Antarctic Peninsula under environmental and life history scenarios"
author: "Mardones, M; Jarvis Mason, E.T.;  Santa Cruz, F.; Watters, G.; Cárdenas, C.A"
date:  "13 December, 2024"
bibliography: LBSPR.bib
csl: apa.csl
link-citations: yes
linkcolor: blue
output:
  bookdown::html_document2:
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





``` r
## Installing the Package
#The LBSPR package is now available on CRAN:
#install.packages("LBSPR")
#install.packages("devtools")
#devtools::install_github("AdrianHordyk/LBSPR")
###load the package
library(LBSPR)
library(devtools)#to install_github
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(ggpubr)
library(kableExtra)
library(hrbrthemes)
library(ggthemes)
library(tidyverse)
```



# CONTEXT

This supplementary material with essential formulas, datasets, and code snippets is delivered for the application of the LBSPR model in studying krill populations within the West Antarctic Peninsula (WAP) region. The LBSPR model, pivotal for unraveling krill dynamics, is reinforced by mathematical expressions elucidating key parameters and environmental influences. Curated datasets provide valuable insights into factors shaping krill habitats. With a strong emphasis on reproducibility, this resource ensures that findings can undergo independent validation, fostering trust in scientific outcomes. Transparent documentation of methods and assumptions promotes understanding and scrutiny of the analysis process by peers.

\newpage

# METHODOLOGY

## Monitoring Data (SISO Program)

For this analysis, data from the monitoring of the krill fishery were used, which have been systematically collected on board fishing vessels by the CCAMLR SISO (Scheme of International Scientific Observation) program. Krill sizes compositions were obtained from the entire area 48.1, which was combined in each management stratum defined at 2.1 section (Figure \@ref(fig:Figure1)). For this analysis, 1,266,267 krill individuals were measured from fishery activity, the majority (~75%) from the Bransfield, Elephant and SouthWest strata. 

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/Strata_Nochina.png" alt="Sizes compositions from SISO program monitoring krill fishery by strata (BS=Brainsfield Strait, EI= Elephant Island, GS= Gerlache Strait, JOIN= Joinville Island, SSWI= South West). Red line represent recruit size" width="60%" />
<p class="caption">(\#fig:Figure1)Sizes compositions from SISO program monitoring krill fishery by strata (BS=Brainsfield Strait, EI= Elephant Island, GS= Gerlache Strait, JOIN= Joinville Island, SSWI= South West). Red line represent recruit size</p>
</div>

The information gaps (years without sizes composition data) are not calculated because there is no autocorrelation between years, but singular estimators over time.

## Correlation Analysis

Regarding difference between `Spearman` and `pearson`, both, the
Spearman test and the Pearson test are statistical methods used to
assess the correlation between two variables. The main difference
between them is that the Pearson test evaluates the linear correlation
between two continuous variables, while the Spearman test evaluates the
monotonic correlation between two continuous or order variables.

In the case of the Pearson test, the degree of association between two
continuous variables is measured through a correlation coefficient that
varies between -1 and 1. A value of 1 indicates a perfectly positive
correlation, a value of -1 indicates a perfectly negative correlation,
and a value of 0 indicates no correlation between the two variables.

On the other hand, the Spearman test is based on the range of the
variables, instead of the actual values. In other words, this test
evaluates the correlation between two order variables, where the values
of each variable are ranked from lowest to highest, and ranges are used
instead of actual values. Spearman's correlation coefficient also varies
between -1 and 1, but it measures the monotonic correlation between two
variables, that is, if one variable increases, the other variable also
increases or decreases[@McCulloch2001].

In summary, the main difference between the Spearman test and the
Pearson test is that the former is used to assess the monotonic
correlation between two order variables, while the latter is used to
assess the linear correlation between two continuous variables.

-   Spearman's Rank Correlation Coefficient

This coefficient is used to see if there is any significant relationship
between the two datasets, and operates under the assumption that the
data being used is ordinal, which here means that the numbers do not
indicate quantity, but rather they signify a position of place of the
subject's standing (e.g. 1st, 2nd, 3rd, etc.)

$\begin{aligned} r_s = \frac{\sum_{i=1}^n (x_i - \bar{x})(y_i - \bar{y})}{\sqrt{\sum_{i=1}^n (x_i - \bar{x})^2}\sqrt{\sum_{i=1}^n (y_i - \bar{y})^2}} \end{aligned}$

-   Pearson Product-Moment Coefficient

This is the most widely used correlation analysis formula, which
measures the strength of the 'linear' relationships between the raw data
from both variables, rather than their ranks. This is an dimensionless
coefficient, meaning that there are no data-related boundaries to be
considered when conducting analyses with this formula, which is a reason
why this coefficient is the first formula researchers try.

$\begin{aligned} r = 1- \frac{6\sum_{i=1}^n D_{i}^n}{n (n^2 - 1)}\end{aligned}$


We can identify through a correlation matrix the data of our set, whether it is positive or negative. The main outputs are displayed in Figure \@ref(fig:Figure2). The comprehensive analysis with methodological routine can be follow in author Github repository [Lenght-Enviromental_Krill](https://github.com/MauroMardones/Krill_Length_Cor.git)

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/corr_length.png" alt="Pearson correlation between krill length and temporal (years) and environmental factors in the WAP" width="80%" />
<p class="caption">(\#fig:Figure2)Pearson correlation between krill length and temporal (years) and environmental factors in the WAP</p>
</div>
## GLMM models

Random effects are a way to model variability in data that comes from
factors that cannot be directly measured or controlled. In the context
of statistical models, random effects refer to variables that are
assumed to have an unknown probability distribution, and are included in
the model to explain some of the variation in the data.

For example, in a study comparing the test scores of students at
different schools, random effects refer to differences between schools
that cannot be explained by variables measured in the study, such as
students' socioeconomic status. These differences may be caused by
factors such as the quality of teaching, school culture, or geographic
location.

Random effects are often modeled by using mixed effects models, which
combine random and fixed effects in the same model. Fixed effects are
those that are assumed to be constant for all study units and are
directly measured, while random effects are those that are assumed to
vary randomly across study units and cannot be directly measured.

In short, random effects are a way of modeling variability in data that
cannot be directly explained by the variables measured in the study, and
are included in the model to improve the precision of the estimates and
reduce the potential for bias [@McCulloch2001].

In this analysis we try test spatial componentent in `cellid` variable like
random effects with `lme4` package [@Bates2015] 

We fitted two linear mixed-effects models to evaluate the influence of environmental variables on krill length:

1. Model 1:  
\[
\text{Length} = \beta_0 + \beta_1 \cdot \text{Year} + \beta_2 \cdot \text{Chl} + \beta_3 \cdot \text{SIC} + \beta_4 \cdot \text{SST} + \text{Random Effect}_{\text{cellid}}
\]

2. Model 2:  
\[
\text{Length} = \beta_0 + \beta_1 \cdot \text{Year} + \beta_2 \cdot \text{Chl} + \beta_3 \cdot \text{SIC} + \beta_4 \cdot \text{SST} + \text{Random Effect}_{\text{Year}}
\]

Here, \(\beta_0\) is the intercept, \(\beta_1, \beta_2, \beta_3, \beta_4\) are the fixed-effect coefficients for the variables *Year*, *Chl*, *SIC*, and *SST*, respectively, and the random effects are specified at the *cellid* and *Year* levels in the respective models.

Main outpus displayed in Figure \@ref(fig:Figure3).

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/glmm.png" alt="Random effects of location (cell) and fixed effects of environmental variables on krill length, with scattered points indicating variability in impact in estimates of the main variable." width="80%" />
<p class="caption">(\#fig:Figure3)Random effects of location (cell) and fixed effects of environmental variables on krill length, with scattered points indicating variability in impact in estimates of the main variable.</p>
</div>
## LBSPR Modeling

### Biological and fishery parameters krill


``` r
#Create a new LB_pars Object
#To create a new LB_pars object you use the new function:
MyPars <- new("LB_pars")

#You can see the elements or slots of the LB_pars object using the
#slotNames function:
#slotNames(MyPars)
```

The model needs specifications related to both biological and fishery parameters according to species evaluated. In a descriptive way, the main parameter sets are described as follows;

**Biology**

-   von Bertalanffy asymptotic length `Linf`
-   M/K ratio (natural mortality)divided by von Bertalanffy K coefficient) `MK`
-   Length at 50% maturity (`L50`)
-   Length at 95% maturity (`L95`)

**Fishery**

-   Length at 50% selectivity (`SL50`)
-   Length at 95% selectivity (`SL95`)
-   Biological Reference Point (BRP). F/M ratio (`FM`) or Spawning Potential Ratio (`SPR`). If you specify both, the F/M value will be ignored.

**Size Classes**

-Width of the length classes (`BinWidth`)


``` r
#MyPars <- new("LB_pars")
## A blank LB_pars object created
## Default values have been set for some parameters
MyPars@Species <- "Euphausia superba"
MyPars@Linf <- 60 
MyPars@L50 <- 34 
MyPars@L95 <- 55 # verrificar bibliografia
MyPars@MK <- 0.4/0.45


#Explotation
MyPars@SL50 <- 40#numeric() #1
MyPars@SL95 <- 56#numeric() #27
MyPars@SPR <- 0.75 #numeric()# ###cambia el numero 0.4 a en blanco
MyPars@BinWidth <- 1
#MyPars@FM <- 1

MyPars@Walpha <- 1
MyPars@Wbeta <- 3.0637 #r2 = 0.9651

MyPars@BinWidth <-1
MyPars@BinMax <- 70
MyPars@BinMin <- 0
MyPars@L_units <- "mm"
```



``` r
tablepar <-data.frame(Value=c(MyPars@Linf,
                       MyPars@L50 ,
                       MyPars@L95 ,
                       round(MyPars@MK,3) ,
                       MyPars@SL50 ,
                       MyPars@SL95 ,
                       MyPars@SPR ,
                       MyPars@Walpha ,
                       MyPars@Wbeta ,
                       MyPars@BinMax,
                       MyPars@BinMin ,
                       MyPars@L_units),
                      Descrption=c("VB asymptotic length",
                                   "Maturity 50%",
                                   "Maturity 95%",
                                   "M/K Ratio",
                                   "Selectivity 50%",
                                   "Seletivity 95%",
                                   "SPR",
                                   "a (Length-Weight Relation)",
                                   "b (Length-Weight Relation)",
                                   "Bin Min",
                                   "Bin Max",
                                   "Units"))
```

These parameters were taken from previous works about krill life history and fishery [@Thanassekos2014; @Maschette2020], which are described in Table \@ref(tab:Table1).


``` r
kbl(tablepar, 
    longtable = F, 
    booktabs = T, 
    caption = "Krill biological and fishery parameters") %>% 
   kable_styling(latex_options = c("striped", 
                                   "hold_position"),
                 html_font = "arial") 
```

<table class="table" style="color: black; font-family: arial; margin-left: auto; margin-right: auto;">
<caption>(\#tab:Table1)Krill biological and fishery parameters</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Value </th>
   <th style="text-align:left;"> Descrption </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 60 </td>
   <td style="text-align:left;"> VB asymptotic length </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 34 </td>
   <td style="text-align:left;"> Maturity 50% </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 55 </td>
   <td style="text-align:left;"> Maturity 95% </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 0.889 </td>
   <td style="text-align:left;"> M/K Ratio </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 40 </td>
   <td style="text-align:left;"> Selectivity 50% </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 56 </td>
   <td style="text-align:left;"> Seletivity 95% </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 0.75 </td>
   <td style="text-align:left;"> SPR </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> a (Length-Weight Relation) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 3.0637 </td>
   <td style="text-align:left;"> b (Length-Weight Relation) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 70 </td>
   <td style="text-align:left;"> Bin Min </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> Bin Max </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mm </td>
   <td style="text-align:left;"> Units </td>
  </tr>
</tbody>
</table>


### Model Estimation LBSPR

Recent work has shown that under equilibrium conditions (that is, constant F and no recruitment variability) and assuming the von Bertalanffy growth equation, constant natural mortality for all ages, and logistic or jack-knife selectivity, standardization of the composition of lengths of two populations with the same ratio of natural mortality to growth rate (*M/k*) and the same ratio of mortality by fishing to natural mortality (*F/M*) will be identical [@Hordyk2016]. Extension of this model to incorporate length-at-age variability and logistic selectivity confirms that, at equilibrium, the composition of the predicted duration of catch of an exploited population is primarily determined by the ratios of M/k and F/M. The analytical models developed in @Hordyk2014 suggest that with knowledge of the asymptotic von Bertalanffy length $L_{\infty}$ and the coefficient of variation in $CVL_{\infty}$, the ratio of total mortality to the von Bertalanffy growth coefficient (*Z/k*) for a given population can be estimated from a representative sample of the size structure of the catch. If *M/k* (or parameters) is also known, then the results of @Hordyk2016 suggest that it is possible to estimate F/M from the composition of the catch. Often the F/M ratio has been used as a biological reference point when is 1 [@Zhou2012].

The LBSPR model requires the following parameters: an estimate of the M/k ratio, $L_{\infty}$, $CVL_{\infty}$, and knowledge of maturity by length (maturity ogive), both set parameters know in krill. This model uses data on composition by catch sizes to estimate intrinsic productivity or Spawning Potential Ratio (SPR). This concept was extracted from @Goodyear1993, where the ratio of lifetime average egg production per recruit (EPR) was calculated for fished and no fished (virgin condition) resources. An algorithm route for the calculation of the SPR is the following;

$${SPR}=\frac{EPR_{fished}}{EPR_{nofished}}$$ where;

$$
EPR_{fished} =
\sum_{a}
\begin{cases}
 E_{a},  a = 0  \\
 e^{-Z_{{a-1}}a} {E_a}, 0 < a  < a_\le{max}
 \end{cases}       
$$ 

where $Z_a$ = $M+F_a$, and $E$ is egg production at age assumed to be proportional to weight;

$$E_a \in  Mat_a W_a$$ on the other hand, the calculation of the reproductive potential is the same as that of those captured without F; $$
EPR_{nofished} =
\sum_{a}
 E_ae^{-M_a}
$$ Assuming that egg production is proportional to the size of mature fish, relative fecundity-at-size is given by;

$$
Fec_{L,g} = Mat_{L,g} L^b
$$ where b is value of the exponent to reflect differente size fecunditity relationship and g is the fraction of recruits to group.

Assuming reasonable estimates of the M/K ratio, $L_{\infty}$ (or $CVL_{\infty}$), size-at-maturity, the parameters F/M, $SL_{50}$, and $SL_{95}$ can be estimated from a representative sample of the length structure of the catch by minimizing the following multinomial negative loglikelihood function (NLL):

$$
NLL =
argmin\sum_{i}
 O_i ln\frac{\hat{P}_i}{\hat{O}_i}
$$

where $O_i$ and $\hat{O}_i$ are the observed number and proportion in length class i, respectively, and $\hat{P}_i$ is the model estimate of the probability in length class i [@Hordyk2016].


The LBSPR package can be used to generate the expected size composition, the SPR, and relative yield for a given set of biological and exploitation pattern parameters. The output of the `LBSPRsim` function is an object of class `LB_obj`. This is another LBSPR object, and contains all of the information from the `LB_pars` object and the output of the `LBSPRsim` function.

It is important to note that the F/M ratio reported in the LBSPR model refers to the apical F over the adult natural mortality rate. That is, the value for fishing mortality refers to the highest level of F experienced by any single size class [@Hordyk2014]. If the selectivity pattern excludes all but the largest individuals from being exploited, it is possible to have a very high F/M ratio in a sustainable fishery (high SPR) and visceverse. Note that only the life history parameters need to be specified for the estimation model. The exploitation parameters will be estimated [@Hordyk2014].

Two objects are required to fit the LBSPR model to length data: life-history parameters (`LB_pars`) described previously (2.4. section) and size compositions (`LB_lengths`), which contains the length frequency data. We provide a set of global size data and also by strata, with which we will be able to identify the spatial differences of the potential.


``` r
MyLengths <- new("LB_lengths")
slotNames(MyLengths)
```

```
## [1] "LMids"   "LData"   "L_units" "Years"   "NYears"  "Elog"
```

However, it is probably easier to create the `LB_lengths` object by directly reading in a `.csv`. A length frequency data of krill set with multiple years (2001-2020) by strata in 48.1 Subarea.


``` r
#Brainsflied Strata
datdir <- setwd("~/DOCAS/LBSPR_Krill")
Lenbs <- new("LB_lengths", 
             LB_pars=MyPars, 
             file=paste0(datdir, "/lenghtBS_nochina.csv"), 
             dataType="freq",
             sep=",",
             header=T)
#Elephan Island Strata
Lenei <- new("LB_lengths", 
             LB_pars=MyPars, 
             file=paste0(datdir, "/lenghtEI_nochina.csv"),
             dataType="freq",
             sep=",",
             header=T)
#Gerlache Strata
Lengc <- new("LB_lengths", 
             LB_pars=MyPars, 
             file=paste0(datdir, "/lenghtExtra_nochina.csv"), 
             dataType="freq",
             sep=",",
             header=T)
#Join Strata
Lenjo <- new("LB_lengths", 
             LB_pars=MyPars, 
             file=paste0(datdir, "/lenghtJOIN_nochina.csv"), 
             dataType="freq",
             sep=",",
             header=T)
#SSIW Strata
Lenssiw <- new("LB_lengths", 
               LB_pars=MyPars, 
               file=paste0(datdir, "/lenghtSSIW_nochina.csv"), 
               dataType="freq",
               sep=",",
               header=T)
```


### Fit the Model by strata

The LBSPR model is fitted by strata using the `LBSPRfit` function. 


``` r
fitbs <- LBSPRfit(MyPars, Lenbs)
fitei <- LBSPRfit(MyPars, Lenei)
fitgc <- LBSPRfit(MyPars, Lengc)
fitjo <- LBSPRfit(MyPars, Lenjo)
fitssiw <- LBSPRfit(MyPars, Lenssiw)
```



### Sensitivity and perfomance analysis metohodology.

Ten sensitivity scenarios based on the upper and lower range for the asymptotic length von Bertalanffy $L_{\infty}$ (55 to 65 mm) used in the base model (60 mm) were tested to identify the impact of this parameter on the SPR estimation, given that it is the parameters exhibiting the high degree of variability and one of the factors that most determines the estimates. On the other hand, the interdependence between krill and their environment is a well-known and influential factor in population dynamics, ecosystem impacts, and fishery. This interdependence also affects the reproductive potential and consequently, to any management decision that takes into account population parameters of krill. To assess the impact of individual growth variability, three scenarios of the VB k parameter were tested, representing different growth types (low = 0.2, medium = 0.7, and high = 1.2). Figure \@ref(fig:Figure4) displays a theoretical growth curve for krill based on three scenarios that were tested using LBSPR.



``` r
# Definir los parámetros del modelo de Von Bertalanffy
L_inf <- 60
k <- c(0.2, 0.7, 1.2)
t_0 <- 0
# Crear un vector de edades
edades <- seq(0, 8, by = 0.2)
# Calcular los valores teóricos de crecimiento somático
df <- data.frame(edades = rep(edades, 
                              length(k)), 
                 k = rep(k, each = length(edades)))

df <- df %>% 
  mutate(crecimiento = L_inf * (1 - exp(-k * (edades - t_0))))
# Graficar la curva de crecimiento somático
sen <- ggplot(df, aes(x = edades, 
               y = crecimiento, 
               color = factor(k))) +
  geom_line(size=1.5) +
  xlab("Ages") +
  ylab("Length (mm)") +
  ggtitle("") +
  theme(panel.background = element_rect(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_color_viridis_d(option="C",
                       name = "Theoretical growth",
                       labels = c("Low", "Medium", "High"))+
  theme_few()
sen
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/Figure4-1.jpeg" alt="Theoretical growth curves for krill based on three SPR sensitivity scenarios that were tested using LBSPR." width="70%" />
<p class="caption">(\#fig:Figure4)Theoretical growth curves for krill based on three SPR sensitivity scenarios that were tested using LBSPR.</p>
</div>

After applying each scenario using each of the parameter setting, the results of scenarios are compared with the results provided by the methods based setting (Table 1), analyzing in this way the efect of underestimation/overestimation of the parameters $L_{\infty}$ and *k*.


# RESULTS

## Model Perfomance

The results demonstrate good fits of the size structure distribution model for krill across the strata. The model accurately captures the distribution patterns of size classes, indicating its effectiveness in characterizing the population structure. The model successfully captures the variations in size distribution, reflecting the natural variability in krill populations across different strata. These findings suggest that the model can be utilized as a valuable tool for understanding and predicting size structure dynamics in krill populations (Figure \@ref(fig:Figure5)).


``` r
bscom1 <-plotSize(fitbs,
                  scales = "fixed",
         Title="Bransfield strata")+
  thememau
eicom1 <-plotSize(fitei,
         Title="Elephan Island strata")+
  thememau
gscom1 <- plotSize(fitgc,
         Title="Gerlache strata")+
  thememau
jocom1 <- plotSize(fitjo,
         Title="Joinville strata")+
  thememau
sswicom1 <- plotSize(fitssiw,
         Title="SSWI strata")+
  thememau


ggarrange(bscom1,
          eicom1, 
          gscom1,
          jocom1, sswicom1, 
          ncol = 1,
          legend="bottom",
          common.legend = TRUE)
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/Figure5-1.jpeg" alt="Fit of the model to the data of lengths in all strata"  />
<p class="caption">(\#fig:Figure5)Fit of the model to the data of lengths in all strata</p>
</div>


The difference between the observed accumulated size compositions for each stratum and compare it with the expected size composition at a target SPR (75% SPR). In the simulation of the structure in its virgin condition (without fishing), the red bars represent each stratum. Additionally, the overlap with the average structures observed during the years of fishery monitoring can be visualized. The SSWI Gerlache and EI strata exhibit the greatest differences from the simulated structure, possibly due to the significant contribution of juveniles in these strata. Conversely, the BS and JO strata demonstrate the closest resemblance to the simulated structure (Figure \@ref(fig:Figure6)).


``` r
yr <- 5 # first year of data

bscom <- plotTarg(MyPars, Lenbs, yr=yr,
                  Cols = c(2,4),
                  size.axtex = 8,
                  title="BS",
                  targtext = FALSE)+
  thememau 
eicom <- plotTarg(MyPars, Lenei, yr=yr,
                  Cols = c(2,4),
                  size.axtex = 8,
                  title="EI",
                  targtext = FALSE)+
  thememau 
gccom <- plotTarg(MyPars, Lengc, yr=yr,
                  Cols = c(2,4),
                  size.axtex = 8,
                  title="GS",
                  targtext = FALSE)+
  thememau 
jocom <- plotTarg(MyPars, Lenjo, yr=1,
                  Cols = c(2,4),
                  size.axtex = 8,
                  title="JOIN",
                  targtext = FALSE)+
  thememau 
sscom <- plotTarg(MyPars, Lenssiw , yr=yr,
                  Cols = c(2,4),
                  size.axtex = 8,
                  title="SSWI",
                  targtext = FALSE)+
  thememau 

ggarrange(bscom, eicom + rremove("ylab"), 
          gccom + rremove("ylab"),
          jocom, sscom + rremove("ylab"), 
          ncol = 3, nrow = 2,
          legend="bottom",
          common.legend = TRUE)
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/Figure6-1.jpeg" alt="Difference between the observed accumulated size structure for each stratum related SPR objective" width="80%" />
<p class="caption">(\#fig:Figure6)Difference between the observed accumulated size structure for each stratum related SPR objective</p>
</div>

Furthermore, ogive maturity curve specific to each stratum and year, as well as the estimated length selectivity curve, are presented in (Figure \@ref(fig:Figure7)). These curves offer crucial insights into the fishery's impact on the population and its reproductive status, providing measures of the population's vulnerability to fishing mortality. Notably, the Brainsfield stratum stands out with a lower proportion of mature individuals, suggesting a higher prevalence of juveniles. It is important to note that the same maturity parameters were applied across all strata.


``` r
bsmat <- plotMat(fitbs,
        useSmooth = TRUE,
        Title="BS")+
  thememau 
eimat <-  plotMat(fitei,
        useSmooth = TRUE,
        Title="EI")+
  thememau 
gcmat <-  plotMat(fitgc,
        useSmooth = TRUE,
        Title="GS")+
  thememau 
jomat <-  plotMat(fitjo,
        useSmooth = TRUE,
        Title="JOIN")+
  thememau 
ssmat <-  plotMat(fitssiw,
        useSmooth = TRUE,
        Title="SSWI")+
  thememau 

ggarrange(bsmat,
          eimat + rremove("ylab"), 
          gcmat , 
          jomat + rremove("ylab"), 
          ssmat, 
          ncol = 2, nrow = 3,
          legend="right",
          common.legend = TRUE)
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/Figure7-1.jpeg" alt="Maturity curves by strata" width="80%" />
<p class="caption">(\#fig:Figure7)Maturity curves by strata</p>
</div>

## Comparing producivity between years and Stratas

The analysis of the krill population's reproductive potential across different years and strata reveals significant differences. Brainsfield and Gerlache strata exhibit a low reproductive potential below the proposed management target of 75% in the last year, with values of 0.121 and 0.085, respectively, falling even below the limit reference point. This condition arises from the concentration of a substantial number of immature individuals (juveniles) in these strata, which are being exploited by the fishery, thereby hindering their reproductive cycles from completing. On the contrary, the Elephant Island stratum demonstrates higher spawning potential ratio (SPR) levels in recent years, reaching 0.421 in 2019, which aligns closer to the management objective. This discrepancy is attributed to the spatial distribution of krill, as the Elephant Island stratum possesses a larger proportion of adult individuals compared to other strata.


``` r
sprbs <- as.data.frame(cbind(fitbs@Years, 
                             fitbs@SPR,
                             fitbs@Vars[,4])) 
colnames(sprbs) <- c("Year","SPR", "Var")
sprbs$SPRv <- rep("BS", nrow(sprbs))

sprei <- as.data.frame(cbind(fitei@Years, 
                             fitei@SPR,
                             fitei@Vars[,4])) 
colnames(sprei) <- c("Year","SPR", "Var")
sprei$SPRv <- rep("EI", nrow(sprei))

sprgc <- as.data.frame(cbind(fitgc@Years, 
                             fitgc@SPR,
                             fitgc@Vars[,4])) 
colnames(sprgc) <- c("Year","SPR", "Var")
sprgc$SPRv <- rep("GS", nrow(sprgc))

sprjo <- as.data.frame(cbind(fitjo@Years, 
                             fitjo@SPR,
                             fitjo@Vars[,4])) 
colnames(sprjo) <- c("Year","SPR", "Var")
sprjo$SPRv <- rep("JOIN", nrow(sprjo))

sprssiw <- as.data.frame(cbind(fitssiw@Years, 
                             fitssiw@SPR,
                             fitssiw@Vars[,4])) 
colnames(sprssiw) <- c("Year","SPR","Var")
sprssiw$SPRv <- rep("SSWI", nrow(sprssiw))

allspr <- rbind(sprbs, sprei, sprgc, sprjo, sprssiw)

allspr <- allspr %>%
  mutate(SD = sqrt(Var))


allsprwide <- round(pivot_wider(allspr,
                          names_from = SPRv,
                          values_from = c(SPR, Var, SD),
                          values_fill = NA,
                          names_sep = " ",),3)
```

Figure \@ref(fig:Figure8) provides a visual representation of the SPR trends across years and strata, clearly indicating the references (yellow line = 75% SPR Objective and Red line = 20% Limit SPR).


``` r
allsprpl2 <- ggplot(allspr, aes(Year, SPR)) +
  geom_point(alpha = 0.8) +  # Puntos de datos
  geom_errorbar(aes(ymin = SPR - SD, ymax = SPR + SD),  # Barras de error
                width = 0.2,  # Ancho de la barra horizontal
                alpha = 0.5) +
  geom_smooth(method = "lm",
              alpha=0.3,
              se=TRUE,
              col="black",
              level = 0.75)+
  geom_hline(yintercept = 0.75,
             colour = '#006d2c',
             alpha = 0.5,
             linewidth = 1.2,
             linetype = 2) +
  geom_hline(yintercept = 0.20,
             colour = '#bd0026',
             alpha = 0.5,
             linewidth = 1.2) +
  facet_wrap(. ~ SPRv, ncol = 5) +
  xlim(2001, 2021) +
  ylim(0, 0.9) +
  theme_few() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

allsprpl2
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/Figure8-1.jpeg" alt="Krill Intrinsic Productivity (SPR) by strata and by year"  />
<p class="caption">(\#fig:Figure8)Krill Intrinsic Productivity (SPR) by strata and by year</p>
</div>

All the estimated values of SPR and their associated variance by stratum and by year can be identified in  Table \@ref(tab:Table2).


``` r
allsprwide_reorganizado <- allsprwide %>%
  mutate(
    `BS` = ifelse(!is.na(`SPR BS`), paste0(`SPR BS`, " (", `SD BS`, ")"), NA),
    `EI` = ifelse(!is.na(`SPR EI`), paste0(`SPR EI`, " (", `SD EI`, ")"), NA),
    `GS` = ifelse(!is.na(`SPR GS`), paste0(`SPR GS`, " (", `SD GS`, ")"), NA),
    `JOIN` = ifelse(!is.na(`SPR JOIN`), paste0(`SPR JOIN`, " (", `SD JOIN`, ")"), NA),
    `SSWI` = ifelse(!is.na(`SPR SSWI`), paste0(`SPR SSWI`, " (", `SD SSWI`, ")"), NA)
  ) %>%
  select(Year, `BS`, `EI`, `GS`, `JOIN`, `SSWI`) %>% 
  arrange(Year)

# guardo
#write_csv(allsprwide_reorganizado, "SPR_Total.csv")

kbl(allsprwide_reorganizado, 
    longtable = FALSE, 
    booktabs = TRUE, 
    caption = "Estimates of SPR by Strata (values in parentheses represent the standard deviation)") %>% 
    kable_styling(latex_options = c("scale_down", "hold_position"))
```

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:Table2)Estimates of SPR by Strata (values in parentheses represent the standard deviation)</caption>
 <thead>
  <tr>
   <th style="text-align:right;"> Year </th>
   <th style="text-align:left;"> BS </th>
   <th style="text-align:left;"> EI </th>
   <th style="text-align:left;"> GS </th>
   <th style="text-align:left;"> JOIN </th>
   <th style="text-align:left;"> SSWI </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 2001 </td>
   <td style="text-align:left;"> 0.132 (0.024) </td>
   <td style="text-align:left;"> 0.223 (0.062) </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.196 (0.011) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2004 </td>
   <td style="text-align:left;"> 0.15 (0.008) </td>
   <td style="text-align:left;"> 0.219 (0.028) </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.21 (0.045) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2005 </td>
   <td style="text-align:left;"> 0.25 (0.014) </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.253 (0.005) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2006 </td>
   <td style="text-align:left;"> 0.214 (0.009) </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.152 (0.09) </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.349 (0.025) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2007 </td>
   <td style="text-align:left;"> 0.108 (0.078) </td>
   <td style="text-align:left;"> 0.095 (0.008) </td>
   <td style="text-align:left;"> 0.069 (0.014) </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.277 (0.057) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2008 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.069 (0.008) </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2009 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.322 (0.058) </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.232 (0.016) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2010 </td>
   <td style="text-align:left;"> 0.197 (0.004) </td>
   <td style="text-align:left;"> 0.329 (0.02) </td>
   <td style="text-align:left;"> 0.242 (0.014) </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.31 (0.008) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2011 </td>
   <td style="text-align:left;"> 0.185 (0.012) </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.426 (0.009) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2012 </td>
   <td style="text-align:left;"> 0.156 (0.004) </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.148 (0.003) </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.386 (0.006) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2013 </td>
   <td style="text-align:left;"> 0.13 (0.001) </td>
   <td style="text-align:left;"> 0.191 (0.013) </td>
   <td style="text-align:left;"> 0.123 (0.002) </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.176 (0.005) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2014 </td>
   <td style="text-align:left;"> 0.175 (0.002) </td>
   <td style="text-align:left;"> 0.201 (0.007) </td>
   <td style="text-align:left;"> 0.172 (0.006) </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.311 (0.004) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2015 </td>
   <td style="text-align:left;"> 0.19 (0.003) </td>
   <td style="text-align:left;"> 0.389 (0.022) </td>
   <td style="text-align:left;"> 0.177 (0.003) </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.386 (0.026) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2016 </td>
   <td style="text-align:left;"> 0.288 (0.01) </td>
   <td style="text-align:left;"> 0.36 (0.057) </td>
   <td style="text-align:left;"> 0.293 (0.007) </td>
   <td style="text-align:left;"> 0.331 (0.19) </td>
   <td style="text-align:left;"> 0.266 (0.033) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2017 </td>
   <td style="text-align:left;"> 0.165 (0.004) </td>
   <td style="text-align:left;"> 1 (0) </td>
   <td style="text-align:left;"> 0.214 (0.016) </td>
   <td style="text-align:left;"> 0.249 (0.022) </td>
   <td style="text-align:left;"> 0.27 (0.01) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2018 </td>
   <td style="text-align:left;"> 0.165 (0.003) </td>
   <td style="text-align:left;"> 0.334 (0.019) </td>
   <td style="text-align:left;"> 0.162 (0.004) </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.294 (0.018) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2019 </td>
   <td style="text-align:left;"> 0.161 (0.004) </td>
   <td style="text-align:left;"> 0.425 (0.02) </td>
   <td style="text-align:left;"> 0.141 (0.005) </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.349 (0.006) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2020 </td>
   <td style="text-align:left;"> 0.11 (0.004) </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.087 (0.004) </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.24 (0.012) </td>
  </tr>
</tbody>
</table>


\newpage

## Sensitivity and perfomance analysis

The results of the methods in the reference setting are compared to the obtained under overstimation/underestimation in intrinsic productivity regarding asymptotic length von Bertalanffy $L_{\infty}$ parameter in krill. First, it was possible to identify that for all strata, the impact of low $L_{\infty}$ ranges under the base model (60 mm) overestimates the level of reproductive potential krill with values between 42% and 32% (Bransfield and Elephant Island strata respectively). Regarding higher $L_{\infty}$ settings, the model tends to underestimate the reproductive potential with values between -25% and -30% (Gerlache and Joinville strata) Figure \@ref(fig:Figure9),  Table \@ref(tab:Table3).


``` r
for (i in 55:65) {
  nombre_objeto <- paste0("MyPars_", i)  # Crear un nombre único para cada objeto modificado
  objeto_modificado <- MyPars  # Crear una copia del objeto original
  
  objeto_modificado@Linf <- i  # Modificar el valor de @Linf en cada iteración
  
  assign(nombre_objeto, objeto_modificado, envir = .GlobalEnv)  # Asignar el objeto modificado al nombre creado
  
  # Resto del código que deseas ejecutar con el objeto modificado
## A blank LB_pars object created
## Default values have been set for some parameters
      MyPars@Species <- "Euphausia superba"
      MyPars@L50 <- 34 
      MyPars@L95 <- 55 # verrificar bibliografia
      MyPars@MK <- 0.4/0.45
      
      
      #Explotation
      MyPars@SL50 <- 40#numeric() #1
      MyPars@SL95 <- 56#numeric() #27
      MyPars@SPR <- 0.75 #numeric()# ###cambia el numero 0.4 a en blanco
      MyPars@BinWidth <- 1
      #MyPars@FM <- 1
      
      MyPars@Walpha <- 1
      MyPars@Wbeta <- 3.0637 #r2 = 0.9651
      
      MyPars@BinWidth <-1
      MyPars@BinMax <- 70
      MyPars@BinMin <- 0
      MyPars@L_units <- "mm"
  
  # Imprimir el valor de @Linf después de cada cambio
  #print(get(nombre_objeto)@Linf)
}
```




``` r
# ciclo de ajuste para BS
for (i in 55:65) {
  nombre_objeto <- paste0("MyPars_", i)  # Nombre del objeto S4 en cada iteración
  nombre_variable <- paste0("fitbs", i)  # Nuevo nombre para el resultado
  
  # Obtener el objeto S4 correspondiente
  objeto <- get(nombre_objeto)
  
  # Ejecutar LBSPRfit() y asignar el resultado a una nueva variable con nombre único
  assign(nombre_variable, LBSPRfit(objeto, Lenbs), envir = .GlobalEnv)
}


# ciclo de ajuste para Gerlache
for (i in 55:65) {
  nombre_objeto <- paste0("MyPars_", i)  
  nombre_variable <- paste0("fitei", i)  
  objeto <- get(nombre_objeto)
  assign(nombre_variable, LBSPRfit(objeto, Lenei), envir = .GlobalEnv)
}

# ciclo de ajuste para Gerlache
for (i in 55:65) {
  nombre_objeto <- paste0("MyPars_", i)  
  nombre_variable <- paste0("fitgc", i)  
  objeto <- get(nombre_objeto)
  assign(nombre_variable, LBSPRfit(objeto, Lengc), envir = .GlobalEnv)
}

# ciclo de ajuste para JoinVIlle
for (i in 55:65) {
  nombre_objeto <- paste0("MyPars_", i) 
  nombre_variable <- paste0("fitjo", i)  
  objeto <- get(nombre_objeto)
  assign(nombre_variable, LBSPRfit(objeto, Lenjo), envir = .GlobalEnv)
}

# ciclo de ajuste para JSSWI
for (i in 55:65) {
  nombre_objeto <- paste0("MyPars_", i) 
  nombre_variable <- paste0("fitsswi", i) 
  objeto <- get(nombre_objeto)
  assign(nombre_variable, LBSPRfit(objeto, Lenssiw), envir = .GlobalEnv)
}
```




``` r
#Genero lo que quiero comparar BS
valbs <- as.data.frame(cbind(fitbs55@Years,
  fitbs55@SPR,
           fitbs56@SPR,
           fitbs57@SPR,
           fitbs58@SPR,
           fitbs59@SPR,
           fitbs60@SPR,
           fitbs61@SPR,
           fitbs62@SPR,
           fitbs63@SPR,
           fitbs64@SPR,
           fitbs65@SPR))
colnames(valbs) <- c("Years", "Linf55","Linf56", "Linf57", "Linf58", 
                     "Linf59", "Linf60", 
                      "Linf61", "Linf62" , "Linf63",  "Linf64",
                     "Linf65")

valbs_largo <- pivot_longer(valbs, 
                            cols = c(2:12,), 
                            names_to = "Parameter", 
                            values_to = "SPR")
valbs_largo$SPRstra <- rep("BS", nrow(valbs_largo))

#Genero lo que quiero comparar EI
valei <- as.data.frame(cbind(fitei55@Years,
  fitei55@SPR,
           fitei56@SPR,
           fitei57@SPR,
           fitei58@SPR,
           fitei59@SPR,
           fitei60@SPR,
           fitei61@SPR,
           fitei62@SPR,
           fitei63@SPR,
           fitei64@SPR,
           fitei65@SPR))
colnames(valei) <- c("Years", "Linf55","Linf56", "Linf57", "Linf58", 
                     "Linf59", "Linf60", 
                      "Linf61", "Linf62" , "Linf63",  "Linf64",
                     "Linf65")

valei_largo <- pivot_longer(valei, 
                            cols = c(2:12,), 
                            names_to = "Parameter", 
                            values_to = "SPR")
valei_largo$SPRstra <- rep("EI", nrow(valei_largo))
#Genero lo que quiero comparar Gerlache
valgc <- as.data.frame(cbind(fitgc55@Years,
  fitgc55@SPR,
           fitgc56@SPR,
           fitgc57@SPR,
           fitgc58@SPR,
           fitgc59@SPR,
           fitgc60@SPR,
           fitgc61@SPR,
           fitgc62@SPR,
           fitgc63@SPR,
           fitgc64@SPR,
           fitgc65@SPR))
colnames(valgc) <- c("Years", "Linf55","Linf56", "Linf57", "Linf58", 
                     "Linf59", "Linf60", 
                      "Linf61", "Linf62" , "Linf63",  "Linf64",
                     "Linf65")

valgc_largo <- pivot_longer(valgc, 
                            cols = c(2:12,), 
                            names_to = "Parameter", 
                            values_to = "SPR")
valgc_largo$SPRstra <- rep("GS", nrow(valgc_largo))

#Genero lo que quiero comparar JOin
valjo <- as.data.frame(cbind(fitjo55@Years,
  fitjo55@SPR,
           fitjo56@SPR,
           fitjo57@SPR,
           fitjo58@SPR,
           fitjo59@SPR,
           fitjo60@SPR,
           fitjo61@SPR,
           fitjo62@SPR,
           fitjo63@SPR,
           fitjo64@SPR,
           fitjo65@SPR))
colnames(valjo) <- c("Years", "Linf55","Linf56", "Linf57", "Linf58", 
                     "Linf59", "Linf60", 
                      "Linf61", "Linf62" , "Linf63",  "Linf64",
                     "Linf65")

valjo_largo <- pivot_longer(valjo, 
                            cols = c(2:12,), 
                            names_to = "Parameter", 
                            values_to = "SPR")
valjo_largo$SPRstra <- rep("JOIN", nrow(valjo_largo))

#Genero lo que quiero comparar SSWI
valsswi <- as.data.frame(cbind(fitsswi55@Years,
  fitsswi55@SPR,
           fitsswi56@SPR,
           fitsswi57@SPR,
           fitsswi58@SPR,
           fitsswi59@SPR,
           fitsswi60@SPR,
           fitsswi61@SPR,
           fitsswi62@SPR,
           fitsswi63@SPR,
           fitsswi64@SPR,
           fitsswi65@SPR))
colnames(valsswi) <- c("Years", "Linf55","Linf56", "Linf57", "Linf58", 
                     "Linf59", "Linf60", 
                      "Linf61", "Linf62" , "Linf63",  "Linf64",
                     "Linf65")

valsswi_largo <- pivot_longer(valsswi, 
                            cols = c(2:12,), 
                            names_to = "Parameter", 
                            values_to = "SPR")
valsswi_largo$SPRstra <- rep("SSWI", nrow(valsswi_largo))
#agrupo
valtodo <- rbind(valsswi_largo, 
              valjo_largo,
              valbs_largo,
              valei_largo,
              valgc_largo)
```

Figure \@ref(fig:Figure9) provides a plot by strata with scenarios of L~inf~ of the SPR trends across strata (yellow line = 75% SPR Objective and Red line = 20% Limit SPR).



``` r
sensproto <- ggplot(valtodo %>%
         drop_na() %>% 
           filter(SPR < 0.65),
       aes(x = SPR,
           y = Parameter)) +
  geom_boxplot(trim=TRUE,
                position=position_dodge(0.9),
              adjust = 1/5)+
   geom_vline(xintercept = 0.75,
             colour= '#006d2c',
             alpha=0.5,
             linetype=2)+
  geom_vline(xintercept = 0.20, 
             colour= '#bd0026',
             alpha=0.5)+
  geom_jitter(alpha=0.4,
              color="black",
              height = 0,
              width = 0.1)+
  facet_wrap(~SPRstra,
             ncol=11)+
  labs(x="SPR",
       y="VB Parameters tested")+
  theme_minimal()+
  theme(legend.position="top",
        axis.text.x = element_text(angle = 90, hjust = 2),
        panel.background = element_rect(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlim(0, 1)
sensproto
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/Figure9-1.jpeg" alt="Sensitivity analysis by strata about asymptotic length VB"  />
<p class="caption">(\#fig:Figure9)Sensitivity analysis by strata about asymptotic length VB</p>
</div>



``` r
# Tables
tablinf <- valtodo %>% 
  group_by(Parameter, SPRstra) %>% 
  summarise(Median = mean(SPR),
            Variance=sd(SPR)) %>% 
  pivot_wider(names_from= SPRstra, 
              values_from = c(Median, Variance),
              names_sep = " ") %>% 
  rename( "VB scenario"="Parameter") %>% 
   mutate_if(is.numeric, round, 2)


# Reorganizar la tabla combinando Median y Variance
tablinf_reorganizado <- tablinf %>%
  mutate(
    `BS` = paste0(`Median BS`, " (", `Variance BS`, ")"),
    `EI` = paste0(`Median EI`, " (", `Variance EI`, ")"),
    `GS` = paste0(`Median GS`, " (", `Variance GS`, ")"),
    `JOIN` = paste0(`Median JOIN`, " (", `Variance JOIN`, ")"),
    `SSWI` = paste0(`Median SSWI`, " (", `Variance SSWI`, ")")
  ) %>%
  select(`VB scenario`, `BS`, `EI`, `GS`, `JOIN`, `SSWI`)

#write_csv(tablinf_reorganizado, "SPR_Linf.csv")
```

All the estimated values of SPR by L~inf~scenario by stratum can be identified in  Table \@ref(tab:Table4).


``` r
kbl(tablinf_reorganizado, 
    longtable = FALSE, 
    booktabs = TRUE, 
    caption = "\\label{Table3}Estimated by asymptotyc lenght (VB) scenario (values in parentheses represent the standard deviation)") %>% 
    kable_styling(latex_options = c("scale_down", "hold_position"))
```

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:Table4)\label{Table3}Estimated by asymptotyc lenght (VB) scenario (values in parentheses represent the standard deviation)</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> VB scenario </th>
   <th style="text-align:left;"> BS </th>
   <th style="text-align:left;"> EI </th>
   <th style="text-align:left;"> GS </th>
   <th style="text-align:left;"> JOIN </th>
   <th style="text-align:left;"> SSWI </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Linf55 </td>
   <td style="text-align:left;"> 0.28 (0.07) </td>
   <td style="text-align:left;"> 0.46 (0.23) </td>
   <td style="text-align:left;"> 0.27 (0.11) </td>
   <td style="text-align:left;"> 0.47 (0.07) </td>
   <td style="text-align:left;"> 0.47 (0.13) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Linf56 </td>
   <td style="text-align:left;"> 0.25 (0.07) </td>
   <td style="text-align:left;"> 0.42 (0.23) </td>
   <td style="text-align:left;"> 0.24 (0.09) </td>
   <td style="text-align:left;"> 0.42 (0.07) </td>
   <td style="text-align:left;"> 0.42 (0.11) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Linf57 </td>
   <td style="text-align:left;"> 0.23 (0.06) </td>
   <td style="text-align:left;"> 0.39 (0.23) </td>
   <td style="text-align:left;"> 0.22 (0.08) </td>
   <td style="text-align:left;"> 0.38 (0.07) </td>
   <td style="text-align:left;"> 0.38 (0.1) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Linf58 </td>
   <td style="text-align:left;"> 0.21 (0.06) </td>
   <td style="text-align:left;"> 0.36 (0.23) </td>
   <td style="text-align:left;"> 0.2 (0.08) </td>
   <td style="text-align:left;"> 0.35 (0.06) </td>
   <td style="text-align:left;"> 0.34 (0.09) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Linf59 </td>
   <td style="text-align:left;"> 0.19 (0.05) </td>
   <td style="text-align:left;"> 0.34 (0.23) </td>
   <td style="text-align:left;"> 0.18 (0.07) </td>
   <td style="text-align:left;"> 0.32 (0.06) </td>
   <td style="text-align:left;"> 0.32 (0.08) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Linf60 </td>
   <td style="text-align:left;"> 0.17 (0.05) </td>
   <td style="text-align:left;"> 0.32 (0.23) </td>
   <td style="text-align:left;"> 0.17 (0.06) </td>
   <td style="text-align:left;"> 0.29 (0.06) </td>
   <td style="text-align:left;"> 0.29 (0.07) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Linf61 </td>
   <td style="text-align:left;"> 0.16 (0.04) </td>
   <td style="text-align:left;"> 0.3 (0.21) </td>
   <td style="text-align:left;"> 0.15 (0.06) </td>
   <td style="text-align:left;"> 0.27 (0.06) </td>
   <td style="text-align:left;"> 0.27 (0.07) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Linf62 </td>
   <td style="text-align:left;"> 0.15 (0.04) </td>
   <td style="text-align:left;"> 0.27 (0.18) </td>
   <td style="text-align:left;"> 0.14 (0.05) </td>
   <td style="text-align:left;"> 0.25 (0.05) </td>
   <td style="text-align:left;"> 0.25 (0.06) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Linf63 </td>
   <td style="text-align:left;"> 0.14 (0.04) </td>
   <td style="text-align:left;"> 0.25 (0.17) </td>
   <td style="text-align:left;"> 0.13 (0.05) </td>
   <td style="text-align:left;"> 0.23 (0.05) </td>
   <td style="text-align:left;"> 0.23 (0.06) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Linf64 </td>
   <td style="text-align:left;"> 0.13 (0.04) </td>
   <td style="text-align:left;"> 0.24 (0.15) </td>
   <td style="text-align:left;"> 0.13 (0.05) </td>
   <td style="text-align:left;"> 0.22 (0.05) </td>
   <td style="text-align:left;"> 0.22 (0.05) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Linf65 </td>
   <td style="text-align:left;"> 0.12 (0.04) </td>
   <td style="text-align:left;"> 0.22 (0.14) </td>
   <td style="text-align:left;"> 0.12 (0.04) </td>
   <td style="text-align:left;"> 0.21 (0.05) </td>
   <td style="text-align:left;"> 0.21 (0.05) </td>
  </tr>
</tbody>
</table>

Regarding the three growth levels of the species (high, medium, low) referred to the parameter $k$, it was possible to identify that high and medium growth types result in very low SPR (spawning potential ratio) estimates compared to slow and medium growth. In fact, with high individual growth, the model estimates that SPR levels would be very close to the target management levels (75% SPR) and far from the limit reference level of 20% (Figure \@ref(fig:Figure10),  Table \@ref(tab:Table5)).


``` r
for (i in c(2,0.5,0.3)) {
  nombre_objeto <- paste0("MyPars_", i)  # Crear un nombre único para cada objeto modificado
  objeto_modificado <- MyPars  # Crear una copia del objeto original
  
  objeto_modificado@MK <- i  # Modificar el valor de @Linf en cada iteración
  
  assign(nombre_objeto, objeto_modificado, envir = .GlobalEnv)  # Asignar el objeto modificado al nombre creado
  
  # Resto del código que deseas ejecutar con el objeto modificado
## A blank LB_pars object created
## Default values have been set for some parameters
      MyPars@Species <- "Euphausia superba"
      MyPars@L50 <- 34 
      MyPars@L95 <- 55 # verrificar bibliografia
      MyPars@Linf <- 60
      
      
      #Explotation
      MyPars@SL50 <- 40#numeric() #1
      MyPars@SL95 <- 56#numeric() #27
      MyPars@SPR <- 0.75 #numeric()# ###cambia el numero 0.4 a en blanco
      MyPars@BinWidth <- 1
      #MyPars@FM <- 1
      
      MyPars@Walpha <- 1
      MyPars@Wbeta <- 3.0637 #r2 = 0.9651
      
      MyPars@BinWidth <-1
      MyPars@BinMax <- 70
      MyPars@BinMin <- 0
      MyPars@L_units <- "mm"
  
  # Imprimir el valor de @Linf después de cada cambio
  #print(get(nombre_objeto)@Linf)
}
```


``` r
# ciclo de ajuste para BS
for (i in c(2,0.5,0.3)) {
  nombre_objeto <- paste0("MyPars_", i)  # Nombre del objeto S4 en cada iteración
  nombre_variable <- paste0("fitbs", i)  # Nuevo nombre para el resultado
  
  # Obtener el objeto S4 correspondiente
  objeto <- get(nombre_objeto)
  
  # Ejecutar LBSPRfit() y asignar el resultado a una nueva variable con nombre único
  assign(nombre_variable, LBSPRfit(objeto, Lenbs), envir = .GlobalEnv)
}


# ciclo de ajuste para Gerlache
for (i in c(2,0.5,0.3)) {
  nombre_objeto <- paste0("MyPars_", i)  
  nombre_variable <- paste0("fitei", i)  
  objeto <- get(nombre_objeto)
  assign(nombre_variable, LBSPRfit(objeto, Lenei), envir = .GlobalEnv)
}

# ciclo de ajuste para Gerlache
for (i in c(2,0.5,0.3)) {
  nombre_objeto <- paste0("MyPars_", i)  
  nombre_variable <- paste0("fitgc", i)  
  objeto <- get(nombre_objeto)
  assign(nombre_variable, LBSPRfit(objeto, Lengc), envir = .GlobalEnv)
}

# ciclo de ajuste para JoinVIlle
for (i in c(2,0.5,0.3)) {
  nombre_objeto <- paste0("MyPars_", i) 
  nombre_variable <- paste0("fitjo", i)  
  objeto <- get(nombre_objeto)
  assign(nombre_variable, LBSPRfit(objeto, Lenjo), envir = .GlobalEnv)
}

# ciclo de ajuste para JSSWI
for (i in c(2,0.5,0.3)) {
  nombre_objeto <- paste0("MyPars_", i) 
  nombre_variable <- paste0("fitsswi", i) 
  objeto <- get(nombre_objeto)
  assign(nombre_variable, LBSPRfit(objeto, Lenssiw), envir = .GlobalEnv)
}
```



``` r
#Genero lo que quiero comparar BS
valprobs <- as.data.frame(cbind(fitbs0.3@Years,
                                fitbs0.3@SPR,
                             fitbs0.5@SPR,
                             fitbs2@SPR))
colnames(valprobs) <- c("Years", "High","Med", "Low")

valprobs_largo <- pivot_longer(valprobs, 
                            cols = c(2:4,), 
                            names_to = "Parameter", 
                            values_to = "SPR")
valprobs_largo$SPRstra <- rep("BS", nrow(valprobs_largo))

#Genero lo que quiero comparar EI
valproei <- as.data.frame(cbind(fitei0.3@Years,
                              fitei0.3@SPR,
                             fitei0.5@SPR,
                             fitei2@SPR))
colnames(valproei) <- c("Years", "High","Med", "Low")

valproei_largo <- pivot_longer(valproei, 
                            cols = c(2:4,), 
                            names_to = "Parameter", 
                            values_to = "SPR")
valproei_largo$SPRstra <- rep("EI", nrow(valproei_largo))
#Genero lo que quiero comparar Gerlache
valprogc <- as.data.frame(cbind(fitgc0.3@Years,
                              fitgc0.3@SPR,
                             fitgc0.5@SPR,
                             fitgc2@SPR))
colnames(valprogc) <- c("Years", "High","Med", "Low")

valprogc_largo <- pivot_longer(valprogc, 
                            cols = c(2:4,), 
                            names_to = "Parameter", 
                            values_to = "SPR")
valprogc_largo$SPRstra <- rep("GS", nrow(valprogc_largo))
#Genero lo que quiero comparar JOin
valprojo <- as.data.frame(cbind(fitjo0.3@Years,
                              fitjo0.3@SPR,
                             fitjo0.5@SPR,
                             fitjo2@SPR))
colnames(valprojo) <- c("Years", "High","Med", "Low")

valprojo_largo <- pivot_longer(valprojo, 
                            cols = c(2:4,), 
                            names_to = "Parameter", 
                            values_to = "SPR")
valprojo_largo$SPRstra <- rep("JOIN", nrow(valprojo_largo))

#Genero lo que quiero comparar SSWI
valprosswi <- as.data.frame(cbind(fitsswi0.3@Years,
                              fitsswi0.3@SPR,
                             fitsswi0.5@SPR,
                             fitsswi2@SPR))
colnames(valprosswi) <- c("Years", "High","Med", "Low")

valprosswi_largo <- pivot_longer(valprosswi, 
                            cols = c(2:4,), 
                            names_to = "Parameter", 
                            values_to = "SPR")
valprosswi_largo$SPRstra <- rep("SSWI", nrow(valprosswi_largo))
#agrupo
valprotodo <- rbind(valprosswi_largo, 
              valproei_largo,
              valprobs_largo,
              valprojo_largo,
              valprogc_largo)
```



``` r
tablk <- valprotodo %>% 
  group_by(Parameter, SPRstra) %>% 
  summarise(Median = mean(SPR),
            Variance=var(SPR),
            SD = sqrt(Variance)) %>% 
  pivot_wider(names_from= SPRstra, 
              values_from = c(Median, Variance, SD),
              names_sep = " ") %>% 
  rename( "Growth scenario"="Parameter") %>% 
   mutate_if(is.numeric, round, 2)

# Reorganización de las columnas
tablk_reorganizado <- tablk %>%
  mutate(
    `BS` = paste0(`Median BS`, " (", `SD BS`, ")"),
    `EI` = paste0(`Median EI`, " (", `SD EI`, ")"),
    `GS` = paste0(`Median GS`, " (", `SD GS`, ")"),
    `JOIN` = paste0(`Median JOIN`, " (", `SD JOIN`, ")"),
    `SSWI` = paste0(`Median SSWI`, " (", `SD SSWI`, ")")
  ) %>%
   select(`Growth scenario`, `BS`, `EI`, `GS`, `JOIN`, `SSWI`)  # Selección de columnas necesarias

# Creación de la tabla en LaTeX
kbl(tablk_reorganizado,
    longtable = FALSE,
    booktabs = TRUE,
    caption = "\\label{Table7}Estimates of SPR by Growth Scenario (values in parentheses represent the standard deviation)") %>%
  kable_styling(latex_options = c("scale_down", "hold_position"))
```

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:Table5)\label{Table7}Estimates of SPR by Growth Scenario (values in parentheses represent the standard deviation)</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Growth scenario </th>
   <th style="text-align:left;"> BS </th>
   <th style="text-align:left;"> EI </th>
   <th style="text-align:left;"> GS </th>
   <th style="text-align:left;"> JOIN </th>
   <th style="text-align:left;"> SSWI </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> High </td>
   <td style="text-align:left;"> 0.05 (0.02) </td>
   <td style="text-align:left;"> 0.1 (0.08) </td>
   <td style="text-align:left;"> 0.04 (0.02) </td>
   <td style="text-align:left;"> 0.08 (0.03) </td>
   <td style="text-align:left;"> 0.09 (0.03) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Low </td>
   <td style="text-align:left;"> 0.49 (0.08) </td>
   <td style="text-align:left;"> 0.61 (0.19) </td>
   <td style="text-align:left;"> 0.47 (0.14) </td>
   <td style="text-align:left;"> 0.71 (0.05) </td>
   <td style="text-align:left;"> 0.65 (0.11) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Med </td>
   <td style="text-align:left;"> 0.08 (0.03) </td>
   <td style="text-align:left;"> 0.17 (0.13) </td>
   <td style="text-align:left;"> 0.08 (0.03) </td>
   <td style="text-align:left;"> 0.15 (0.05) </td>
   <td style="text-align:left;"> 0.15 (0.05) </td>
  </tr>
</tbody>
</table>

``` r
#write_csv(tablk_reorganizado, "SPR_K.csv")
```


``` r
sensproto2 <- ggplot(valprotodo %>% 
                       filter(SPR < 0.99) %>% 
                       drop_na(),
                    aes(x = Parameter,
                        y = SPR)) +
  geom_boxplot()+
  geom_hline(yintercept = 0.75,
             colour= '#006d2c',
             alpha=0.5,
             linetype="dashed")+
  geom_hline(yintercept = 0.20, 
             colour= '#bd0026',
             alpha=0.5)+
  geom_jitter(alpha=0.4,
             color="black",
             width=0.1)+
  facet_wrap(~SPRstra,
            ncol=5)+
  labs(y="SPR",
       x="Growth Scenario")+
  theme_minimal()+
  theme(legend.position="none",
        panel.background = element_rect(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_viridis_d(option="F",
                       name="Strata")
sensproto2
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/Figure10-1.jpeg" alt="Sensitivity analysis by strata about krill growth type" width="70%" />
<p class="caption">(\#fig:Figure10)Sensitivity analysis by strata about krill growth type</p>
</div>

## Analysis to four Subareas Outputs (48.1, 48.2 and 48.3)

Reading data


``` r
datdir <- setwd("~/DOCAS/LBSPR_Krill")
Len481 <- new("LB_lengths", 
             LB_pars=MyPars, 
             file=paste0(datdir, "/lenght481_nochina.csv"), 
             dataType="freq",
             sep=",",
             header=T)

Len482 <- new("LB_lengths", 
             LB_pars=MyPars, 
             file=paste0(datdir, "/lenght482_nochina.csv"),
             dataType="freq",
             sep=",",
             header=T)

Len483 <- new("LB_lengths", 
             LB_pars=MyPars, 
             file=paste0(datdir, "/lenght483_nochina.csv"), 
             dataType="freq",
             sep=",",
             header=T)
```

The LBSPR model is fitted by strata using the `LBSPRfit` function.


``` r
fit481 <- LBSPRfit(MyPars, Len481)
fit482 <- LBSPRfit(MyPars, Len482)
fit483 <- LBSPRfit(MyPars, Len483)
```

Now, to paper propose, we estimate `SPR` to all subareas (48.1, 48.2 and 48.3)


``` r
# Get values estimated

spr481 <- as.data.frame(cbind(fit481 @Years, 
                             fit481 @SPR,
                             fit481 @Vars[,4])) 
colnames(spr481) <- c("Year","SPR", "Var")
spr481$SPRv <- rep("481", nrow(spr481))

spr482 <- as.data.frame(cbind(fit482 @Years, 
                             fit482 @SPR,
                             fit482 @Vars[,4])) 
colnames(spr482) <- c("Year","SPR", "Var")
spr482$SPRv <- rep("482", nrow(spr482))

spr483 <- as.data.frame(cbind(fit483 @Years, 
                             fit483 @SPR,
                             fit483 @Vars[,4])) 
colnames(spr483) <- c("Year","SPR", "Var")
spr483$SPRv <- rep("483", nrow(spr483))


spr48 <- rbind(spr481,
               spr482,
               spr483)

spr48 <- spr48 %>%
  mutate(SD = sqrt(Var))

spr48wide <- round(pivot_wider(spr48,
                          names_from = SPRv,
                          values_from = c(SPR, Var, SD),
                          values_fill = NA,
                          names_sep = " ",),3)
```

Figure \@ref(fig:Figure11) show the SPR trends across years to 48.1, 48.2 and 48.3, indicating the references (yellow line = 75% SPR Objective and Red line = 20% Limit SPR).


``` r
spr48plot <- ggplot(spr48, aes(Year, SPR)) +
  geom_point(alpha = 0.8) +  # Puntos de datos
  geom_errorbar(aes(ymin = SPR - SD, ymax = SPR + SD),  # Barras de error
                width = 0.2,  # Ancho de la barra horizontal
                alpha = 0.5) +
  geom_smooth(method = "lm",
              alpha=0.3,
              se=TRUE,
              col="black",
              level = 0.75)+
  geom_hline(yintercept = 0.75,
             colour = '#006d2c',
             alpha = 0.5,
             linetype = 2) +
  geom_hline(yintercept = 0.20,
             colour = '#bd0026',
             alpha = 0.5) +
  facet_wrap(. ~ SPRv, ncol = 5) +
  xlim(2001, 2021) +
  ylim(0, 0.9) +
  theme_few() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
spr48plot
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/Figure11-1.jpeg" alt="Krill Intrinsic Productivity (SPR) by subarea and by year"  />
<p class="caption">(\#fig:Figure11)Krill Intrinsic Productivity (SPR) by subarea and by year</p>
</div>
Now, just 48.1 (Figure \@ref(fig:Figure12)).


``` r
spr481plot <- ggplot(spr48 %>% 
                       filter(SPRv == "481"), aes(Year, SPR)) +
  geom_point(alpha = 0.8) +  # Puntos de datos
  geom_errorbar(aes(ymin = SPR - SD, ymax = SPR + SD),  # Barras de error
                width = 0.2,  # Ancho de la barra horizontal
                alpha = 0.5) +
  geom_smooth(method = "lm",
              alpha=0.3,
              se=TRUE,
              col="black",
              level = 0.75)+
  geom_hline(yintercept = 0.75,
             colour = '#006d2c',
             alpha = 0.5,
             linetype = 2) +
  geom_hline(yintercept = 0.20,
             colour = '#bd0026',
             alpha = 0.5) +
  xlim(2001, 2021) +
  ylim(0, 0.9) +
  theme_few() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
spr481plot
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/Figure12-1.jpeg" alt="Krill Intrinsic Productivity (SPR) in 48.1 by year"  />
<p class="caption">(\#fig:Figure12)Krill Intrinsic Productivity (SPR) in 48.1 by year</p>
</div>
Main outputs in Table \@ref(tab:Table6).

``` r
spr48wide_reorganizado <- spr48wide %>%
  mutate(
    `481` = ifelse(!is.na(`SPR 481`), paste0(`SPR 481`, " (", `SD 481`, ")"), NA),
    `482` = ifelse(!is.na(`SPR 482`), paste0(`SPR 482`, " (", `SD 482`, ")"), NA),
    `483` = ifelse(!is.na(`SPR 483`), paste0(`SPR 483`, " (", `SD 483`, ")"), NA)
  ) %>%
  select(Year, `481`, `482`, `483`) %>% 
  arrange(Year)

#write_csv(spr48wide_reorganizado, "SPR_all_48.csv")

kbl(spr48wide_reorganizado, 
    longtable = FALSE, 
    booktabs = TRUE, 
    caption = "\\label{Table6}Estimates of SPR by SubArea (parentheses represent the standard deviation)") %>% 
    kable_styling(latex_options = c("scale_down", "hold_position"))
```

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
<caption>(\#tab:Table6)\label{Table6}Estimates of SPR by SubArea (parentheses represent the standard deviation)</caption>
 <thead>
  <tr>
   <th style="text-align:right;"> Year </th>
   <th style="text-align:left;"> 481 </th>
   <th style="text-align:left;"> 482 </th>
   <th style="text-align:left;"> 483 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 2001 </td>
   <td style="text-align:left;"> 0.169 (0.017) </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2002 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.199 (0.008) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2004 </td>
   <td style="text-align:left;"> 0.158 (0.009) </td>
   <td style="text-align:left;"> 0.274 (0.004) </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2005 </td>
   <td style="text-align:left;"> 0.253 (0.004) </td>
   <td style="text-align:left;"> 0.313 (0.019) </td>
   <td style="text-align:left;"> 0.109 (0.001) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2006 </td>
   <td style="text-align:left;"> 0.262 (0.013) </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0.185 (0.007) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2007 </td>
   <td style="text-align:left;"> 0.199 (0.01) </td>
   <td style="text-align:left;"> 0.112 (0.003) </td>
   <td style="text-align:left;"> 0.082 (0.003) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2008 </td>
   <td style="text-align:left;"> 0.069 (0.008) </td>
   <td style="text-align:left;"> 0.165 (0.008) </td>
   <td style="text-align:left;"> 0.189 (0.004) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2009 </td>
   <td style="text-align:left;"> 0.264 (0.019) </td>
   <td style="text-align:left;"> 0.328 (0.009) </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2010 </td>
   <td style="text-align:left;"> 0.225 (0.003) </td>
   <td style="text-align:left;"> 0.426 (0.004) </td>
   <td style="text-align:left;"> 0.154 (0.007) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2011 </td>
   <td style="text-align:left;"> 0.409 (0.009) </td>
   <td style="text-align:left;"> 0.235 (0.004) </td>
   <td style="text-align:left;"> 0.143 (0.002) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2012 </td>
   <td style="text-align:left;"> 0.412 (0.006) </td>
   <td style="text-align:left;"> 0.089 (0.004) </td>
   <td style="text-align:left;"> 0.187 (0.002) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2013 </td>
   <td style="text-align:left;"> 0.135 (0.001) </td>
   <td style="text-align:left;"> 0.259 (0.007) </td>
   <td style="text-align:left;"> 0.292 (0.003) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2014 </td>
   <td style="text-align:left;"> 0.215 (0.002) </td>
   <td style="text-align:left;"> 0.25 (0.009) </td>
   <td style="text-align:left;"> 0.155 (0.002) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2015 </td>
   <td style="text-align:left;"> 0.215 (0.002) </td>
   <td style="text-align:left;"> 0.151 (0.005) </td>
   <td style="text-align:left;"> 0.229 (0.002) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2016 </td>
   <td style="text-align:left;"> 0.285 (0.008) </td>
   <td style="text-align:left;"> 0.31 (0.011) </td>
   <td style="text-align:left;"> 0.188 (0.003) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2017 </td>
   <td style="text-align:left;"> 0.217 (0.004) </td>
   <td style="text-align:left;"> 0.158 (0.003) </td>
   <td style="text-align:left;"> 0.178 (0.007) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2018 </td>
   <td style="text-align:left;"> 0.224 (0.002) </td>
   <td style="text-align:left;"> 0.186 (0.006) </td>
   <td style="text-align:left;"> 0.31 (0.016) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2019 </td>
   <td style="text-align:left;"> 0.219 (0.004) </td>
   <td style="text-align:left;"> 0.335 (0.007) </td>
   <td style="text-align:left;"> 0.143 (0.004) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2020 </td>
   <td style="text-align:left;"> 0.114 (0.003) </td>
   <td style="text-align:left;"> 0.217 (0.003) </td>
   <td style="text-align:left;"> 0.224 (0.002) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2021 </td>
   <td style="text-align:left;"> 0.225 (0.002) </td>
   <td style="text-align:left;"> 0.379 (0.004) </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2022 </td>
   <td style="text-align:left;"> 0.876 (0.006) </td>
   <td style="text-align:left;"> 0.271 (0.004) </td>
   <td style="text-align:left;"> 0.244 (0.005) </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2023 </td>
   <td style="text-align:left;"> 0.359 (0.003) </td>
   <td style="text-align:left;"> 0.358 (0.006) </td>
   <td style="text-align:left;"> 0.311 (0.004) </td>
  </tr>
</tbody>
</table>

sensitivity analysis to L~inf~ to 48.1


``` r
# ciclo de ajuste para 481
for (i in 55:65) {
  nombre_objeto <- paste0("MyPars_", i) 
  nombre_variable <- paste0("fit481", i) 
  objeto <- get(nombre_objeto)
  assign(nombre_variable, LBSPRfit(objeto, Len481), envir = .GlobalEnv)
}
```



``` r
#Genero lo que quiero comparar BS
val481 <- as.data.frame(cbind(fit48155@Years,
  fit48155@SPR,
           fit48156@SPR,
           fit48157@SPR,
           fit48158@SPR,
           fit48159@SPR,
           fit48160@SPR,
           fit48161@SPR,
           fit48162@SPR,
           fit48163@SPR,
           fit48164@SPR,
           fit48165@SPR))
colnames(val481) <- c("Years", "Linf55","Linf56", "Linf57", "Linf58", 
                     "Linf59", "Linf60", 
                      "Linf61", "Linf62" , "Linf63",  "Linf64",
                     "Linf65")

val481_largo <- pivot_longer(val481, 
                            cols = c(2:12,), 
                            names_to = "Parameter", 
                            values_to = "SPR")
val481_largo$SPRstra <- rep("481", nrow(val481_largo))
```


Now, boxplot just 48.1 (Figure \@ref(fig:Figure13)).



``` r
sensprototest481 <- ggplot(val481_largo %>%
         drop_na() %>% 
           filter(SPR < 0.65),
       aes(x = SPR,
           y = Parameter)) +
  geom_boxplot(trim=TRUE,
                position=position_dodge(0.9),
              adjust = 1/5)+
   geom_vline(xintercept = 0.75,
             colour= '#006d2c',
             alpha=0.5,
             linetype=2)+
  geom_vline(xintercept = 0.20, 
             colour= '#bd0026',
             alpha=0.5)+
  geom_jitter(alpha=0.4,
              color="black",
              height = 0,
              width = 0.1)+
  labs(x="SPR",
       y="VB Parameters tested")+
  theme_minimal()+
  theme(legend.position="top",
        axis.text.x = element_text(angle = 90, hjust = 2),
        panel.background = element_rect(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlim(0, 1)
sensprototest481
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/Figure13-1.jpeg" alt="Sensitivity analysis to 48.1 about asymptotic length VB"  />
<p class="caption">(\#fig:Figure13)Sensitivity analysis to 48.1 about asymptotic length VB</p>
</div>


Now, sensitivity analysis to *k* to 48.1


``` r
# ciclo de ajuste para 481
for (i in c(2,0.5,0.3)) {
  nombre_objeto <- paste0("MyPars_", i)  # Nombre del objeto S4 en cada iteración
  nombre_variable <- paste0("fit481", i)  # Nuevo nombre para el resultado
  
  # Obtener el objeto S4 correspondiente
  objeto <- get(nombre_objeto)
  
  # Ejecutar LBSPRfit() y asignar el resultado a una nueva variable con nombre único
  assign(nombre_variable, LBSPRfit(objeto, Len481), envir = .GlobalEnv)
}
```




``` r
#Genero lo que quiero comparar BS
valprofit481 <- as.data.frame(cbind(fit4810.3@Years,
                                fit4810.3@SPR,
                             fit4810.5@SPR,
                             fit4812@SPR))
colnames(valprofit481) <- c("Years", "High","Med", "Low")

valprofit481_largo <- pivot_longer(valprofit481, 
                            cols = c(2:4,), 
                            names_to = "Parameter", 
                            values_to = "SPR")
valprofit481_largo$SPRstra <- rep("481", nrow(valprofit481_largo))
```

And, plot to  k scenarios (Figure \@ref(fig:Figure14)).



``` r
senspr481 <- ggplot(valprofit481_largo %>% 
                       filter(SPR < 0.99) %>% 
                       drop_na()%>%
                       mutate(Parameter = factor(Parameter, 
                                                 levels = c("Low", "Med", "High"))), # Orden deseado
                    aes(x = Parameter, 
                        y = SPR)) +
  geom_boxplot()+
  geom_hline(yintercept = 0.75,
             colour= '#006d2c',
             alpha=0.5,
             linetype="dashed")+
  geom_hline(yintercept = 0.20, 
             colour= '#bd0026',
             alpha=0.5)+
  geom_jitter(alpha=0.4,
             color="black",
             width=0.1)+
  labs(y="SPR",
       x="Growth Scenario")+
  theme_minimal()+
  theme(legend.position="none",
        panel.background = element_rect(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_viridis_d(option="F",
                       name="Strata")
senspr481
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/Figure14-1.jpeg" alt="Sensitivity analysis to 48.1 about k"  />
<p class="caption">(\#fig:Figure14)Sensitivity analysis to 48.1 about k</p>
</div>

Finally, join plot to 48.1 analysis; `spr481plot`, `sensprototest481`, `senspro481` (Figure \@ref(fig:Figure15)).



``` r
ggarrange(spr481plot,
          sensprototest481,
          senspr481,
          ncol=3)
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/Figure15-1.jpeg" alt="Plot to SPR in different context to 48.1 subarea"  />
<p class="caption">(\#fig:Figure15)Plot to SPR in different context to 48.1 subarea</p>
</div>

# CODE REPOSITORY

The data, codes and another documents of this exercise can be found in the following link [LBSPR-Krill](https://github.com/MauroMardones/LBSPR_Krill)


# REFERENCES
