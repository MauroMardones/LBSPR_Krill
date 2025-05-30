---
title: "Supporting Information 1"
subtitle: "Disparate estimates of intrinsic productivity for Antarctic krill (Euphausia superba) across small spatial scales, under a rapidly changing ocean."
author: "Mardones, M; Jarvis Mason, E.T.;  Santa Cruz, F.; Watters, G.; Cárdenas, C.A"
date:  "29 May, 2025"
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

\begin{figure}[H]

{\centering \includegraphics[width=0.6\linewidth]{index_files/figure-html/Strata_Nochina} 

}

\caption{Sizes compositions from SISO program monitoring krill fishery by strata (BS=Brainsfield Strait, EI= Elephant Island, GS= Gerlache Strait, JOIN= Joinville Island, SSWI= South West). Red line represent recruit size}(\#fig:Figure1)
\end{figure}

The information gaps (years without sizes composition data) are not calculated because there is no autocorrelation between years, but singular estimators over time.


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

\begin{table}[!h]
\centering
\caption{(\#tab:Table1)Krill biological and fishery parameters}
\centering
\begin{tabular}[t]{ll}
\toprule
Value & Descrption\\
\midrule
\cellcolor{gray!10}{60} & \cellcolor{gray!10}{VB asymptotic length}\\
34 & Maturity 50\%\\
\cellcolor{gray!10}{55} & \cellcolor{gray!10}{Maturity 95\%}\\
0.889 & M/K Ratio\\
\cellcolor{gray!10}{40} & \cellcolor{gray!10}{Selectivity 50\%}\\
\addlinespace
56 & Seletivity 95\%\\
\cellcolor{gray!10}{0.75} & \cellcolor{gray!10}{SPR}\\
1 & a (Length-Weight Relation)\\
\cellcolor{gray!10}{3.0637} & \cellcolor{gray!10}{b (Length-Weight Relation)}\\
70 & Bin Min\\
\addlinespace
\cellcolor{gray!10}{0} & \cellcolor{gray!10}{Bin Max}\\
mm & Units\\
\bottomrule
\end{tabular}
\end{table}


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

\begin{figure}[H]

{\centering \includegraphics[width=0.7\linewidth]{index_files/figure-latex/Figure4-1} 

}

\caption{Theoretical growth curves for krill based on three SPR sensitivity scenarios that were tested using LBSPR.}(\#fig:Figure4)
\end{figure}

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

\begin{figure}[H]

{\centering \includegraphics{index_files/figure-latex/Figure5-1} 

}

\caption{Fit of the model to the data of lengths in all strata}(\#fig:Figure5)
\end{figure}


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

\begin{figure}[H]

{\centering \includegraphics[width=0.8\linewidth]{index_files/figure-latex/Figure6-1} 

}

\caption{Difference between the observed accumulated size structure for each stratum related SPR objective}(\#fig:Figure6)
\end{figure}

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

\begin{figure}[H]

{\centering \includegraphics[width=0.8\linewidth]{index_files/figure-latex/Figure7-1} 

}

\caption{Maturity curves by strata}(\#fig:Figure7)
\end{figure}

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

\begin{figure}[H]

{\centering \includegraphics{index_files/figure-latex/Figure8-1} 

}

\caption{Krill Intrinsic Productivity (SPR) by strata and by year}(\#fig:Figure8)
\end{figure}

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

\begin{table}[!h]
\centering
\caption{(\#tab:Table2)Estimates of SPR by Strata (values in parentheses represent the standard deviation)}
\centering
\resizebox{\ifdim\width>\linewidth\linewidth\else\width\fi}{!}{
\begin{tabular}[t]{rlllll}
\toprule
Year & BS & EI & GS & JOIN & SSWI\\
\midrule
2001 & 0.132 (0.024) & 0.223 (0.062) & NA & NA & 0.196 (0.011)\\
2004 & 0.15 (0.008) & 0.219 (0.028) & NA & NA & 0.21 (0.045)\\
2005 & 0.25 (0.014) & NA & NA & NA & 0.253 (0.005)\\
2006 & 0.214 (0.009) & NA & 0.152 (0.09) & NA & 0.349 (0.025)\\
2007 & 0.108 (0.078) & 0.095 (0.008) & 0.069 (0.014) & NA & 0.277 (0.057)\\
\addlinespace
2008 & NA & 0.069 (0.008) & NA & NA & NA\\
2009 & NA & 0.322 (0.058) & NA & NA & 0.232 (0.016)\\
2010 & 0.197 (0.004) & 0.329 (0.02) & 0.242 (0.014) & NA & 0.31 (0.008)\\
2011 & 0.185 (0.012) & NA & NA & NA & 0.426 (0.009)\\
2012 & 0.156 (0.004) & NA & 0.148 (0.003) & NA & 0.386 (0.006)\\
\addlinespace
2013 & 0.13 (0.001) & 0.191 (0.013) & 0.123 (0.002) & NA & 0.176 (0.005)\\
2014 & 0.175 (0.002) & 0.201 (0.007) & 0.172 (0.006) & NA & 0.311 (0.004)\\
2015 & 0.19 (0.003) & 0.389 (0.022) & 0.177 (0.003) & NA & 0.386 (0.026)\\
2016 & 0.288 (0.01) & 0.36 (0.057) & 0.293 (0.007) & 0.331 (0.19) & 0.266 (0.033)\\
2017 & 0.165 (0.004) & 1 (0) & 0.214 (0.016) & 0.249 (0.022) & 0.27 (0.01)\\
\addlinespace
2018 & 0.165 (0.003) & 0.334 (0.019) & 0.162 (0.004) & NA & 0.294 (0.018)\\
2019 & 0.161 (0.004) & 0.425 (0.02) & 0.141 (0.005) & NA & 0.349 (0.006)\\
2020 & 0.11 (0.004) & NA & 0.087 (0.004) & NA & 0.24 (0.012)\\
\bottomrule
\end{tabular}}
\end{table}


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
# Crear una nueva columna para identificar los puntos de Linf60
valtodo <- valtodo %>%
  mutate(Highlight = ifelse(Parameter == "Linf60", "Highlighted", "Normal"))

# Crear el gráfico con jitter coloridos y boxplots vacíos con borde negro
sensproto <- ggplot(valtodo %>%
                      drop_na() %>%
                      filter(SPR < 0.65),
                    aes(x = SPR,
                        y = Parameter)) +
  geom_boxplot(aes(fill = ifelse(Parameter == "Linf60", "red", NA)),  # Destacar el boxplot de Linf60 en rojo
               trim = TRUE,
               color = "black",  # Borde negro
               position = position_dodge(0.9),
               adjust = 1/5) +
  geom_vline(xintercept = 0.75,
             colour = '#006d2c',
             alpha = 0.5,
             linetype = 2) +
  geom_vline(xintercept = 0.20,
             colour = '#bd0026',
             alpha = 0.5) +
  geom_jitter(aes(fill = ifelse(Parameter == "Linf60", "red", NA)),  # Destacar los puntos de Linf60 en rojo
              shape = 21,  # Forma de puntos vacíos con borde
              color = "black",  # Borde negro para todos los puntos
              alpha = 0.5,  # Transparencia
              height = 0,
              width = 0.1) +
  facet_wrap(~SPRstra,
             ncol = 11) +
  labs(x = "SPR",
       y = "VB Parameters tested") +
  theme_minimal() +
  theme(legend.position = "none",  # Ocultamos la leyenda
        axis.text.x = element_text(angle = 90, hjust = 2),
        panel.background = element_rect(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlim(0, 1) +
  scale_fill_identity()  # Usar el color definido en aes(fill) directamente

sensproto
```

\begin{figure}[H]

{\centering \includegraphics{index_files/figure-latex/Figure9-1} 

}

\caption{Sensitivity analysis by strata about asymptotic length VB}(\#fig:Figure9)
\end{figure}



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

\begin{table}[!h]
\centering
\caption{(\#tab:Table4)\label{Table3}Estimated by asymptotyc lenght (VB) scenario (values in parentheses represent the standard deviation)}
\centering
\resizebox{\ifdim\width>\linewidth\linewidth\else\width\fi}{!}{
\begin{tabular}[t]{llllll}
\toprule
VB scenario & BS & EI & GS & JOIN & SSWI\\
\midrule
Linf55 & 0.28 (0.07) & 0.46 (0.23) & 0.27 (0.11) & 0.47 (0.07) & 0.47 (0.13)\\
Linf56 & 0.25 (0.07) & 0.42 (0.23) & 0.24 (0.09) & 0.42 (0.07) & 0.42 (0.11)\\
Linf57 & 0.23 (0.06) & 0.39 (0.23) & 0.22 (0.08) & 0.38 (0.07) & 0.38 (0.1)\\
Linf58 & 0.21 (0.06) & 0.36 (0.23) & 0.2 (0.08) & 0.35 (0.06) & 0.34 (0.09)\\
Linf59 & 0.19 (0.05) & 0.34 (0.23) & 0.18 (0.07) & 0.32 (0.06) & 0.32 (0.08)\\
\addlinespace
Linf60 & 0.17 (0.05) & 0.32 (0.23) & 0.17 (0.06) & 0.29 (0.06) & 0.29 (0.07)\\
Linf61 & 0.16 (0.04) & 0.3 (0.21) & 0.15 (0.06) & 0.27 (0.06) & 0.27 (0.07)\\
Linf62 & 0.15 (0.04) & 0.27 (0.18) & 0.14 (0.05) & 0.25 (0.05) & 0.25 (0.06)\\
Linf63 & 0.14 (0.04) & 0.25 (0.17) & 0.13 (0.05) & 0.23 (0.05) & 0.23 (0.06)\\
Linf64 & 0.13 (0.04) & 0.24 (0.15) & 0.13 (0.05) & 0.22 (0.05) & 0.22 (0.05)\\
\addlinespace
Linf65 & 0.12 (0.04) & 0.22 (0.14) & 0.12 (0.04) & 0.21 (0.05) & 0.21 (0.05)\\
\bottomrule
\end{tabular}}
\end{table}

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

\begin{table}[!h]
\centering
\caption{(\#tab:Table5)\label{Table7}Estimates of SPR by Growth Scenario (values in parentheses represent the standard deviation)}
\centering
\resizebox{\ifdim\width>\linewidth\linewidth\else\width\fi}{!}{
\begin{tabular}[t]{llllll}
\toprule
Growth scenario & BS & EI & GS & JOIN & SSWI\\
\midrule
High & 0.05 (0.02) & 0.1 (0.08) & 0.04 (0.02) & 0.08 (0.03) & 0.09 (0.03)\\
Low & 0.49 (0.08) & 0.61 (0.19) & 0.47 (0.14) & 0.71 (0.05) & 0.65 (0.11)\\
Med & 0.08 (0.03) & 0.17 (0.13) & 0.08 (0.03) & 0.15 (0.05) & 0.15 (0.05)\\
\bottomrule
\end{tabular}}
\end{table}

``` r
#write_csv(tablk_reorganizado, "SPR_K.csv")
```


``` r
# Crear el gráfico con los elementos resaltados para Parameter == "Med"
sensproto2 <- ggplot(valprotodo %>% 
                       filter(SPR < 0.99) %>% 
                       drop_na() %>% 
                       mutate(Parameter = factor(Parameter, 
                                                 levels = c("Low", "Med", "High"))), # Orden deseado
                     aes(x = Parameter, y = SPR)) +
  # Boxplot con color condicional
  geom_boxplot(aes(fill = Parameter), alpha = 0.6) +
  # Líneas horizontales de referencia
  geom_hline(yintercept = 0.75, colour= '#006d2c', alpha=0.5, linetype="dashed") +
  geom_hline(yintercept = 0.20, colour= '#bd0026', alpha=0.5) +
  # Jitter con color condicional
  geom_jitter(aes(color = Parameter), alpha = 0.5, width = 0.1) +
  # Facet por SPRstra
  facet_wrap(~SPRstra, ncol = 5) +
  # Etiquetas
  labs(y = "SPR", x = "Growth Scenario") +
  # Tema minimalista
  theme_minimal() +
  theme(legend.position = "none",
        panel.background = element_rect(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  # Especificar colores manualmente para Parameter
  scale_fill_manual(values = c("Low" = "white", "Med" = "red", "High" = "white")) +
  scale_color_manual(values = c("Low" = "black", "Med" = "red", "High" = "black"))

# Mostrar el gráfico
sensproto2
```

\begin{figure}[H]

{\centering \includegraphics{index_files/figure-latex/Figure10-1} 

}

\caption{Sensitivity analysis by strata about krill growth type}(\#fig:Figure10)
\end{figure}

## Analysis to three Subareas (48.1, 48.2 and 48.3)

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

\begin{figure}[H]

{\centering \includegraphics{index_files/figure-latex/Figure11-1} 

}

\caption{Krill Intrinsic Productivity (SPR) by subarea and by year}(\#fig:Figure11)
\end{figure}

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

\begin{table}[!h]
\centering
\caption{(\#tab:Table6)\label{Table6}Estimates of SPR by SubArea (parentheses represent the standard deviation)}
\centering
\resizebox{\ifdim\width>\linewidth\linewidth\else\width\fi}{!}{
\begin{tabular}[t]{rlll}
\toprule
Year & 481 & 482 & 483\\
\midrule
2001 & 0.169 (0.017) & NA & NA\\
2002 & NA & NA & 0.199 (0.008)\\
2004 & 0.158 (0.009) & 0.274 (0.004) & NA\\
2005 & 0.253 (0.004) & 0.313 (0.019) & 0.109 (0.001)\\
2006 & 0.262 (0.013) & NA & 0.185 (0.007)\\
\addlinespace
2007 & 0.199 (0.01) & 0.112 (0.003) & 0.082 (0.003)\\
2008 & 0.069 (0.008) & 0.165 (0.008) & 0.189 (0.004)\\
2009 & 0.264 (0.019) & 0.328 (0.009) & NA\\
2010 & 0.225 (0.003) & 0.426 (0.004) & 0.154 (0.007)\\
2011 & 0.409 (0.009) & 0.235 (0.004) & 0.143 (0.002)\\
\addlinespace
2012 & 0.412 (0.006) & 0.089 (0.004) & 0.187 (0.002)\\
2013 & 0.135 (0.001) & 0.259 (0.007) & 0.292 (0.003)\\
2014 & 0.215 (0.002) & 0.25 (0.009) & 0.155 (0.002)\\
2015 & 0.215 (0.002) & 0.151 (0.005) & 0.229 (0.002)\\
2016 & 0.285 (0.008) & 0.31 (0.011) & 0.188 (0.003)\\
\addlinespace
2017 & 0.217 (0.004) & 0.158 (0.003) & 0.178 (0.007)\\
2018 & 0.224 (0.002) & 0.186 (0.006) & 0.31 (0.016)\\
2019 & 0.219 (0.004) & 0.335 (0.007) & 0.143 (0.004)\\
2020 & 0.114 (0.003) & 0.217 (0.003) & 0.224 (0.002)\\
2021 & 0.225 (0.002) & 0.379 (0.004) & NA\\
\addlinespace
2022 & 0.876 (0.006) & 0.271 (0.004) & 0.244 (0.005)\\
2023 & 0.359 (0.003) & 0.358 (0.006) & 0.311 (0.004)\\
\bottomrule
\end{tabular}}
\end{table}


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
```


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




``` r
sensprototest481 <- ggplot(val481_largo %>%
         drop_na() %>% 
           filter(SPR < 0.65),
       aes(x = SPR,
           y = Parameter)) +
  geom_boxplot(aes(fill = ifelse(Parameter == "Linf60", "red", "white")),  # Destacar el boxplot de Linf60 en rojo
               trim=TRUE,
               position=position_dodge(0.9),
               adjust = 1/5) +
  geom_vline(xintercept = 0.75,
             colour= '#006d2c',
             alpha=0.5,
             linetype=2) +
  geom_vline(xintercept = 0.20, 
             colour= '#bd0026',
             alpha=0.5) +
  geom_jitter(aes(color = ifelse(Parameter == "Linf60", "red", "black")),  # Destacar los puntos de Linf60 en rojo
              alpha=0.4,
              height = 0,
              width = 0.1) +
  labs(x="SPR",
       y="VB Parameters tested") +
  theme_minimal() +
  theme(legend.position="top",
        axis.text.x = element_text(angle = 90, hjust = 2),
        panel.background = element_rect(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlim(0, 1) +
  scale_fill_identity() +  # Usar el color definido en aes(fill) directamente
  scale_color_identity()   # Usar el color definido en aes(color) directamente
```


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




``` r
senspr481 <- ggplot(valprofit481_largo %>% 
                       filter(SPR < 0.99) %>% 
                       drop_na() %>%
                       mutate(Parameter = factor(Parameter, 
                                                 levels = c("Low", "Med", "High"))), # Orden deseado
                    aes(x = Parameter, 
                        y = SPR)) +
  geom_boxplot(aes(fill = ifelse(Parameter == "Med", "red", "white")),  # Destacar el boxplot de "Med" en rojo
               color = "black") +  # Color del borde del boxplot
  geom_hline(yintercept = 0.75,
             colour= '#006d2c',
             alpha=0.5,
             linetype="dashed") +
  geom_hline(yintercept = 0.20, 
             colour= '#bd0026',
             alpha=0.5) +
  geom_jitter(aes(color = ifelse(Parameter == "Med", "red", "black")),  # Destacar los puntos de "Med" en rojo
              alpha=0.4,
              width=0.1) +
  labs(y="SPR",
       x="Growth Scenario") +
  theme_minimal() +
  theme(legend.position="none",
        panel.background = element_rect(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_identity() +  # Usar el color definido en aes(fill) directamente
  scale_color_identity()   # Usar el color definido en aes(color) directamente
```

Finally, join plot to 48.1 analysis; `spr481plot`, `sensprototest481`, `senspro481` (Figure \@ref(fig:Figure15)).



``` r
ggarrange(spr481plot,
          sensprototest481,
          senspr481,
          ncol=3,
          labels = c("A", "B", "C"))
```

\begin{figure}[H]

{\centering \includegraphics{index_files/figure-latex/Figure15-1} 

}

\caption{Plot to SPR in different context to 48.1 subarea}(\#fig:Figure15)
\end{figure}

# CODE REPOSITORY

The data, codes and another documents of this exercise can be found in the following link [LBSPR-Krill](https://github.com/MauroMardones/LBSPR_Krill)


# REFERENCES
