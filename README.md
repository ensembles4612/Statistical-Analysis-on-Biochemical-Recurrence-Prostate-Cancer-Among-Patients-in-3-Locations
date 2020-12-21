# Project Overview: Statistical Analysis on Biochemical Recurrence Prostate Cancer Among Patients in 3 Locations
 
## Project Highlights

* Computed Patient Baseline Characteristics table and Kaplan-Meier Curve graphs regarding 3 locations for prostate cancer patients
* Proved probabilities of freedom from BCR among three locations were significantly different after stratification by conducting Stratified Logrank Test
* built a AIC Stepwise Cox’s proportional hazards regression model with imputation and interaction terms, obtaining which demographic factors and genes have effects on survival
* Found which factors and genes have greater effects on survival using Random Survival Forest method

## Code

* **Data Source**: 3 datasets are from built-in libraries: library(prostateCancerCamcap), library(prostateCancerTaylor), library(prostateCancerStockholm)
* **Language Used**: R
* **Libraries Used**: ggplot2, tidyr, devtools, dplyr, party, survival, tcltk, tcltk2, forcats, compareGroups, rms, npsurv, MASS, survminer, randomForestSRC

## Project Objective

* To investigate whether the target genes as well as demographic factors have effects on biochemical recurrence(BCR) free time among prostate cancer patients in three public datasets, Taylor, Cambridge and Stockholm (see Section 4)
* To test if BCR free times are significantly different among three datasets (see Section 2 and Section 3)
* To investigate which factors and genes have greater effects on survival using Random Survival Forests (Section 7)

## Project Background

Biochemical recurrence(BCR) continues to occur in a large proportion of prostate cancer patients. Yet, established clinical predictors(PSA, Gleason score) often provide poor prognosis. It’s noticeable that more research is now focusing on the effect of genes such as integrated genomics and gene signatures on prostate cancer BCR.

## Data Cleaning 

Some main cleaning I did as follows: 
* Cleaned 3 datasets (3 locations) respectively to prepare them in the same format for horizontal join such as 
  * extracting Gene expressions, survival and demographic variables
  * making variables the same names with the same data types for the 3 datasets
  * dropping the rows having NAs in Time and Event indicator
* Selected variables in the same order for all datasets and fully joined them horizontally 
* Created new variables such as 'dataset' with 3 levels (Cam, Taylor and stoch) to distinguish data from 3 different locations in the fully joined dataset
* Cleaned categorical variables with multiple levels such as combining values with the same meaning to one level

After cleaning, we learned the following from the final joined dataset:

* Columns included: 
  * Age: patient’s age
  * Time:  patient’s biochemical recurrence(BRC) free time
  * PSA: patient’s PSA(prostate-specific antigen) level
  * Event: two levels, 0 (right censoring), 1 (relapse)
  * Stage: it measures how far the cancer has spread. Nine levels
  * iCluster: seven levels
  * Gleason: Gleason Score. Six levels
  * ECE: extra-capsular extension. Two levels, Yes and No
  * 23 columns for 23 target genes: SLC2A1-14, HK1-3, GCK, OAX1, ENO2, NCAM1, SYP, CHGA, FOLH1, KLK3

* Total number of patients (rows): 343 patients. 111 from Cam dataset, 92 from Stoch dataset and 140 from Taylor dataset.
* Factors that were missing or contained missing values, which were imputed later when fitting some Cox's proportional hazards regression models in Section 4: 
  * Factor—“Age” was unavailable in both Stoch and Taylor datasets
  * Factors—“PSA” and “ECE” were unavailable in Taylor dataset.
  * 2 NAs in “Stage”, 14 NAs in “iCluster” and 4 NAs in “Gleason”. Also, we assumed the deleted missing values in “Time” and “Event” were random distributed.
  
  

## Section 1: Conducting univariate analysis among three datasets

After cleaning and joining the three datasets(Taylor, Cambridge and Stockholm), I conducted univariate analysis among three datasets as a table (Table 1) using compareGroups library. Below is part of the table. 

<img width="500" height="600" src="https://github.com/ensembles4612/Statistical-Analysis-on-Biochemical-Recurrence-Prostate-Cancer-Among-Patients-in-3-Locations/blob/main/table1.png">


## Section 2: Kaplan-Meier Curves

* KM curves among three datasets (Graph 1)

It’s clear from Graph 1 that patients in dataset—Stoch have the highest survival rate before BCR free time turns 90 months. The patients in dataset—Cam have the second highest survival rate until BCR free time turns 60 months, after which the survival rate plumps to 0 when the BCR free time is about 67 months. The patients in dataset—Taylor have the lowest survival time at the beginning and becomes the second in the end. Even though the curves have 95% CIs, it’s still hard to tell if they are significantly different among three datasets.

<img width="600" height="300" src ="https://github.com/ensembles4612/Statistical-Analysis-on-Biochemical-Recurrence-Prostate-Cancer-Among-Patients-in-3-Locations/blob/main/graph1.png">

* KM curves among different Gleason Scores (Graph 2)  

From Graph 2, it’s uneasy to tell if the survival rates of patients having different Gleason Scores are significantly different. Also, not enough observations for Gleason Scores are 5 and 10.

<img width="600" height="300" src="https://github.com/ensembles4612/Statistical-Analysis-on-Biochemical-Recurrence-Prostate-Cancer-Among-Patients-in-3-Locations/blob/main/graph2.png">

## Section 3: Stratified Logrank Test among three datasets 

It’s uneasy to see whether the statistically difference exists among three datasets from Graph 1 since the KM curves cross each other and so on, so conducting a logrank test is a good way to investigate if the probabilities of freedom from BCR among three datasets are significantly different.

Also, a cofounder might exist. After investigation, I found “iCluster” was a cofounder. Firstly, “iCluster” was noticed to be related to the survival under a logrank test. Secondly, it was associated with “dataset” under the Pearson's Chi-squared test.

Therefore, I conducted a logrank test stratified by "iCluster" among the three datasets. From the R output below, we can see the survival rates among these three datasets are significantly different after being stratified by “iCluster”.
```
survdiff(formula = Surv(Time, Event) ~ dataset + strata(iCluster), 
    data = datacamtaylorstock, rho = 0)

n=329, 14 observations deleted due to missingness.

                 N Observed Expected (O-E)^2/E (O-E)^2/V
dataset=Cam    104       18     13.1     1.840      2.83
dataset=stoch   85       41     51.1     2.010      7.73
dataset=taylor 140       36     30.8     0.889      2.21

 Chisq= 8  on 2 degrees of freedom, p= 0.02 
```

## Section 4: Cox’s proportional hazards regression model

**Imputation of missing values**: 

Since “Age”, “PSA” and “ECE” were unavailable in all three datasets, I deleted them. Besides, there were still 2 missing values in “Stage”, 14 Missing values in “iCluster” and 4 missing values in “Gleason”. I imputed all the NA’s in “iCluster” and “Gleason” using Random Forest algorithm and delete the 2 rows with NA’s in “Stage”. 

In order to know which demographic factors and genes have influence on BCR free time, I fit several Cox’s proportional hazards regression models using AIC Stepwise method such as the 2 models below:
* **Model 1**: AIC Stepwise Cox’s model using the pooled dataset with imputation: Likelihood ratio test=114.8 on 20 df, p=3e-15
  ```
  Surv(Time, Event, type = "right") ~ iCluster + Gleason + ENO2 + SYP + CHGA + SLC2A1 + SLC2A6 + SLC2A12 + SLC2A13  + SLC2A14 + HK1 # the model after AIC stepwise
  Likelihood ratio test=114.8  on 20 df, p=3e-15
  n= 341, number of events= 98 
  ```
* **Model 2**: adding interaction terms between “dataset” and each gene: Likelihood ratio test=182.5 on 55 df, p=1e-15 
 ```
  Surv(Time, Event, type = "right") ~ iCluster + Gleason + dataset + ENO2 + NCAM1 + SYP + CHGA + FOLH1 + SLC2A1 +  SLC2A3 + SLC2A4 + SLC2A5 + SLC2A6 + SLC2A7 +SLC2A8 + SLC2A10 + SLC2A12 + GCK + HK1 + dataset:NCAM1 + dataset:SYP + dataset:CHGA + dataset:FOLH1 + dataset:SLC2A1 + dataset:SLC2A3 + dataset:SLC2A5 +dataset:SLC2A6 + dataset:SLC2A7 + dataset:SLC2A8 + dataset:SLC2A12 + dataset:GCK + dataset:HK1 #26(13*2) interaction terms stay in the model
  Likelihood ratio test=182.5  on 55 df, p=1e-15
  n= 341, number of events= 98
 ```

## Section 7: Investigating which factors and genes have greater effects using Random Survival Forests

I used Random Survival Forest technique to find which factors and genes have greater effects on survival since it applies to right-censored data. I grew 400 bootstrapped trees and implemented log-rank splitting. 

* Graph below shows the Leave-one-out cross validation error drops fast and then stables at around 29% as the number of trees increases. 

<img width="300" height="300" src="https://github.com/ensembles4612/Statistical-Analysis-on-Biochemical-Recurrence-Prostate-Cancer-Among-Patients-in-3-Locations/blob/main/Random%20Survival%20Forests-trees.png">

* We can see the importance of each variable with ranking below. “Gleason” has the greatest effect on survival followed by genes—FOLH1, SLC2A1, SYP, SLC2A8 and so forth.

<img width="300" height="300" src="https://github.com/ensembles4612/Statistical-Analysis-on-Biochemical-Recurrence-Prostate-Cancer-Among-Patients-in-3-Locations/blob/main/Random%20Survival%20Forests-importance.png">

