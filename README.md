# Hip Fracture Prediction in Older Chinese Adults

This repository contains the R code used for the development and internal validation of a multidomain prediction model for incident hip fracture in older Chinese adults using data from the China Health and Retirement Longitudinal Study (CHARLS).

The repository accompanies the manuscript:

**Development and internal validation of a multidomain prediction model for incident hip fracture in older Chinese adults: a prospective cohort study**

## Study overview

This study used the 2011 CHARLS survey as baseline, with follow-up assessments in 2013, 2015, 2018, and 2020.

Participants were eligible if they were aged ≥60 years and had no history of hip fracture at baseline. The final analytic cohort included 6,373 participants, among whom 324 developed incident hip fracture during follow-up.

The final prediction model included 10 baseline predictors:

- Age
- History of falls
- Activities of daily living (ADL) limitations
- Arthritis
- CESD-10 depressive symptom score
- Hearing status
- Global cognition score
- 2.5-m walking time
- Peak expiratory flow
- Dental status

## Analytical workflow

The analysis was conducted in R and included the following major steps:

1. Cohort construction and data cleaning
2. Candidate predictor definition and harmonization
3. Multiple imputation by chained equations
4. Design matrix construction and Elastic Net/LASSO tuning
5. Predictor stability assessment and final model specification
6. Cox proportional hazards modelling and proportional hazards diagnostics
7. Internal validation using bootstrap optimism correction
8. Evaluation of discrimination, calibration, prediction error, and clinical utility
9. Exploratory geographic hold-out validation
10. Sensitivity and subgroup analyses

## Repository structure

The scripts should be run approximately in the following order:

| Script | Description |
|---|---|
| `STEP1.Data Cleaning & Cohort Definition.R` | Construction of the baseline and prospective analytic cohorts |
| `STEP2.Candidate Predictor Refining & Data Harmonization.R` | Definition, recoding, and harmonization of candidate predictors |
| `STEP3.Multiple Imputation by Chained Equations (MICE).R` | Multiple imputation of missing baseline predictors |
| `STEP4.Design Matrix Construction & Elastic Net_LASSO Hyperparameter Tuning.R` | Design matrix construction and penalized Cox model tuning |
| `STEP5.Stability analysis, model scale.R` | Predictor stability assessment and definition of the parsimonious final predictor set |
| `STEP6.Cox model PH assumption diagnostics and bootstrap.R` | Final Cox model estimation, proportional hazards diagnostics, and bootstrap procedures |
| `STEP7.Inner CV Selection & Internal Validation.R` | Cross-validation and internal validation analyses |
| `STEP8.IPCW Brier IPA Evaluation and Decision Curve Analysis.R` | Time-dependent prediction error, IPA, and decision curve analysis |
| `STEP9.Geographic Internal-External.R` | Exploratory geographic hold-out analysis |
| `STEP10.Sensitivity_and_Subgroup_Performance.R` | Sensitivity analyses and subgroup performance evaluation |

## Model evaluation

Model performance was evaluated using complementary measures including:

- Harrell's C-index
- Time-dependent AUC at 3, 5, and 7 years
- Calibration intercept and slope
- Brier score
- Index of Prediction Accuracy (IPA)
- Decision curve analysis

Internal validation of the locked final model was performed using bootstrap optimism correction with 500 bootstrap resamples within each of 10 multiply imputed datasets.

Sensitivity analyses included:

- A 2-year landmark analysis
- Alternative event-time assignments within the observed fracture-reporting interval
- An interval-censored proportional hazards analysis
- Subgroup analyses by age, sex, and place of residence

An exploratory geographic hold-out analysis was additionally conducted using Sichuan Province as the held-out sample.

## Data availability

CHARLS data are not redistributed in this repository.

Researchers may apply for access to the original CHARLS data through the official CHARLS data platform:

http://charls.pku.edu.cn/

Users of this repository are responsible for obtaining appropriate permission to access and use CHARLS data in accordance with the CHARLS data-use policies.

Because the original participant-level data cannot be distributed through this repository, local file paths and data-import steps may need to be adapted before running the scripts.

## Reproducibility

The scripts document the analytical workflow used in the study, including cohort construction, preprocessing, multiple imputation, predictor selection, model development, internal validation, performance assessment, sensitivity analyses, subgroup analyses, and visualization.

No individual-level CHARLS data are included in this repository.

## Software

Analyses were performed in R.

Major R packages used in the analytical workflow include packages for:

- survival analysis
- multiple imputation
- penalized regression
- time-dependent prediction performance
- calibration
- bootstrap validation
- decision curve analysis

Exact package requirements can be identified from the `library()` calls within the individual analysis scripts.

## Citation

If you use or adapt this code, please cite the accompanying article:

> Kang CW, Wu LX, Gong JB, Dai SY, Yu TB.  
> **Development and internal validation of a multidomain prediction model for incident hip fracture in older Chinese adults: a prospective cohort study.**  
> Manuscript under review.

The citation will be updated after publication.

## License

This repository contains research code associated with the accompanying manuscript.  
Please contact the authors before substantial reuse or redistribution of the code.

## Contact

For questions regarding the study or analytical code, please contact the corresponding authors through the contact information provided in the manuscript.
