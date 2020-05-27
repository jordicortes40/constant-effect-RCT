# Constant effect in randomized clinical trials with quantitative outcome. Bibliographic review.

**Under construction**

This code was used to conduct the work explained in the Phd dissertation: *Constant effect in randomized clinical trials with quantitative outcome. Bibliographic review.*

## Scripts

This section contains a brief explanation of the scripts contained in the *code* folder:

### Main scripts

- _**read_data**_. Install and load libraries. Read and clean the data.
- _**descriptive**_. Descriptive of the dataset containing information on 208 randomized clinical trials.
- _**MA_main_analysis_rma**_. Main analysis based on random effects models.
- _**SA_I_heuristic**_. Sensitivity analysis I. Heuristic procedure based on removing studies one by one to achieve a negligible heterogeneity.
- _**SA_II_simulation**_. Senitivity analysis II. Simulation study to determine the conditions under wich, the etimations of the random effects model would be obtained.
- _**SA_III_usual_tests**_. Sensitivity analysis III. Based on classic tests to compare variances.
- _**SA_IV_mixture_distribution**_. Senssitivity analysis IV. Based on fitting a mixture distribution to the p-values derived from the previous analysis.
- -**summary_table.R**_. Main results of all previous analysis.


### Ancilliary scripts

- _**functions**_. Ancilliary functions used in the main scripts.
- _**rma_models**_. Fit of all models for the main analysis.
- _**subgroups**_. Subgroup analyses.
- _**rma_models_reduced_data**_. Fit of all models for the sensitivity analysis II.

## Main results

According to the performed analyses and the bibliograpic review, the percentage of studies with non-constant effect in randomized controlled trials with quantitative outcome is around 20%.

