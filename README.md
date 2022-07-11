# Respiratory Virus Interactions

Code to infer the strength and duration of interactions between flu and RSV

Directory Structure
-------------------
* src
    * format_data
    * explore_data
    * functions
    * checks
    * process_results
    * bootstrap
    * age_structured_SA
    * vaccination_simulation_study

Formatting data
-----------------------------

1. Data from Hong Kong should be downloaded and placed in data/raw/hk/ (see "Data Sources" below). Data on influenza positivity for each year should be saved as csv files with names of the form "hk_[YEAR]_flu.csv," while data on positiviy of "other respiratory viruses" should be saved as "hk_[YEAR]_rsv.csv." Weekly ILI data for each year should be saved as "hk_[YEAR]_ili.csv." Data on population size should be saved as "pop_dat_hk.csv".
2. Run "format_HK_data.R" to compile ILI and all virologic data.
3. Run "split_into_seasons_HK.R" to assign seasons and organize data by virus pair.
4. Run "get_climate_data.R" to get normalized temperature and absolute humidity data from Hong Kong (and France, for the simulation study of vaccination).

Model Code
----------

* "resp_interaction_model.c":
* "resp_interaction_model.R":
* "functions/functions_flu_RSV.R":
* "functions/test_code.R":
* "functions/setup_global_likelihood.R"

Fitting the interaction model to data
-------------------------------------

1. Run "fit_traj_matching_round1.R" to obtain initial fits for all season-specific parameters. The parameter "fit_shared" should be set to FALSE. The ranges of values fit will be used to generate starting values for the next round of estimations.
    * Run "process_results/01_check_missing_files_traj_matching_r1.R" to check whether estimation failed for any virus, season, or starting parameter set.
    * Run "process_results/02_compile_results_traj_matching_r1.R" to format and compile individual results into comprehensive output files.
2. Run "fit_traj_matching_round2.R" to fit all parameters, both shared and season-specific, and obtain maximum likelihood estimates. Parameter "search_type" should be set to "round1_CIs".
    * Run "checks/check_no_states_below_0.R" to ensure that none of the best-fit parameter values lead to impossible (i.e., negative) values for any of the state variables.
    * At this point, very few (1-2) fit parameter sets are statistically supported (have log-likelihood value falling within qchisq(0.95, df = 46) / 2 of the log-likelihood of the MLE). This suggests that the MLE has not yet been reached, and that a second run of round 2 may be necessary.
    * Run "get_start_ranges_from_round2.R" to get range of starting parameter values for second run of round 2. The parameter "method" should be set to "perc".
3. Rerun "fit_traj_matching_round2.R", this time with "search_type" set to "round2_CIs".
    * Run "checks/check_no_states_below_0.R" to ensure that none of the best-fit parameter values lead to impossible (i.e., negative) values for any of the state variables.
    * Run "get_MLEs.R" to obtain the maximum likelihood estimates of each parameter and save them.
    * Run "get_start_ranges_from-round2.R", this time with "method" set to "ci", in order to get parameter start ranges for parametric bootstrapping.
4. Run parametric bootstrapping to get 99% confidence intervals for all parameters.
    * First, run "bootstrap_01_generate_synthetic.R" to generate several synthetic datasets at the MLEs.
    * Next, run "bootstrap_02_fit.R" to fit the model to each synthetic dataset.
    * Finally, run "bootstrap_03_process_and_CIs.R" to compile results and calculate the 99% confidence intervals.

Code to explore data/model fit
------------------------------

* "explore_data/calculate_outbreak_metrics.R":
* "explore_data/explore_allcause_mortality_seasonality.R":
* "process_results/analyze_fit_parameter_vals_round1.R":
* "process_results/analyze_traj_matching_round1.R":
* "process_results/analyze_traj_matching_round2.R":
* "process_results/compare_AHvnoAH_and_lags.R":
* "checks/calculate_simulated_AR.R":
* "checks/explore_interaction_impact.R":

Check for lag on climate data, inclusion of humidity data
---------------------------------------------------------

Age-Structured Sensitivity Analysis
-----------------------------------

Simulation Study of Vaccine Impact
----------------------------------

Data Sources
------------

* Virologic data: https://www.chp.gov.hk/en/statistics/data/10/641/642/2274.html
* ILI data: https://www.chp.gov.hk/en/static/24015.html
* Population data: https://data.gov.hk/en-data/dataset/hk-censtatd-tablechart-popn/resource/75ca854e-d06d-48f4-a565-e136f0d46db3
