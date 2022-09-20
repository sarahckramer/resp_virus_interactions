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

* "resp_interaction_model.c": C code used to run the model, including parameter transformations, observation model, and deterministic skeleton
* "resp_interaction_model.R": Code to read in season- and virus-virus pair-specific data, create a pomp model, and run various model checks
* "functions/functions_flu_RSV.R": Functions to create pomp objects used for running and fitting the model
* "functions/test_code.R": Contains various functions used to check that model code is behaving as expected
* "functions/setup_global_likelihood.R": Code to read in data from all seasons and load functions for evaluating the global log-likelihood

Fitting the interaction model to data
-------------------------------------

1. Run "fit_traj_matching_round1.R" to obtain initial fits for all season-specific parameters. The parameter "fit_shared" should be set to FALSE. The ranges of values fit will be used to generate starting values for the next round of estimations.
    * Run "process_results/01_check_missing_files_traj_matching_r1.R" to check whether estimation failed for any virus, season, or starting parameter set.
    * Run "process_results/02_compile_results_traj_matching_r1.R" to format and compile individual results into comprehensive output files.
2. Run "fit_traj_matching_round2.R" to fit all parameters, both shared and season-specific, and obtain maximum likelihood estimates. Parameter "search_type" should be set to "round1_CIs"; no_jobs should be set to 125 and sobol_size to 500.
    * Run "checks/check_no_states_below_0.R" to ensure that none of the best-fit parameter values lead to impossible (i.e., negative) values for any of the state variables.
    * If only one parameter set is statistically supported (has log-likelihood value falling within qchisq(0.95, df = 46) / 2 of the log-likelihood of the MLE) at this point, we assume that the MLE has not yet been reached, and run an additional round of fits.
    * Run "get_start_ranges_from_round2.R" to get range of starting parameter values for second run of round 2. The parameter "method" should be set to "perc".
3. Rerun "fit_traj_matching_round2.R", this time with "search_type" set to "round2_CIs".
    * Run "checks/check_no_states_below_0.R" to ensure that none of the best-fit parameter values lead to impossible (i.e., negative) values for any of the state variables.
    * Again, if only one parameter set is supported, repeat the procedure under point 2 above and run again.
4. Once multiple parameter sets are supported, run "get_start_ranges_from_round2.R," this time with "method" set to "ci." Then, perform one final round of model fits by running "fit_traj_matching_round2.R," again with "search_type" set to "round2_CIs".
    * Run "checks/check_no_states_below_0.R" to ensure that none of the best-fit parameter values lead to impossible (i.e., negative) values for any of the state variables.
    * Run "get_MLEs.R" to obtain the maximum likelihood estimates of each parameter and save them.
    * Run "get_start_ranges_from-round2.R", with "method" set to "ci", in order to get parameter start ranges for parametric bootstrapping.
5. Run parametric bootstrapping to get 99% confidence intervals for all parameters.
    * First, run "bootstrap_01_generate_synthetic.R" to generate several synthetic datasets at the MLEs.
    * Next, run "bootstrap_02_fit.R" with no_jobs set to 1 and "sobol_size" set to 10 to fit the model to each synthetic dataset.
    * Finally, run "bootstrap_03_process_and_CIs.R" to compile results and calculate the 99% confidence intervals.
6. Run profile likelihood on theta_lambda1 in order to check that model is converging to the MLE.
    * Run "fit_traj_matching_round2.R" with no_jobs set to 5, "sobol_size" set to 50, and "prof_lik" set to TRUE.
    * Run "process_results/analyze_traj_matching_proflik.R" to determine which values of theta_lambda1 are supported by the analysis and plot the results. Before running, set "res_dir" to the location of the profile likelihood results.

Code to explore data/model fit
------------------------------

* "explore_data/calculate_outbreak_metrics.R":
  * Code to calculate several outbreak metrics for the influenza and RSV data, including attack rates, week of peak activity, and duration of outbreaks
* "explore_data/explore_allcause_mortality_seasonality.R":
  * Code to check for obvious seasonal patterns in all-cause mortality in Hong Kong
* "explore_data/explore_data_smoothness.R":
  * Code to calculate lag-one autocorrelation for all available data
* "process_results/analyze_fit_parameter_vals_round1.R":
  * Code to explore and plot best-fit values from round 1 fits. Before running, set "res_dir" to the location of the results of the "round 1" fits (step 1 above).
* "process_results/analyze_traj_matching_round1.R":
  * Code to explore fit values from round 1, and to check for convergence, correlations between fit parameter values, and agreement between data and simulations at the best-fit values. Before running, set "res_dir" to the location of the results of the "round 1" fits (step 1 above).
* "process_results/analyze_traj_matching_round2.R":
  * Code to explore fit values from round 2; to compare fit values to those from round 1; and to check for convergence, correlations between fit parameter values, and agreement between data and simulations at the best-fit values. Before running, set "res_dir_h1" and "res_dir_b" to the location of the results of the final "round 2" fits (step 4 above) for H1/RSV and B/RSV, respectively, and set "res_dir_round1" to the location of the results of the "round 1" fits.
* "checks/calculate_PAF.R":
  * Code to estimate how the seasonal attack rate of one virus would differ in the absence of the other virus
* "checks/calculate_simulated_AR.R":
  * Code to calculate attack rates of influenza and RSV when model is simulated at the MLE
* "checks/explore_interaction_impact.R":
  * Code to calculate the log-likelihood and run simulations at several combinations of interaction parameter values; intended mostly for exploration

Check for lag on climate data, inclusion of humidity data
---------------------------------------------------------

1. Run "fit_traj_matching_round2.R" with "search_type" set to "round1_CIs" as described above under "Fitting the interaction model to data," but with the parameter "no_ah" set to TRUE.
2. Continue the model fitting process as described above, until the "final" estimates are obtained.
3. Run "compare_AHvnoAH.R" (with lines 13-17 indicating the locations of all relevant results files) to evaluate whether there is a significant difference in model fit between models including and excluding an effect of absolute humidity.

Age-Structured Sensitivity Analysis
-----------------------------------

1. Run "age_structured_SA/m-generate_covariate.R" to generate synthetic, age-structured covariate data (ILI rates and number of tests performed).
2. Run "age_structured_SA/m-run_model.R" to generate synthetic, age-structured case data at the MLE of the model fits to the A(H1N1)-RSV virus-virus pair.
3. To fit the model to the age-structured synthetic data, uncomment line 47 in "fit_traj_matching_round1.R," line 32 in "02_compile_results_traj_matching_r1.R," and line 77 in "fit_traj_matching_round2.R," so that "age_structured" is equal to TRUE. Then, simply follow the same fitting procedure as outlined in numbers 1-4 under "Fitting the interaction model to data" above.

Simulation Study of Vaccine Impact
----------------------------------

1. Run "vaccination_simulation_study/choose_temperate_parameters.R" to select and save initial conditions yielding temperate-like outbreaks.
2. Run "vaccination_simulation_study/run_vaccination_simulation_study.R" to get simulations for all seasons, vaccine coverage levels, and vaccine timings. By default, this is the scenario where vaccination confers the same protection against RSV as does natural influenza infection. To change this, uncomment lines 55/177. To run sensitivity analyses regarding vaccine efficacy against influenza and duration of protection against RSV, lines 15 or 56/178, respectively, can be changed accordingly.
3. Run "vaccination_simulation_study/analyze_vaccination_simulation_study.R" to explore and plot results for both the temperate and subtropical scenarios.

Related code:

* If desired, "vaccination_simulation_study/explore_drivers_of_seasonal_differences.R" can be run to explore potential drivers of the varying results by season.

Generate publication-ready figures for manuscript
-------------------------------------------------

* Figures 1, 3, and 4 from the manuscript can be generated by running "generate_figures.R", which also calculates the $R^2$ between the data and the model simulations at the MLE.
* All supplementary figures (except the schematic for the simulation study of LAIV) can be generated by running "generate_figures_SUPP.R."

Data Sources
------------

* Virologic data: https://www.chp.gov.hk/en/statistics/data/10/641/642/2274.html
* ILI data: https://www.chp.gov.hk/en/static/24015.html
* Population data: https://www.censtatd.gov.hk/en/web_table.html?id=1A
