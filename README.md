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

1. Data from Hong Kong should be downloaded and placed in data/raw/hk/ (see "Data Sources" below). Data on influenza positivity for each year should be saved as csv files with names of the form "hk\_[YEAR]\_flu.csv", while data on positiviy of "other respiratory viruses" should be saved as "hk\_[YEAR]\_rsv.csv". Weekly ILI data for each year should be saved as "hk\_[YEAR]\_ili.csv". Data on population size should be saved as "pop\_dat\_hk.csv".
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

1. Run "fit_traj_matching_round1.R" to obtain initial fits for all season-specific parameters. The parameter "sobol_size" should be set to 500, "search_type" should be set to "broad", "sens" should be set to "main", and "fit_canada" should be set to FALSE. The ranges of values fit will be used to generate starting values for the next round of estimations.
    * Run "process_results/01_check_missing_files_traj_matching_r1.R" with "fit_canada" set to FALSE to check whether estimation failed for any virus, season, or starting parameter set.
    * Run "process_results/02_compile_results_traj_matching_r1.R" to format and compile individual results into comprehensive output files. Again, "sens" should be set to "main" and "fit_canada" to FALSE.
2. Run "fit_traj_matching_round2.R" to fit all parameters, both shared and season-specific, and obtain maximum likelihood estimates. Parameter "search_type" should be set to "round1_CIs", "which_round" to 1, "sobol_size" should be set to 500, "int_eff" to "susc", "prof_lik" to FALSE, "sens" to "main", and "fit_canada" to FALSE.
    * Run "get_start_ranges_from_round2.R", with line 10 set to be the location of the current results. This will ensure that none of the best-fit parameter values lead to impossible (i.e., negative) values for any of the state variables, and obtain the range of starting parameter values for the next round of fitting.
    * If only one parameter set is statistically supported (has log-likelihood value falling within qchisq(0.95, df = 46) / 2 of the log-likelihood of the MLE) at this point, we assume that the MLE has not yet been reached, and run an additional round of fits.
3. Rerun "fit_traj_matching_round2.R", this time with "search_type" set to "round2_CIs" and "which_round" set to 2.
    * Again, run "get_start_ranges_from_round2.R".
    * Repeat this step until multiple parameter sets are supported, each time increasing "which_round" by one.
4. Once multiple parameter sets are supported, perform one final round of model fits by running "fit_traj_matching_round2.R," again with "search_type" set to "round2_CIs".
    * Run "get_start_ranges_from_round2.R". This will calculate the parameter start ranges for parametric bootstrapping, as well as save the maximum likelihood estimates for all parameters.
5. Run parametric bootstrapping to get 95% confidence intervals for all parameters.
    * First, run "bootstrap_01_generate_synthetic.R" to generate several synthetic datasets at the MLEs.
    * Next, run "bootstrap_02_fit.R" with "sobol_size" set to 10 and "final_round" set to the number of rounds of fitting run in steps 2-4. All other parameters should be set to the same values as in step 2. This will fit the model to each synthetic dataset.
    * Finally, run "bootstrap_03_process_and_CIs.R" to compile results and calculate the 99% confidence intervals. Line 12 should be the location of the results.
6. Run profile likelihood on $\theta_{\lambda1}$ and $\theta_{\lambda2}$ in order to check that model is converging to the MLE.
    * Run "fit_traj_matching_round2.R" with "sobol_size" set to 100, "search_type" set to "round2_CIs", "final_round" set to the number of rounds of fitting run in steps 2-4, and "prof_lik" set to TRUE.
    * Comment out line 70, uncomment line 71, and repeat the above step.

Code to explore data/model fit
------------------------------

* "explore_data/calculate_outbreak_metrics.R":
  * Code to calculate several outbreak metrics for the influenza and RSV data, including attack rates, week of peak activity, and duration of outbreaks
* "explore_data/explore_allcause_mortality_seasonality.R":
  * Code to check for obvious seasonal patterns in all-cause mortality in Hong Kong
* "explore_data/explore_data_smoothness.R":
  * Code to calculate lag-one autocorrelation for all available data
* "process_results/analyze_traj_matching_round1.R":
  * Code to explore fit values from round 1, and to check for convergence, correlations between fit parameter values, and agreement between data and simulations at the best-fit values. Before running, set "res_dir" to the location of the results of the "round 1" fits (step 1 above).
* "process_results/analyze_traj_matching_round2.R":
  * Code to explore fit values from round 2; to compare fit values to those from round 1; and to check for convergence, correlations between fit parameter values, and agreement between data and simulations at the best-fit values. Before running, set "res_dir" sto the location of the results of the final "round 2" fits (step 4 above), and set "res_dir_round1" to the location of the results of the "round 1" fits.
* "process_results/compare_sensitivity.R":
  * Code to compare the parameter estimates and log-likelihoods of the main model and several sensitivity analyses. Lines 13-23 should be set to the locations of the various results.
* "checks/calculate_PAF.R":
  * Code to estimate how the seasonal attack rate of one virus would differ in the absence of the other virus
* "checks/calculate_simulated_AR.R":
  * Code to calculate attack rates of influenza and RSV when model is simulated at the MLE
* "checks/calculate_simulated_metrics.R":
  * Code to compare predicted and observed outbreak metrics for influenza and RSV at the MLE
* "checks/check_determ_vs_stoch.R":
  * Code to expore the extent to which stochasticity leads to substantial variations from deterministic trajectory at the MLE
* "checks/stat_quantify_similarity.R":
  * Code to quantify the correlation between observed rates of influenza, RSV, and rhinovirus

Age-Structured Sensitivity Analysis
-----------------------------------

1. Run "age_structured_SA/m-generate_covariate.R" to generate synthetic, age-structured covariate data (ILI rates and number of tests performed).
2. Run "age_structured_SA/m-run_model.R" to generate synthetic, age-structured case data at the MLE of the model fits.
3. To fit the model to the age-structured synthetic data, uncomment line 53 in "fit_traj_matching_round1.R," line 36 in "02_compile_results_traj_matching_r1.R," and line 82 in "fit_traj_matching_round2.R," so that "age_structured" is equal to TRUE. Then, simply follow the same fitting procedure as outlined in numbers 1-4 under "Fitting the interaction model to data" above.

Additional Sensitivity Analyses
-------------------------------

Various sensitivity analyses can be conducted by running steps 2-5 above with "sens" set to:

* "sinusoidal_forcing": Fits a model using a sine wave to capture seasonal changes in the force of infection of influenza and RSV, rather than explicitly including climate data. Rough start ranges for the parameters $b_1$ and $b_2$, describing the extent to which the strength of forcing varies over the year, can be obtained by running "process_results/get_start_ranges_sinusoidal.R", which fits sine waves to the force of infection for both viruses at the MLE from the main analysis. The resulting values are then used in lines 310-311 in "fit_traj_matching_round2.R".
* "no_ah": Fits the model exluding an effect of absolute humidity on the force of infection of influenza and RSV, such that both are modulated by temperature only.
* "no_rsv_immune": Fits a model assuming, as in Waterlow et al. (2022), that the entire model population is susceptible to RSV at the beginning of each season. Before running steps 2-5 above, step 1 should also be repeated, this time with "sens" also set to "no_rsv_immune".
* "no_int": Fits a model in which no interaction occurs between influenza and RSV.
* "h3_covar": Fits a model allowing H3N2 incidence to modulate susceptibility to RSV. In order to explore different lags on H3N2 incidence, lines 115-125 in "resp_interaction_model.R" can be updated.
* "less_circ_h3": Fits the model only for those seasons with little H3N2 circulation (2017-18 and 2018-19).
* "rhino_covar": Fits a model allowing rhinovirus incidence to modulate susceptibility to influenza.

Simulation Study of Vaccine Impact
----------------------------------

1. Run "vaccination_simulation_study/choose_temperate_parameters.R" to select and save initial conditions yielding temperate-like outbreaks.
2. Run "vaccination_simulation_study/run_vaccination_simulation_study.R" to get simulations for all seasons, vaccine coverage levels, and vaccine timings. When "sens" is set to "main", this is the scenario where vaccination confers the same protection against RSV as does natural influenza infection. To change this, set "sens" to "thetalambdavacc0.50", "deltavacc1month", "deltavacc6months", "vacceff60", and "vacceff95", in turn. This will run sensitivity analyses with different values for the impact of vaccination on RSV susceptibility, the duration of the effect on RSV susceptibility, and the efficacy of vaccination against influenza.
3. Run "vaccination_simulation_study/analyze_vaccination_simulation_study.R" to explore and plot results for both the temperate and subtropical scenarios.

Generate publication-ready figures for manuscript
-------------------------------------------------

* Figures 1, 3, and 4 from the manuscript can be generated by running "generate_figures.R", which also calculates the $R^2$ between the data and the model simulations at the MLE.
* All supplementary figures (except the schematic for the simulation study of LAIV) can be generated by running "generate_figures_SUPP.R."

Data Sources
------------

* Virologic data: https://www.chp.gov.hk/en/statistics/data/10/641/642/2274.html
* ILI data: https://www.chp.gov.hk/en/static/24015.html
* Population data: https://www.censtatd.gov.hk/en/web_table.html?id=1A
