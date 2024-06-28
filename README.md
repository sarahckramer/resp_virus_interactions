# Respiratory Virus Interactions

Code to infer the strength and duration of interactions between flu and RSV

Directory Structure
-------------------
* src
    * demo
    * format_data
    * explore_data
    * functions
    * checks
    * process_results
    * bootstrap
    * age_structured_SA
    * vaccination_simulation_study

Demo
----

A short demo showing how the model of influenza and RSV cocirculation is fit to data is available as "src/demo/demo_run.R". Synthetic data used for the demo can be found in the same folder as "demo_data.csv". The demo should take about 40 minutes to run, and outputs the fitted parameter values for a limited number of runs, both when the model is fit to data one season at a time, and when the model is fit to all seasons simulataneously.

The demo has been tested on the following operating systems:

* Windows 10

R, the programming language used for this work, can be downloaded at: https://www.r-project.org/

For this work, I used R version 4.4.0. Code was written and run using the RStudio IDE, version 2024.4.1.748, which can be obtained [here](https://posit.co/downloads/).

For the demo, the following R packages must be installed:

* tidyverse (version 2.0.0)
* testthat (version 3.2.1.1)
* pomp (version 5.8)
* nloptr (version 2.0.3)

The "pomp" package has additional dependencies, which vary by operating system. Detailed installation instructions can be found [here](https://kingaa.github.io/pomp/install.html).

Formatting data (Hong Kong)
---------------------------

1. Data from Hong Kong should be downloaded and placed in data/raw/hk/ (see "Data Sources" below). Data on influenza positivity for each year should be saved as csv files with names of the form "hk\_[YEAR]\_flu.csv", while data on positivity of "other respiratory viruses" should be saved as "hk\_[YEAR]\_rsv.csv". Weekly ILI data for each year should be saved as "hk\_[YEAR]\_ili.csv". Data on population size should be saved as "pop_dat_hk.csv".
2. Run "format_HK_data.R" to compile ILI and all virologic data.
3. Run "split_into_seasons_HK.R" to assign seasons and organize data by virus pair.
4. Run "get_climate_data.R" to get normalized temperature and absolute humidity data from Hong Kong (and France, for the simulation study of vaccination).

Formatting data (Canada)
------------------------
1. Data from Canada should be downloaded and placed in data/raw/canada/ (see "Data Sources" below). Data on influenza and RSV positivity for each season should be saved as csv files with names of the form "detection-[SEASON].csv" (so, for example, "detection-201011.csv" for the 2010-11 season); ILI data should be saved as "consultation-[SEASON].csv". Data on population size should be saved as "pop_dat_Canada.csv".
2. Run "format_Canada_data.R" to compile ILI and virologic data.

Model Code
----------

* "resp_interaction_model.c": C code used to run the model, including parameter transformations, observation model, and deterministic skeleton
* "resp_interaction_model.R": Code to read in location- and season-specific data, create a pomp model, and run various model checks
* "functions/functions_flu_RSV.R": Functions to create pomp objects used for running and fitting the model
* "functions/test_code.R": Contains various functions used to check that model code is behaving as expected
* "functions/setup_global_likelihood.R": Code to read in data from all seasons and load functions for evaluating the global log-likelihood

Fitting the interaction model to data
-------------------------------------

1. Run "fit_traj_matching_round1.R" to obtain initial fits for all season-specific parameters. The parameter "sobol_size" should be set to 500, and "search_type" should be set to "broad". For fitting in Hong Kong, "sens" should be set to "main" and "fit_canada" should be set to FALSE. For Canada, set "sens" to "sinusoidal_forcing" and "fit_canada" to TRUE. The ranges of values fit will be used to generate starting values for the next round of estimations. The code will output an RDS file for each of the 500 runs for each season, which contains the best-fit parameter values, the corresponding log-likelihood, and whether the algorithm converged.
    * Run "process_results/01_check_missing_files_traj_matching_r1.R" with "fit_canada" set to the appropriate value to check whether estimation failed for any virus, season, or starting parameter set. The code will print a warning for each missing results file.
    * Run "process_results/02_compile_results_traj_matching_r1.R" to format and compile individual results into comprehensive output files. Ensure that both "sens" and "fit_canada" are set appropriately. Outputs include compiled results from all runs and seasons, only the top fits from each season, and only the maximum likelihood estimates from each season.
2. Run "fit_traj_matching_round2.R" to fit all parameters, both shared and season-specific, and obtain maximum likelihood estimates. Parameter "search_type" should be set to "round1_CIs", "which_round" to 1, "sobol_size" should be set to 500, "int_eff" to "susc", and "prof_lik" to FALSE. Set "sens" and "fit_canada" depending on the location being fit. The parameter "run_parallel" can either be set to TRUE or FALSE, depending on whether runs for several starting parameter sets should be run in parallel. For each run, an RDS file containing the best-fit parameter values, the log-likelihood, and information on convergence is returned, as in step 1. If "run_parallel" is set to TRUE, all results run on a given core will be returned in a single RDS file; the number of cores is controlled by the parameter "no_jobs".
    * Run "get_start_ranges_from_round2.R", with line 10 set to be the location of the current results. This will ensure that none of the best-fit parameter values lead to impossible (i.e., negative) values for any of the state variables, and obtain the range of starting parameter values for the next round of fitting, which are output as an RDS file.
    * If only one parameter set is statistically supported (has log-likelihood value falling within qchisq(0.95, df = # of parameters estimated) / 2 of the log-likelihood of the MLE) at this point, we assume that the MLE has not yet been reached, and run an additional round of fits.
3. Rerun "fit_traj_matching_round2.R", this time with "search_type" set to "round2_CIs" and "which_round" set to 2. Output are as described in step 2.
    * Again, run "get_start_ranges_from_round2.R".
    * Repeat this step until multiple parameter sets are supported, each time increasing "which_round" by one.
4. Once multiple parameter sets are supported, perform one final round of model fits by running "fit_traj_matching_round2.R," again with "search_type" set to "round2_CIs". Again, output are as in step 2.
    * Run "get_start_ranges_from_round2.R". This will calculate the parameter start ranges for parametric bootstrapping, as well as save the maximum likelihood estimates for all parameters.
5. Run parametric bootstrapping to get 95% confidence intervals for all parameters.
    * First, run "bootstrap_01_generate_synthetic.R" to generate several synthetic datasets at the MLEs. Both "sens" and "fit_canada" should be set appropriately. Output is an RDS file containing a list of synthetic case counts for influenza and RSV; 500 synthetic data sets are generated for each season.
    * Next, run "bootstrap_02_fit.R" with "sobol_size" set to 10 and "final_round" set to the number of rounds of fitting run in steps 2-4. All other parameters should be set to the same values as in step 2. This will fit the model to each synthetic dataset. As in step 2, the resulting RDS files contain the best-fit parameter values, the log-likelihoods, and information on covergence for each run.
    * Finally, run "bootstrap_03_process_and_CIs.R" to compile results and calculate the 95% confidence intervals. Line 12 should be the location of the results. Outputs are csv files containing information on the maximum likelihood estimates and 95% confidence intervals for all model parameters.
6. Run profile likelihood on $\theta_{\lambda1}$ in order to check that model is converging to the MLE.
    * Run "fit_traj_matching_round2.R" with "sobol_size" set to 100, "search_type" set to "round2_CIs", "final_round" set to the number of rounds of fitting run in steps 2-4, and "prof_lik" set to TRUE. Outputs contain the same information as in step 2 above.

Code to explore data/model fit
------------------------------

* "explore_data/calculate_outbreak_metrics.R":
  * Code to calculate several outbreak metrics for the influenza and RSV data, including attack rates, week of peak activity, and duration of outbreaks
* "explore_data/explore_allcause_mortality_seasonality.R":
  * Code to check for obvious seasonal patterns in all-cause mortality in Hong Kong
* "explore_data/explore_data_smoothness.R":
  * Code to calculate lag-one autocorrelation for all available data
* "explore_data/stat_quantify_similarity.R":
  * Code to calculate the extent of the correlation between rates of different pathogens
* "process_results/analyze_traj_matching_round1.R":
  * Code to explore fit values from round 1, and to check for convergence, correlations between fit parameter values, and agreement between data and simulations at the best-fit values. Before running, set "res_dir" to the location of the results of the "round 1" fits (step 1 above), and set "sens" and "fit_canada" appropriately.
* "process_results/analyze_traj_matching_round2.R":
  * Code to explore fit values from round 2; to compare fit values to those from round 1; and to check for convergence, correlations between fit parameter values, and agreement between data and simulations at the best-fit values. Before running, set "res_dir" to the location of the results of the final "round 2" fits (step 4 above), and set "res_dir_round1" to the location of the results of the "round 1" fits. The parameter "fit_canada" in line 13 should also be set appropriately.
* "process_results/get_start_ranges_sinusoidal.R":
  * Code to obtain reasonable starting parameter ranges for the amplitude of sinusoidal seasonal forcing, based on fits of sine waves through the extent of climate forcing in Hong Kong at the MLE.
* "process_results/compare_sensitivity.R":
  * Code to compare the parameter estimates and log-likelihoods of the main model and several sensitivity analyses. Lines 13-24 should be set to the locations of the various results.
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
* "functions/functions_evaluate_res.R":
  * Contains a function to obtain deterministic or stochastic simulations from the model at the MLE

Age-Structured Sensitivity Analysis
-----------------------------------

1. Run "age_structured_SA/m-generate_covariate.R" to generate synthetic, age-structured covariate data (ILI rates and number of tests performed).
2. Run "age_structured_SA/m-run_model.R" to generate synthetic, age-structured case data at the MLE of the model fits.
3. To fit the model to the age-structured synthetic data, uncomment line 56 in "fit_traj_matching_round1.R," line 37 in "02_compile_results_traj_matching_r1.R," and line 110 in "fit_traj_matching_round2.R," so that "age_structured" is equal to TRUE. Then, simply follow the same fitting procedure as outlined in numbers 1-4 under "Fitting the interaction model to data" above.

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

1. Run "vaccination_simulation_study/run_vaccination_simulation_study.R" to get simulations for all seasons, vaccine coverage levels, and vaccine timings. When "sens" is set to "main", this is the scenario where vaccination confers the same protection against RSV as does natural influenza infection, as fit in Hong Kong. To change this, set "sens_sim" to "vacc_can", "deltavacc1month", "deltavacc6months", "vacceff60", and "vacceff95", in turn. This will run sensitivity analyses with different values for the impact of vaccination on RSV susceptibility, the duration of the effect on RSV susceptibility, and the efficacy of vaccination against influenza. By default, the interaction parameters for natural infection are taken from the MLE fit in Hong Kong; to instead use values as fit in Canada, set "sens_sim" to "fit_can". Results for each season and vaccination coverage level are saved as RDS files.
2. Run "vaccination_simulation_study/analyze_vaccination_simulation_study.R" to explore and plot results for both the temperate and subtropical scenarios. Compiled results for the temperate and subtropical scenarios are output as RDS files; plots of the average and season-specific impact of LAIV are saved as SVG and PDF files, respectively.

Generate publication-ready figures for manuscript
-------------------------------------------------

* Figures 1, 2b-c, 3, and 4 from the manuscript can be generated by running "generate_figures.R", which also calculates the $R^2$ between the data and the model simulations at the MLE.
* All supplementary figures (except the schematic for the simulation study of LAIV) can be generated by running "generate_figures_SUPP.R."

Data Sources
------------

* Hong Kong:
  * Virologic data: https://www.chp.gov.hk/en/statistics/data/10/641/642/2274.html
  * ILI data: https://www.chp.gov.hk/en/static/24015.html
  * Population data: https://www.censtatd.gov.hk/en/web_table.html?id=1A

* Canada:
  * Virologic and ILI data: https://search.open.canada.ca/opendata/?od-search-portal=Open%20Data&search_text=fluwatch (separate pages exist for each season)
  * Population data: https://www150.statcan.gc.ca/t1/tbl1/en/tv.action?pid=1710000901
