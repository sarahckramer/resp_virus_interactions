# Respiratory Virus Interactions

Code to infer the strength and duration of interactions between flu and RSV

Directory Structure
-------------------
* src
    * format_data
    * explore_data
    * functions
    * checks

Formatting and assessing data
-----------------------------

Model Code
----------

Fitting the interaction model to data
-------------------------------------

1. Run "fit_traj_matching_round1.R" to obtain initial fits for all season-specific parameters. The parameter "fit_shared" should be set to FALSE. The ranges of values fit will be used to generate starting values for the next round of estimations.
    * Run "01_check_missing_files_traj_matching_r1.R" to check whether estimation failed for any virus, season, or starting parameter set. If so, these can be rerun.
    * Run "02_compile_results_traj_matching_r1.R" to format and compile individual results into comprehensive output files.
2. Run "fit_traj_matching_round2.R" to fit all parameters, both shared and season-specific, and obtain maximum likelihood estimates. Parameter "search_type" should be set to "round1_CIs".
    * Run "check_no_states_below_0.R" to ensure that none of the best-fit parameter values lead to impossible (i.e., negative) values for any of the state variables.
    * At this point, very few (1-2) fit parameter sets are statistically supported (have log-likelihood value falling within qchisq(0.95, df = 46) / 2 of the log-likelihood of the MLE). This suggests that the MLE has not yet been reached, and that a second run of round 2 may be necessary.
    * Run "get_start_ranges_from_round2.R" to get range of starting parameter values for second run of round 2. The parameter "method" should be set to "perc".
3. Rerun "fit_traj_matching_round2.R", this time with "search_type" set to "round2_CIs".
    * Run "check_no_states_below_0.R" to ensure that none of the best-fit parameter values lead to impossible (i.e., negative) values for any of the state variables.
    * Run "get_MLEs.R" to obtain the maximum likelihood estimates of each parameter and save them.
    * Run "get_start_ranges_from-round2.R", this time with "method" set to "ci", in order to get parameter start ranges for parametric bootstrapping.
4. Run parametric bootstrapping to get 99% confidence intervals for all parameters.
    * First, run "bootstrap_01_generate_synthetic.R" to generate several synthetic datasets at the MLEs.
    * Next, run "bootstrap_02_fit.R" to fit the model to each synthetic dataset.
    * Finally, run "bootstrap_03_process_and_CIs.R" to compile results and calculate the 99% confidence intervals.

Model Checks
------------

Check for lag on climate data, inclusion of humidity data
---------------------------------------------------------

Data Sources
------------
