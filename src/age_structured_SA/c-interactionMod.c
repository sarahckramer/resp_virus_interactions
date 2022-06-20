// Implementation of age-structured model of flu-RSV interaction
// Extension on Sarah's homogeneous model, with other parameters taken from Waterlow et al. Epidemics 2021
// Model equations are in pdf document written by Sarah 

//start_globs
// Expit transform for parameters constrained in interval [a,b]
static double expitCons(double x, double a, double b) {
double out = (a + b * exp(x)) / (1.0 + exp(x));
if(ISNAN(out)) out = (b + a * exp(-x)) / (1.0 + exp(-x)); // If x=+Inf, must return b
return out;
}

// Logit transform for parameters constrained in interval [a,b]
static double logitCons(double x, double a, double b) {
x = (x <= a) ? a : (x >= b ? b : x);
double out = log((x - a) / (b - x));
return out;
}

// Probability of transition for event of rate r during time step delta_t
// p = 1 - exp(-r * delta_t)
static double pTrans(double r, double delta_t) {

// r: event (instantaneous rate)
// delta_t: time step
// Returns: probability of transition
double out = 1.0 - exp(-r * delta_t);
return out;
}
//end_globs

//start_rinit

//Pointers to state variables
double *X_SS_vec = (double *) &X_SS_1;
double *X_IS_vec = (double *) &X_IS_1;
double *X_TS_vec = (double *) &X_TS_1;
double *X_RS_vec = (double *) &X_RS_1;
double *X_SI_vec = (double *) &X_SI_1;
double *X_II_vec = (double *) &X_II_1;
double *X_TI_vec = (double *) &X_TI_1;
double *X_RI_vec = (double *) &X_RI_1;
double *X_ST_vec = (double *) &X_ST_1;
double *X_IT_vec = (double *) &X_IT_1;
double *X_TT_vec = (double *) &X_TT_1;
double *X_RT_vec = (double *) &X_RT_1;
double *X_SR_vec = (double *) &X_SR_1;
double *X_IR_vec = (double *) &X_IR_1;
double *X_TR_vec = (double *) &X_TR_1;
double *X_RR_vec = (double *) &X_RR_1;
double *H1tot_vec = (double *) &H1tot_1;
double *H2tot_vec = (double *) &H2tot_1;
double *H1_vec = (double *) &H1_1;
double *H2_vec = (double *) &H2_1;

// Pointers to age-specific parameters
double *N_vec = (double *) &N_1; // age-specific population size
int i; 

// Calculate initial values of state variables 
// Initial fractions assumed constant across age groups
for(i = 0; i < nA; i++) {
 	X_SS_vec[i] = nearbyint((1.0 - I10 - I20 - R10 - R20 - R120) * N_vec[i]);
    X_IS_vec[i] = nearbyint(I10 * N_vec[i]);
    X_TS_vec[i] = 0;
    X_RS_vec[i] = nearbyint(R10 * N_vec[i]);
    X_SI_vec[i] = nearbyint(I20 * N_vec[i]);
    X_II_vec[i] = 0;
    X_TI_vec[i] = 0;
    X_RI_vec[i] = 0;
    X_ST_vec[i] = 0;
    X_IT_vec[i] = 0;
    X_TT_vec[i] = 0;
    X_RT_vec[i] = 0;
    X_SR_vec[i] = nearbyint(R20 * N_vec[i]);
    X_IR_vec[i] = 0;
    X_TR_vec[i] = 0;
    X_RR_vec[i] = nearbyint(R120 * N_vec[i]);
    H1tot_vec[i] = 0;
    H2tot_vec[i] = 0;
    H1_vec[i] = 0; 
    H2_vec[i] = 0;
 }
//end_rinit

//start_skel

int i, j;
double delta2; 

//Pointers to state variables
double *X_SS_vec = (double *) &X_SS_1;
double *X_IS_vec = (double *) &X_IS_1;
double *X_TS_vec = (double *) &X_TS_1;
double *X_RS_vec = (double *) &X_RS_1;
double *X_SI_vec = (double *) &X_SI_1;
double *X_II_vec = (double *) &X_II_1;
double *X_TI_vec = (double *) &X_TI_1;
double *X_RI_vec = (double *) &X_RI_1;
double *X_ST_vec = (double *) &X_ST_1;
double *X_IT_vec = (double *) &X_IT_1;
double *X_TT_vec = (double *) &X_TT_1;
double *X_RT_vec = (double *) &X_RT_1;
double *X_SR_vec = (double *) &X_SR_1;
double *X_IR_vec = (double *) &X_IR_1;
double *X_TR_vec = (double *) &X_TR_1;
double *X_RR_vec = (double *) &X_RR_1;
double *H1tot_vec = (double *) &H1tot_1;
double *H2tot_vec = (double *) &H2tot_1;
double *H1_vec = (double *) &H1_1;
double *H2_vec = (double *) &H2_1;

//Pointers to derivatives
double *DX_SS_vec = (double *) &DX_SS_1;
double *DX_IS_vec = (double *) &DX_IS_1;
double *DX_TS_vec = (double *) &DX_TS_1;
double *DX_RS_vec = (double *) &DX_RS_1;
double *DX_SI_vec = (double *) &DX_SI_1;
double *DX_II_vec = (double *) &DX_II_1;
double *DX_TI_vec = (double *) &DX_TI_1;
double *DX_RI_vec = (double *) &DX_RI_1;
double *DX_ST_vec = (double *) &DX_ST_1;
double *DX_IT_vec = (double *) &DX_IT_1;
double *DX_TT_vec = (double *) &DX_TT_1;
double *DX_RT_vec = (double *) &DX_RT_1;
double *DX_SR_vec = (double *) &DX_SR_1;
double *DX_IR_vec = (double *) &DX_IR_1;
double *DX_TR_vec = (double *) &DX_TR_1;
double *DX_RR_vec = (double *) &DX_RR_1;
double *DH1tot_vec = (double *) &DH1tot_1;
double *DH2tot_vec = (double *) &DH2tot_1;
double *DH1_vec = (double *) &DH1_1;
double *DH2_vec = (double *) &DH2_1;

// Pointers to age-specific parameters
double *N_vec = (double *) &N_1; // age-specific population size
double *tau_vec = (double *) &tau_1; // age-specific susceptibility to RSV

// Other age-specific parameters
double lambda1_vec[nA]; // Age-specific force of infection of flu
double lambda2_vec[nA]; // Age-specific force of infection of RSV
double p1_vec[nA]; // Age-specific prevalence of flu
double p2_vec[nA]; // Age-specific prevalence of RSV

// Re-create contact matrix
double CM_mat[nA][nA]; 
double *CM_vec = (double *) &CM_1; // Contact matrix stacked as vector (line by line), size nA*nA
for(i = 0; i < nA; i++) {
    for(j = 0; j < nA; j++) {
        CM_mat[i][j] = CM_vec[nA * i + j];
    }
}

// Calculate prevalence of infection
for(i = 0; i < nA; i++) {
    p1_vec[i] = (X_IS_vec[i] + X_II_vec[i] + X_IT_vec[i] + X_IR_vec[i]) / N_vec[i];
    p2_vec[i] = (X_SI_vec[i] + X_II_vec[i] + X_TI_vec[i] + X_RI_vec[i]) / N_vec[i];
}

// Calculate force of infection
for(i = 0; i < nA; i++) {
    lambda1_vec[i] = 0.0;
    lambda2_vec[i] = 0.0;
    for(j = 0; j < nA; j++) {
        lambda1_vec[i] += CM_mat[i][j] * p1_vec[j];
        lambda2_vec[i] += CM_mat[i][j] * p2_vec[j];
    }
    lambda1_vec[i] *= b1 * exp(eta_ah1 * ah + eta_temp1 * temp);
    lambda2_vec[i] *= tau_vec[i] * b2 * exp(eta_ah2 * ah + eta_temp2 * temp);
}

delta2 = d2 * delta1; // 1 / duration of refractory period (virus2 -> virus1)

// Check values of parameters, for debugging
if(debug) {
	/*for(i = 0; i < nA; i++) {
		for(j = 0; j < nA; j++) {
			Rprintf("%.2f ", CM_mat[i][j]);
		}
		Rprintf("\n");
	}
    Rprintf("\n\n");*/
    for(i = 0; i< nA; i++) {
        //Rprintf("tau[%d] = %.2f\n", i + 1, tau_vec[i]);
        //Rprintf("N[%d] = %.0f\n", i + 1, N_vec[i]);
        Rprintf("lambda[%d] = %.5f\n", i + 1, lambda2_vec[i]);
        Rprintf("p[%d] = %.5f\n", i + 1, p2_vec[i]);
    }
    Rprintf("\n\n");
}

// ODEs
for(i = 0; i < nA; i++) {
    DX_SS_vec[i] = -(lambda1_vec[i] + lambda2_vec[i]) * X_SS_vec[i]; 
    DX_IS_vec[i] = lambda1_vec[i] * X_SS_vec[i] - (gamma1 + theta_lambda1 * lambda2_vec[i]) * X_IS_vec[i]; 
    DX_TS_vec[i] = gamma1 * X_IS_vec[i] - (delta1 + theta_lambda1 * lambda2_vec[i]) * X_TS_vec[i];
    DX_RS_vec[i] = delta1 * X_TS_vec[i] - lambda2_vec[i] * X_RS_vec[i];
    
    DX_SI_vec[i] = lambda2_vec[i] * X_SS_vec[i] - (theta_lambda2 * lambda1_vec[i] + gamma2) * X_SI_vec[i];
    DX_II_vec[i] = theta_lambda1 * lambda2_vec[i] * X_IS_vec[i] + theta_lambda2 * lambda1_vec[i] * X_SI_vec[i] - (gamma1 + gamma2) * X_II_vec[i]; 
    DX_TI_vec[i] = theta_lambda1 * lambda2_vec[i] * X_TS_vec[i] + gamma1 * X_II_vec[i] - (delta1 + gamma2) * X_TI_vec[i];
    DX_RI_vec[i] = lambda2_vec[i] * X_RS_vec[i] + delta1 * X_TI_vec[i] - gamma2 * X_RI_vec[i];
    
    DX_ST_vec[i] = gamma2 * X_SI_vec[i] - (theta_lambda2 * lambda1_vec[i] + delta2) * X_ST_vec[i]; 
    DX_IT_vec[i] = gamma2 * X_II_vec[i] + theta_lambda2 * lambda1_vec[i] * X_ST_vec[i] - (gamma1 + delta2) * X_IT_vec[i]; 
    DX_TT_vec[i] = gamma2 * X_TI_vec[i] + gamma1 * X_IT_vec[i] - (delta1 + delta2) * X_TT_vec[i];
    DX_RT_vec[i] = gamma2 * X_RI_vec[i] + delta1 * X_TT_vec[i] - delta2 * X_RT_vec[i];
    
    DX_SR_vec[i] = delta2 * X_ST_vec[i] - lambda1_vec[i] * X_SR_vec[i]; 
    DX_IR_vec[i] = delta2 * X_IT_vec[i] + lambda1_vec[i] * X_SR_vec[i] - gamma1 * X_IR_vec[i]; 
    DX_TR_vec[i] = delta2 * X_TT_vec[i] + gamma1 * X_IR_vec[i] - delta1 * X_TR_vec[i]; 
    DX_RR_vec[i] = delta2 * X_RT_vec[i] + delta1 * X_TR_vec[i];

    DH1tot_vec[i] = gamma1 * p1_vec[i]; // Incidence rate of virus 1 infections (total)
    DH2tot_vec[i] = gamma2 * p2_vec[i]; // Incidence rate of virus 2 infections (total)
    DH1_vec[i] = gamma1 * (X_IS_vec[i] + theta_rho2 * (X_II_vec[i] + X_IT_vec[i]) + X_IR_vec[i]) / N_vec[i]; // Incidence rate of reported virus 1 infections
    DH2_vec[i] = gamma2 * (X_SI_vec[i] + theta_rho1 * (X_II_vec[i] + X_TI_vec[i]) + X_RI_vec[i]) / N_vec[i]; // Incidence rate of reported virus 2 infections 
}
//end_skel
