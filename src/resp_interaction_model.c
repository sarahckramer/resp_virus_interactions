// Implementation of the SIR model for circulation of respiratory viruses

//start_globs
// Expit transform for parameters constrained in interval [a,b]
static double expitCons(double x, double a, double b) {
  double out = (a + b * exp(x)) / (1.0 + exp(x));
  if(ISNAN(out) | isinf(out)) out = (b + a * exp(-x)) / (1.0 + exp(-x)); // If x=+Inf, must return b
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

//start_toest
double sum_init = 0.0;

T_Ri1 = logitCons(Ri1, 1.0, Ri1_max); 
T_Ri2 = logitCons(Ri2, 1.0, Ri2_max); 
T_gamma1 = log(gamma1); 
T_gamma2 = log(gamma2);
T_delta = logitCons(delta, delta_min, 7 / 1);
//T_delta1 = log(delta1);
//T_delta2 = log(delta2); 
T_theta_lambda1 = logit(theta_lambda1);
T_theta_lambda2 = logit(theta_lambda2);
T_theta_rho1 = logit(theta_rho1);
T_theta_rho2 = logit(theta_rho2);
T_rho1 = logit(rho1); 
T_rho2 = logit(rho2);
T_N = N;

sum_init = I10 + I20 + R10 + R20 + R120;
T_I10 = log(I10 / (1.0 - sum_init));
T_I20 = log(I20 / (1.0 - sum_init));
T_R10 = log(R10 / (1.0 - sum_init));
T_R20 = log(R20 / (1.0 - sum_init));
T_R120 = log(R120 / (1.0 - sum_init));

//end_toest

//start_fromest
double sum_init = 0.0; 

Ri1 = expitCons(T_Ri1, 1.0, Ri1_max);
Ri2 = expitCons(T_Ri2, 1.0, Ri2_max);
gamma1 = exp(T_gamma1);
gamma2 = exp(T_gamma2);
delta = expitCons(T_delta, delta_min, 7 / 1);
//delta1 = exp(T_delta1);
//delta2 = exp(T_delta2);
theta_lambda1 = expit(T_theta_lambda1);
theta_lambda2 = expit(T_theta_lambda2);
theta_rho1 = expit(T_theta_rho1);
theta_rho2 = expit(T_theta_rho2);
rho1 = expit(T_rho1);
rho2 = expit(T_rho2);
N = T_N;

sum_init = exp(T_I10) + exp(T_I20) + exp(T_R10) + exp(T_R20) + exp(T_R120);
I10 = exp(T_I10) / (1.0 + sum_init);
I20 = exp(T_I20) / (1.0 + sum_init);
R10 = exp(T_R10) / (1.0 + sum_init);
R20 = exp(T_R20) / (1.0 + sum_init);
R120 = exp(T_R120) / (1.0 + sum_init);

//end_fromest

//start_dmeas
double fP1, fP2, ll;

double rho1_w = fmin2(1.0, rho1 * H1 / i_ARI); // Probability of detecting virus 1
double rho2_w = fmin2(1.0, rho2 * H2 / i_ARI); // Probability of detecting virus 2

fP1 = dbinom(n_P1, n_T, rho1_w, 1); // First likelihood component, natural scale
fP2 = dbinom(n_P2, n_T, rho2_w, 1); // Second likelihood component, natural scale

//fP1 = ISNA(n_P1) ? 0.0 : dbinom(n_P1, n_T, rho1_w, 1);
//fP2 = ISNA(n_P2) ? 0.0 : dbinom(n_P2, n_T, rho2_w, 1);
//if (isnan(fP1)) {
//  fP1 = dbinom(n_P1, n_T, 0.0, 1);
//}
//if (isnan(fP2)) {
//  fP2 = dbinom(n_P2, n_T, 0.0, 1);
//}

// If rho_w == 1, the resulting observation probability might be 0 (-Inf on log-scale)
// Replace by a big, but finite penalty if that's the case 
ll = fmax2(fP1 + fP2, -1e3);

if(debug) {
  Rprintf("t=%.1f, rho1_w=%.1f, rho2_w=%.1f, n_T=%.1f, fP1=%.1f, fP2=%.1f, sum=%.1f, ll=%.f\n", t, rho1_w, rho2_w, n_T, fP1, fP2, fP1 + fP2, ll);
} 
lik = (give_log) ? ll : exp(ll);

//end_dmeas

//start_rmeas
double rho1_w = fmin2(1.0, rho1 * H1 / i_ARI); // Probability of detecting virus 1
double rho2_w = fmin2(1.0, rho2 * H2 / i_ARI); // Probability of detecting virus 2

n_P1 = rbinom(n_T, rho1_w); // Generate of tests positive to virus 1
n_P2 = rbinom(n_T, rho2_w); // Generate of tests positive to virus 2
//end_rmeas

//start_rinit
X_SS = nearbyint((1.0 - I10 - I20 - R10 - R20 - R120) * N);
X_IS = nearbyint(I10 * N);
X_TS = 0;
X_RS = nearbyint(R10 * N);
X_SI = nearbyint(I20 * N);
X_II = 0;
X_TI = 0;
X_RI = 0;
X_ST = 0;
X_IT = 0;
X_TT = 0;
X_RT = 0;
X_SR = nearbyint(R20 * N);
X_IR = 0;
X_TR = 0;
X_RR = nearbyint(R120 * N);
H1_tot = 0;
H2_tot = 0;
H1 = 0; 
H2 = 0;

if ((X_SS + X_IS + X_RS + X_SI + X_SR + X_RR) != N) {
  X_SS = nearbyint(N - X_IS - X_RS - X_SI - X_SR - X_RR);
}

if ((X_SS + X_IS + X_RS + X_SI + X_SR + X_RR) != N) {
  Rprintf("SS=%f, IS=%f, RS=%f, SI=%f, SR=%f, RR=%f, sum=%f, N=%f\n", X_SS, X_IS, X_RS, X_SI, X_SR, X_RR, X_SS + X_IS + X_RS + X_SI + X_SR + X_RR, N);
}

if(debug) {
  Rprintf("%f, %f, %f, %f, %f, %f, %f\n", Ri1, Ri2, I10, I20, R10, R20, R120);
}

//end_rinit

//start_skel
double p1 = (X_IS + X_II + X_IT + X_IR) / N; // Prevalence of infection with virus 1
double p2 = (X_SI + X_II + X_TI + X_RI) / N; // Prevalence of infection with virus 2

double beta1 = Ri1 / (1.0 - (R10 + R120)) * gamma1; // Initial reproduction no of virus 1 (R10+R120: initial prop immune to v1)
double beta2 = Ri2 / (1.0 - (R20 + R120)) * gamma2; // Initial reproduction no of virus 2 (R20+R120: initial prop immune to v2)
double lambda1 = beta1 * p1; // Force of infection with virus 1
double lambda2 = beta2 * p2; // Force of infection with virus 2

// ODEs
DX_SS = -(lambda1 + lambda2) * X_SS; 
DX_IS = lambda1 * X_SS - (gamma1 + theta_lambda1 * lambda2) * X_IS; 
DX_TS = gamma1 * X_IS - (delta + theta_lambda1 * lambda2) * X_TS;
DX_RS = delta * X_TS - lambda2 * X_RS; 
DX_SI = lambda2 * X_SS - (theta_lambda2 * lambda1 + gamma2) * X_SI;
DX_II = theta_lambda1 * lambda2 * X_IS + theta_lambda2 * lambda1 * X_SI - (gamma1 + gamma2) * X_II; 
DX_TI = theta_lambda1 * lambda2 * X_TS + gamma1 * X_II - (delta + gamma2) * X_TI;
DX_RI = lambda2 * X_RS + delta * X_TI - gamma2 * X_RI; 
DX_ST = gamma2 * X_SI - (theta_lambda2 * lambda1 + delta) * X_ST; 
DX_IT = gamma2 * X_II + theta_lambda2 * lambda1 * X_ST - (gamma1 + delta) * X_IT; 
DX_TT = gamma2 * X_TI + gamma1 * X_IT - (2 * delta) * X_TT;
DX_RT = gamma2 * X_RI + delta * X_TT - delta * X_RT;
DX_SR = delta * X_ST - lambda1 * X_SR; 
DX_IR = delta * X_IT + lambda1 * X_SR - gamma1 * X_IR; 
DX_TR = delta * X_TT + gamma1 * X_IR - delta * X_TR; 
DX_RR = delta * X_RT + delta * X_TR;

//DX_SS = -(lambda1 + lambda2) * X_SS; 
//DX_IS = lambda1 * X_SS - (gamma1 + theta_lambda1 * lambda2) * X_IS; 
//DX_TS = gamma1 * X_IS - (delta1 + theta_lambda1 * lambda2) * X_TS;
//DX_RS = delta1 * X_TS - lambda2 * X_RS; 
//DX_SI = lambda2 * X_SS - (theta_lambda2 * lambda1 + gamma2) * X_SI;
//DX_II = theta_lambda1 * lambda2 * X_IS + theta_lambda2 * lambda1 * X_SI - (gamma1 + gamma2) * X_II; 
//DX_TI = theta_lambda1 * lambda2 * X_TS + gamma1 * X_II - (delta1 + gamma2) * X_TI;
//DX_RI = lambda2 * X_RS + delta1 * X_TI - gamma2 * X_RI; 
//DX_ST = gamma2 * X_SI - (theta_lambda2 * lambda1 + delta2) * X_ST; 
//DX_IT = gamma2 * X_II + theta_lambda2 * lambda1 * X_ST - (gamma1 + delta2) * X_IT; 
//DX_TT = gamma2 * X_TI + gamma1 * X_IT - (delta1 + delta2) * X_TT;
//DX_RT = gamma2 * X_RI + delta1 * X_TT - delta2 * X_RT;
//DX_SR = delta2 * X_ST - lambda1 * X_SR; 
//DX_IR = delta2 * X_IT + lambda1 * X_SR - gamma1 * X_IR; 
//DX_TR = delta2 * X_TT + gamma1 * X_IR - delta1 * X_TR; 
//DX_RR = delta2 * X_RT + delta1 * X_TR; 

DH1_tot = gamma1 * p1; // Incidence rate of virus 1 infections (total)
DH2_tot = gamma2 * p2; // Incidence rate of virus 2 infections (total)
DH1 = gamma1 * (X_IS + theta_rho2 * (X_II + X_IT) + X_IR) / N; // Incidence rate of reported virus 1 infections
DH2 = gamma2 * (X_SI + theta_rho1 * (X_II + X_TI) + X_RI) / N; // Incidence rate of reported virus 2 infections 
//end_skel

//start_rsim
//end_rsim
