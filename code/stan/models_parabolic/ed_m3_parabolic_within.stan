// #include /pre/license.stan
data {
  int<lower=1> S; // number of unique subjects
  int<lower=1> N; // number of sessions
  int<lower=1> T; // max trials per session
  array[N] int<lower=1, upper=S> subj;      // subject id for each session
  array[N] int<lower=0, upper=1> phase;     // 0=periovulatory, 1=perimenstrual
  array[N] int<lower=1, upper=T> Tsubj;     // trials per session
  array[N, T] real<lower=0> effort_a;
  array[N, T] real<lower=0> amount_a;
  array[N, T] real<lower=0> effort_b;
  array[N, T] real<lower=0> amount_b;
  array[N, T] int choice;
}
parameters {
  // Group-level (on probit scale)
  vector[3] mu_pr;
  // Subject-to-subject SDs (random intercept SDs) on probit scale
  vector<lower=0>[3] sigma_subj;
  // Subject-level random intercepts (CENTERED: drawn from normal distribution)
  vector[S] kE_subj_pr;
  vector[S] kR_subj_pr;
  vector[S] a_subj_pr;
  // Fixed phase effects (on probit scale)
  vector[3] beta_phase;
}
transformed parameters {
  // Session-level parameters after adding subject intercept + phase effect
  vector<lower=0, upper=1>[N]  kE;
  vector<lower=0, upper=6>[N]  kR;
  vector<lower=-10, upper=10>[N] a;
  for (n in 1:N) {
    int s = subj[n];
    // probit-scale linear predictors
    real eta_kE = mu_pr[1] + sigma_subj[1] * kE_subj_pr[s] + beta_phase[1] * phase[n];
    real eta_kR = mu_pr[2] + sigma_subj[2] * kR_subj_pr[s] + beta_phase[2] * phase[n];
    real eta_a  = mu_pr[3] + sigma_subj[3] * a_subj_pr[s]  + beta_phase[3] * phase[n];
    // transform to your original constrained scales
    kE[n] = Phi_approx(eta_kE);
    kR[n] = Phi_approx(eta_kR) * 6;
    a[n]  = (Phi_approx(eta_a) * 20) - 10;
  }
}
model {
  // Priors
  mu_pr       ~ normal(0, 1);
  sigma_subj  ~ normal(0, 0.2);
  // CENTERED: these are drawn from normal(0, 1)
  kE_subj_pr ~ normal(0, 1);
  kR_subj_pr ~ normal(0, 1);
  a_subj_pr  ~ normal(0, 1);
  beta_phase ~ normal(0, 0.5);
  
  for (n in 1:N) {
    real sv_a;
    real sv_b;
    for (t in 1:Tsubj[n]) {
      sv_a = (kR[n] * amount_a[n, t]) - (kE[n] * square(effort_a[n, t]));
      sv_b = (kR[n] * amount_b[n, t]) - (kE[n] * square(effort_b[n, t]));
      choice[n, t] ~ bernoulli_logit(a[n] + (sv_a - sv_b));
    }
  }
}
generated quantities {
  // Group means by phase (0 vs 1), in your original scales
  real<lower=0, upper=1> mu_kE_phase0;
  real<lower=0, upper=1> mu_kE_phase1;
  real<lower=0, upper=6> mu_kR_phase0;
  real<lower=0, upper=6> mu_kR_phase1;
  real<lower=-10, upper=10> mu_a_phase0;
  real<lower=-10, upper=10> mu_a_phase1;
  // Phase differences (phase1 - phase0)
  real delta_kE;
  real delta_kR;
  real delta_a;
  // Log-lik per session
  array[N] real log_lik;
  // PPC per session x trial
  array[N, T] real y_pred;
  // Initialize
  for (n in 1:N)
    for (t in 1:T)
      y_pred[n, t] = -1;
  // Compute phase-specific group means (subject random intercept is mean 0)
  mu_kE_phase0 = Phi_approx(mu_pr[1] + beta_phase[1] * 0);
  mu_kE_phase1 = Phi_approx(mu_pr[1] + beta_phase[1] * 1);
  mu_kR_phase0 = Phi_approx(mu_pr[2] + beta_phase[2] * 0) * 6;
  mu_kR_phase1 = Phi_approx(mu_pr[2] + beta_phase[2] * 1) * 6;
  mu_a_phase0  = (Phi_approx(mu_pr[3] + beta_phase[3] * 0) * 20) - 10;
  mu_a_phase1  = (Phi_approx(mu_pr[3] + beta_phase[3] * 1) * 20) - 10;
  // Phase differences
  delta_kE = mu_kE_phase1 - mu_kE_phase0;
  delta_kR = mu_kR_phase1 - mu_kR_phase0;
  delta_a  = mu_a_phase1  - mu_a_phase0;
  for (n in 1:N) {
    real sv_a;
    real sv_b;
    log_lik[n] = 0;
    for (t in 1:Tsubj[n]) {
      sv_a = (kR[n] * amount_a[n, t]) - (kE[n] * square(effort_a[n, t]));
      sv_b = (kR[n] * amount_b[n, t]) - (kE[n] * square(effort_b[n, t]));
      log_lik[n] += bernoulli_logit_lpmf(choice[n, t] | a[n] + (sv_a - sv_b));
      y_pred[n, t] = bernoulli_rng(inv_logit(a[n] + (sv_a - sv_b)));
    }
  }
}

