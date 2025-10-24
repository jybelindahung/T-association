data {
  int<lower=1> N;
  array[N] vector[2] y;
  vector<lower=0>[N] se1;
  vector<lower=0>[N] se2;
  real<lower=0> ig_shape;
  real<lower=0> ig_scale;
  int<lower=0, upper=1> use_uniform;    // 1 = Uniform(-0.995,0.995), 0 = Fisher-z prior
  int<lower=0, upper=1> use_half_t_psi; // 1 = half-t, 0 = inverse-gamma
  real<lower=0> half_t_df;              // degrees of freedom
  real<lower=0> half_t_scale;           // scale for half-t
}

parameters {
  real beta1;
  real beta2;
  real<lower=1e-2, upper=10> psi1; // SD
  real<lower=1e-2, upper=10> psi2;

  real<lower=-0.99, upper=0.99> rho_uniform;
  real rho_z;
}

transformed parameters {
  real psi1_2 = psi1^2;
  real psi2_2 = psi2^2;
  real rho;

  if (use_uniform == 1) {
    rho = rho_uniform;
  } else {
    rho = tanh(rho_z) * 0.9999;
  }
}

model {
  // Regression priors
  beta1 ~ normal(0, 2.5);
  beta2 ~ normal(0, 2.5);

  // Priors for psi1 and psi2
  if (use_half_t_psi == 1) {
    // Half-t prior: psi ~ Student_t(df, 0, scale) truncated at 0
    psi1 ~ student_t(half_t_df, 0, half_t_scale);
    psi2 ~ student_t(half_t_df, 0, half_t_scale);
  } else {
    // Inverse-gamma prior: psi^2 ~ IG(shape, scale), with Jacobian
    target += inv_gamma_lpdf(psi1^2 | ig_shape, ig_scale) + log(2*psi1);
    target += inv_gamma_lpdf(psi2^2 | ig_shape, ig_scale) + log(2*psi2);
  }

  // Priors for correlation
  if (use_uniform == 1) {
    rho_uniform ~ uniform(-0.99, 0.99);
  } else {
    rho_z ~ normal(0, 1);
  }

  // Likelihood
  vector[2] mu;
  mu[1] = beta1;
  mu[2] = beta2;

  for (i in 1:N) {
    real eps = 1e-6;
    matrix[2,2] Sigma;
    Sigma[1,1] = psi1^2 + se1[i]^2 + eps;
    Sigma[2,2] = psi2^2 + se2[i]^2 + eps;
    Sigma[1,2] = rho * sqrt(Sigma[1,1]) * sqrt(Sigma[2,2]);
    Sigma[2,1] = Sigma[1,2];
    y[i] ~ multi_normal(mu, Sigma);
  }
}
