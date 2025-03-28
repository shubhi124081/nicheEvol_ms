multiprobit2_eff_off_nicheEvol <- "
functions
{
  int sum2d(int[,] a) {
  int s = 0;
  for ( i in 1:size(a))
  s += sum(a[i]);
  return s;
  }
}

data
{
  int<lower=1> N; // Number of observations
  int<lower=1> J; // Number of dependent variables
  int<lower=0> K; // Number of independent variables
  int<lower=0> E; // Number of effort variables - these must be the last
  matrix[N,K] x; // independent variables
  int<lower=0, upper = 1> y[N, J]; // dependent variables
  matrix[J,J] L_dist; // beta distance matrix
  int<lower=1> U; // upper limit for truncated inverse gamma
  real<lower=0> a; // inverse gamma hyper parameter
  real<lower=0> b; // inverse gamma hyper parameter
  matrix[J,J] I; // identity matrix for effort variables
  real rho_mean; // rho mean for prior
  real alpha_mean; //rho sd for prior
  real<lower=0> rho_sd; // rho sd for prior
  real<lower=0> alpha_sd; // alpha sd for prior
 // real<lower=0> alpha; // keeping alpha fixed for the moment
}

transformed data
{
  int<lower=0> N_pos;
  int<lower=1, upper=N> n_pos[sum2d(y)];
  int<lower=1, upper=J> d_pos[size(n_pos)];
  int<lower=0> N_neg;
  int<lower=1, upper=N> n_neg[(N * J) - size(n_pos)];
  int<lower=1, upper=J> d_neg[size(n_neg)];

  N_pos = size(n_pos);
  N_neg = size(n_neg);
  {
  int i;
  int j;
  i = 1;
  j = 1;
  for (n in 1:N) {
  for (d in 1:J) {
  if (y[n, d] == 1){
  n_pos[i] = n;
  d_pos[i] = d;
  i += 1;
  } else {
  n_neg[j] = n;
  d_neg[j] = d;
  j += 1;
  }
  }
  }
  }
}

parameters
{
  matrix[K, J] beta;
  vector<lower=0>[N_pos] z_pos;
  vector<upper=0>[N_neg] z_neg;
  vector<lower=0>[J] y_devs;
  vector<lower=0>[J] st_devs;
  real<lower=0> rho;
  real<lower=0> alpha;
  vector[J] beta_mu1; // beta centers for Effort
  vector[J] beta_mu; // beta centers for
}

transformed parameters
{
  vector[J] z[N];
  vector[J] zpred[N];
  for (n in 1:N_pos)
  z[n_pos[n], d_pos[n]] = z_pos[n];
  for (n in 1:N_neg)
  z[n_neg[n], d_neg[n]] = z_neg[n];

  for (n in 1:N_pos)
  zpred[n_pos[n], d_pos[n]] = z_pos[n];
  for (n in 1:N_neg)
  zpred[n_neg[n], d_neg[n]] = z_neg[n];
}

model
{
  matrix[J,J] L_K;
  matrix[J,J] cov;
  matrix[N,J] xbeta = x * beta;

  for (i in 1:(J-1)) {
  L_K[i,i]  = 1 + 0.1; // for numerical stability
  for(j in (i + 1):J){
  L_K[i,j] = alpha^2* exp(-(L_dist[i, j])/rho);
  L_K[j,i] = L_K[i,j];
  }
  }

  for( n in 1:J)
  L_K[n, n] = alpha^2 + 0.0001;

  cov = cholesky_decompose(L_K);

  rho ~ inv_gamma(a, b);
  alpha ~ normal(0,1);
  // rho ~ lognormal(rho_mean, rho_sd);
  // alpha ~ lognormal(alpha_mean, alpha_sd);

  for (k in 1:(K-E))
  beta[k] ~ multi_normal_cholesky(beta_mu, diag_pre_multiply(st_devs, cov));

  for (g in (K-E):K)
  beta[g] ~ multi_normal_cholesky(beta_mu1, diag_pre_multiply(st_devs, I));

  beta_mu ~ normal(0, 1);
  // offset should be 1
  beta_mu1 ~ normal(1, 0.000001);
  st_devs ~ cauchy(0, 2.5);
  y_devs ~ cauchy(0, 2.5);

  for (n in 1:J)
  z[n] ~ normal(xbeta[n], y_devs);

  // prediction model
  //for (n in 1:J)
  //zpred[n] ~ normal(xbeta[n], y_devs);

}
"

GPtrun2_mod_newprior <- "
data
{
  int<lower=1> N; // Number of observations
  int<lower=1> J; // Number of dependent variables i.e. number of regression lines
  int<lower=0> K; // Number of independent variables
  matrix[N,K] x;  // Matrix containing data on independent variable
  matrix[N,J] y;  // Matrix containing data on dependent variable
  //matrix[N,K] x_tilde; //new x data
  //int<lower=1> N_tilde;  // number of new xs
  matrix[J,J] L_dist; // distance matrix
  vector[J] beta_mu;  //beta centers
  matrix[J,J] I; // for the correlation matrix for Y
  // vector[J] st_devs;    // for variance of betas
  int<lower=1> U; // upper limit for truncated inverse gamma
  real<lower=0> a; // inverse gamma hyper parameter
  real<lower=0> b; // inverse gamma hyper parameter
  real rho_mean; // rho mean for prior
  real alpha_mean; //rho sd for prior
  real<lower=0> rho_sd; // rho sd for prior
  real<lower=0> alpha_sd; // alpha sd for prior
}

parameters
{
  matrix[K,J] beta;
  vector<lower=0>[J] st_devs;
  vector<lower = 0>[J] y_devs;
  real<lower = 0> rho;
  real<lower = 0> alpha;
  //real<lower = 0> sigma;

}

model
{
  matrix[N,J] xbeta = x * beta;
  matrix[J, J] L_K;
  matrix[J, J] cov;
  matrix[J, J] L_cov;

  for (i in 1:(J-1)) {

    L_K[i,i]  = alpha^2 + 0.0001; // for numerical stability
    for(j in (i + 1):J){
      L_K[i,j] = alpha^2 * exp(-0.5 * (L_dist[i, j])/rho);
      L_K[j,i] = L_K[i,j];

    }
  }

  for( n in 1:J)
    L_K[n, n] = alpha^2 + 0.0001;

  cov = L_K + diag_pre_multiply(st_devs, I);
  L_cov = cholesky_decompose(cov);

  rho ~ lognormal(rho_mean, rho_sd);
  // rho ~ uniform(0, 50);
  alpha ~ lognormal(alpha_mean, alpha_sd);
  st_devs ~ cauchy(0, 2.5);


  for(k in 1:K)
    beta[K] ~ multi_normal_cholesky(rep_vector(0, J) , L_cov);


  y_devs ~ cauchy(0, 2.5);

  for (i in 1:N)
    y[i] ~ multi_normal_cholesky(xbeta[i],diag_pre_multiply(y_devs, I));
}

generated quantities {
vector[N] log_lik;
for (n in 1:N) log_lik[n] = multi_normal_lpdf(y[n] | x[n] * beta, diag_pre_multiply(y_devs, I));
}
"
