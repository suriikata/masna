data {
  int<lower=1> n; // total number of data points
  int<lower=1> k; // number of predictors
  matrix[n, k] X; // independent variables
  vector[n] y;    // dependent variable
}

parameters {
  real a;
  vector[k] b;         // slope      
  real<lower=0> sigma; // stdev
}

model {
  //priors
  b ~ normal(0, 10);
  sigma ~ cauchy(0, 2.5);

  // model
  y ~ normal(a + X * b, sigma);
}

