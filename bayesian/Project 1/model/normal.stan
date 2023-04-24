data {
  int<lower=1> n;               // number of data points
  vector<lower=0>[n] y;         // data points
  vector<lower=1,upper=3>[n] g; // group identifier
}

parameters {
  real<lower=50,upper=300> mu1;  // mean for the first group
  real<lower=50,upper=300> mu2;  // mean for the second group
  real<lower=50,upper=300> mu3;  // mean for the third group
  real<lower=0,upper=100> sigma; // stdev
}

model {
  // model
  for (i in 1:n) {
    if (g[i] == 1)
      y[i] ~ normal(mu1, sigma);
    else if (g[i] == 2)
      y[i] ~ normal(mu2, sigma);
    else
      y[i] ~ normal(mu3, sigma);
  }
}

generated quantities {
  // calculate difference
  real diff1;
  real diff2;
  real diff3;
  diff1 = mu2 - mu1;
  diff2 = mu3 - mu2;
  diff3 = mu3 - mu1;
}


