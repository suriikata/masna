# libraries --------------------------------------------------------------------
library(cmdstanr)
library(ggplot2)
library(tidyverse)
library(bayesplot)
library(posterior)
library(mcmcse)
library(ggdist)

# modelling and data prep ------------------------------------------------------
# model
model <- cmdstan_model("model/multi-normal.stan")

# data
data <- read.csv("data/02_crime.csv")

# let us focus only on police funding, 25 + high school and 16 + no high school => crime_rate
data <- data %>% select(crime_rate, police_funding, X25_plus_high_school, X16_19_no_high_school)

# normalized columns
data$crime <- (data$crime_rate - min(data$crime_rate)) /
           (max(data$crime_rate) - min(data$crime_rate))
data$police <- (data$police_funding - min(data$police_funding)) /
           (max(data$police_funding) - min(data$police_funding))
data$high25 <- (data$X25_plus_high_school - min(data$X25_plus_high_school)) /
          (max(data$X25_plus_high_school) - min(data$X25_plus_high_school))
data$high16 <- (data$X16_19_no_high_school - min(data$X16_19_no_high_school)) /
          (max(data$X16_19_no_high_school) - min(data$X16_19_no_high_school))


# fit without normalization
# prep for Stan
y <- data$crime_rate
X <- data %>% select(police_funding, X25_plus_high_school, X16_19_no_high_school)
stan_data <- list(n = length(y), k = ncol(X), y = y, X = X)

# fit
fit <- model$sample(
  data = stan_data,
  parallel_chains = 4,
  seed = 1
)

# fit with normalization
y <- data$crime
X <- data %>% select(police, high25, high16)
stan_data <- list(n = length(y), k = ncol(X), y = y, X = X)

# fit normalized
fit_normalized <- model$sample(
  data = stan_data,
  parallel_chains = 4,
  seed = 1
)

# diagnostics ------------------------------------------------------------------
# traceplot
mcmc_trace(fit$draws())
mcmc_trace(fit_normalized$draws())

# MCMC converged well > sampling process reached a stable state 

# compare summaries
fit$summary()
fit_normalized$summary()
 
# all means have rhat of 1, both ess exceed the expected value of 400

# extract samples --------------------------------------------------------------
# convert draws to df
df <- as_draws_df(fit$draws())
df_normalized <- as_draws_df(fit_normalized$draws())


# mcse
mcse(df$`b[1]`)   
mcse(df$`b[2]`)    
mcse(df$`b[3]`)   


# ratios of beta for each variables, use those ratios to define the investements
# by normalizing values we would break up the ratios > TAKE CARE HOW YOU HANDLE THE DATA
# crime coef ----------------------------------------------------------------------------------
mcse(df_normalized$`b[1]`)  # for one dollar increase in police funding, .57 unit increase in crime 
mcse(df_normalized$`b[2]`)  # one precentage increase in 25 with high school, -.13 decrease in crime
mcse(df_normalized$`b[3]`)  # 

# first and third are significant at 95 %, whereas the second isnt 

# what the ratio of our investment should be based on looking at beta coefficients?
# sum up to 1 and calculate ratios  

# compare default vs normalized ------------------------------------------------
# calculate ratios
df_ratios_default <- df %>%
  mutate(police = `b[1]` / (`b[1]` + `b[2]` + `b[3]`),
         high25 = `b[2]` / (`b[1]` + `b[2]` + `b[3]`),
         high16 = `b[3]` / (`b[1]` + `b[2]` + `b[3]`)) %>%
  select(police, high25, high16)

df_ratios_normalized <- df_normalized %>%
  mutate(police = `b[1]` / (`b[1]` + `b[2]` + `b[3]`),
         high25 = `b[2]` / (`b[1]` + `b[2]` + `b[3]`),
         high16 = `b[3]` / (`b[1]` + `b[2]` + `b[3]`)) %>%
  select(police, high25, high16)

# column means
mcse(df_ratios_default$police)
mcse(df_ratios_default$high25)
mcse(df_ratios_default$high16)
mcse(df_ratios_normalized$police)
mcse(df_ratios_normalized$high25)
mcse(df_ratios_normalized$high16)


# crime ratios ------------------------------------------------------------------------
# looking at the default ratios, we would invest 65 % police funding, 35 % ind 16-19 with no high school
# normalized ratios CHANGES > 75 % police funding, 25 % reducing the proportion of ind 16-19 with no high school 
# however, the first approach is better, because it keep the original ratios !!

# violence ratios --------------------------------------------------------------------
# municipality should allocate roughly 80 % police funding and 20 % of ind aged 16-19 with no high school 


# add labels
df_ratios_default$Type <- "Default"
df_ratios_normalized$Type <- "Normalized"

# means
df_ratios <- rbind(df_ratios_default, df_ratios_normalized)

# to long format
df_ratios <- df_ratios %>% gather(Variable,
                                  Value,
                                  c(police, high25, high16))
# plot
ggplot(data = df_ratios, aes(x = Value, y = Variable)) +
  stat_eye(fill = "skyblue", alpha = 0.75) +
  facet_grid(Type ~ .) +
  xlim(0, 1)
