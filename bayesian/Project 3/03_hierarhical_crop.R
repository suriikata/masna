# libraries --------------------------------------------------------------------
library(cmdstanr)
library(mcmcse)
library(ggplot2)
library(bayesplot)
library(posterior)
library(tidyverse)
library(HDInterval)
library(cowplot)

# data prep and model compilation ----------------------------------------------
# load data
data <- read.csv("data/01_crop.csv")

# hierarchical normal model ----------------------------------------------------
model_h <- cmdstan_model("model/hierarchical_normal.stan")

# data prep
stan_data <- list(n = nrow(data),
                  m = max(data$fertilizer),
                  y = data$yield,
                  g = data$fertilizer)

# fit
fit_h <- model_h$sample(
  data = stan_data,
  parallel_chains = 4,
  seed = 1
)

# diagnostics
mcmc_trace(fit_h$draws())
fit_h$summary()

df_1 <- as_draws_df(fit_h$draws("mu[1]"))
df_2 <- as_draws_df(fit_h$draws("mu[2]"))
df_3 <- as_draws_df(fit_h$draws("mu[3]"))


data1 <- data.frame(mu = c(df_1$`mu[1]` , df_2$`mu[2]` , df_3$`mu[3]` )
                   , names = c(rep("Fertilizer 1", length(df_1$`mu[1]`)),rep("Fertilizer 2", length(df_2$`mu[2]`)),rep("Fertilizer 3", length(df_3$`mu[3]`))))


# plot means
library(viridis)
my_palette <- viridis_pal(option = "plasma")(5)

ggplot(data1, aes(mu, fill = names)) + geom_density(alpha = 0.8) + 
  scale_fill_manual(values = my_palette,
                    name = "Legend") +
  ylab("Density") + xlab("Mean") +
  ggtitle("Distribution of Means") +
  theme(legend.position = c(0.85, 0.85), legend.justification = c(1, 1)) +
  theme_minimal()


mcse(df_1$`mu[1]` < df_2$`mu[2]`)
mcse(df_1$`mu[1]`< df_3$`mu[3]`)
mcse(df_2$`mu[2]` < df_3$`mu[3]`)

# probability F3 is supreme
mcse(df_2$`mu[2]` < df_3$`mu[3]`)$est * mcse(df_1$`mu[1]` < df_3$`mu[3]`)$est
# probability F2 is supreme
mcse(df_2$`mu[2]` > df_3$`mu[3]`)$est * mcse(df_2$`mu[2]` > df_1$`mu[1]`)$est
# probability F1 is supreme
mcse(df_1$`mu[1]` > df_3$`mu[3]`)$est * mcse(df_1$`mu[1]` > df_2$`mu[2]`)$est


mcse((df_1$`mu[1]` > df_2$`mu[2]`) & (df_1$`mu[1]` > df_3$`mu[3]`))
mcse((df_1$`mu[1]`< df_3$`mu[3]`) & (df_2$`mu[2]` < df_3$`mu[3]`))
mcse((df_2$`mu[2]` > df_3$`mu[3]`) & (df_2$`mu[2]` > df_1$`mu[1]`))


# samples
df_h <- as_draws_df(fit_h$draws(c("sigma", "mu", "mu_mu", "sigma_mu")))
df_h <- df_h %>% select(-.draw, -.chain, -.iteration)

# visual posterior check
# use only n_dist distributions
df_sample_h <- sample_n(df_h, n_dist)

# prep for plotting
df_generated_h <- data.frame(x = numeric(),
                             y = factor(),
                             iteration = numeric(),
                             fertilizer = numeric())

for (i in 1:n_fert) {
  for (j in 1:n_dist) {
    # mu for fertilizer i is in column i+1
    # sigma is always in the first column
    y <- dnorm(x,
               mean = df_sample_h[j, i + 1][[1]],
               sd = df_sample_h[j, ]$sigma)

    df_generated_h <- rbind(df_generated_h,
                            data.frame(x = x, y = y,
                                       iteration = j, fertilizer = i))
  }
}

# plot
ggplot() +
  geom_density(data = data, aes(x = yield),
               fill = "skyblue", alpha = 0.75, color = NA) +
  geom_line(data = df_generated_h,
            aes(x = x, y = y, group = iteration), alpha = 0.1, size = 1) +
  facet_wrap(. ~ fertilizer, ncol = 4) +
  xlim(0, 180) +
  xlab("Weight") +
  ylab("Density")

# compare group level means ----------------------------------------------------
df_group <- data.frame(Mean = numeric(),
                       HDI5 = numeric(),
                       HDI95 = numeric(),
                       Model = character())

# sample
sample_mean <- mean(data$yield)
df_group <- rbind(df_group, data.frame(Mean = sample_mean,
                                       HDI5 = sample_mean,
                                       HDI95 = sample_mean,
                                       Model = "Sample"))

# simple normal model
normal_mean <- mean(df_n$mu)
normal_90_hdi <- hdi(df_n$mu, credMass = 0.9)
df_group <- rbind(df_group, data.frame(Mean = normal_mean,
                                       HDI5 = normal_90_hdi[1],
                                       HDI95 = normal_90_hdi[2],
                                       Model = "Normal"))

# hierarchical model
hierarchical_mean <- mean(df_h$mu_mu)
hierarchical_90_hdi <- hdi(df_h$mu_mu, credMass = 0.9)
df_group <- rbind(df_group, data.frame(Mean = hierarchical_mean,
                                       HDI5 = hierarchical_90_hdi[1],
                                       HDI95 = hierarchical_90_hdi[2],
                                       Model = "Hierarchical"))

# plot
# set model factors so the colors are the same
df_group$Model <- factor(df_group$Model,
                         levels = c("Normal", "Hierarchical", "Sample"))

ggplot(data = df_group,
       aes(x = Model,
           y = Mean,
           ymin = HDI5,
           ymax = HDI95,
           colour = Model)) +
  geom_point() +
  geom_errorbar(width = 0.2) +
  scale_color_brewer(palette = "Set1") +
  ylim(0, 180) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# compare subject level means --------------------------------------------------
df_subject <- data.frame(Mean = numeric(),
                         Q5 = numeric(),
                         Q95 = numeric(),
                         Model = character(),
                         fertilizer = numeric())

# sample means
df_mu_sample <- data %>%
  group_by(fertilizer) %>%
  summarise(mean_weight = mean(yield))
df_subject <- rbind(df_subject, data.frame(Mean = df_mu_sample$mean_weight,
                                           HDI5 = df_mu_sample$mean_weight,
                                           HDI95 = df_mu_sample$mean_weight,
                                           Model = "Sample",
                                           fertilizer = seq(1:n_fert)))

# subject means
df_mu_s <- df_s %>% select(2:(1 + n_fert))
s_means <- colMeans(df_mu_s)
s_90_hdi <- apply(df_mu_s, 2, hdi, credMass = 0.9)
df_subject <- rbind(df_subject, data.frame(Mean = s_means,
                                           HDI5 = s_90_hdi[1, ],
                                           HDI95 = s_90_hdi[2, ],
                                           Model = "Group",
                                           fertilizer = seq(1:n_fert)))

# hierarchical means
df_mu_h <- df_h %>% select(2:(1 + n_fert))
h_means <- colMeans(df_mu_h)
h_hdi90 <- apply(df_mu_h, 2, hdi, credMass = 0.9)
df_subject <- rbind(df_subject, data.frame(Mean = h_means,
                                           HDI5 = h_hdi90[1, ],
                                           HDI95 = h_hdi90[2, ],
                                           Model = "Hierarchical",
                                           fertilizer = seq(1:n_fert)))

# plot
ggplot(data = df_subject,
       aes(x = Model,
           y = Mean,
           ymin = HDI5,
           ymax = HDI95,
           colour = Model)) +
  geom_hline(yintercept = mean(data$yield), color = "grey75") +
  geom_point() +
  geom_errorbar(width = 0.2) +
  scale_color_brewer(palette = "Set1") +
  ylim(175, 180) +
  facet_wrap(. ~ fertilizer, ncol = 4) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
