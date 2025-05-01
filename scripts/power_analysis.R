#Simulate model----

simulate_binomial_glm <- function(n, effect_promoter, effect_degron, alpha = 0.05) {
  # simulated data will look like
  # promoter degron success
  # 0        0      0 or 1
  # 1        0      0 or 1
  # 0        1      0 or 1
  # 1        1      0 or 1
  # * n (so n is the sample size required for each group)

  promoter <- rep(c(0,1,0,1), each = n)
  degron <- rep(c(0,0,1,1), each = n)
  
  #set probabilty for each treatment group based on intended effect size
  prob_00 <- 0.5
  prob_10 <- 0.5 + effect_promoter
  prob_01 <- 0.5 + effect_degron
  prob_11 <- 0.5 + effect_promoter + effect_degron
  
  #using rbinom to generate n success/failures for each possible group 
  success_00 <- rbinom(n, size = 1, prob = prob_00)
  success_10 <- rbinom(n, size = 1, prob = prob_10)
  success_01 <- rbinom(n, size = 1, prob = prob_01)
  success_11 <- rbinom(n, size = 1, prob = prob_11)
  
  successes = c(success_00, success_10, success_01, success_11)
  
  data = data.frame(success = successes, promoter = as.factor(promoter), degron = as.factor(degron))
  
  # run the model. power analysis is based off a simplified model with no interaction term
  
  model <- glm(success ~ promoter + degron, family = binomial(link = "logit"), data = data)
  
  return(
    list(
      promoter = summary(model)$coefficients[2,4] < alpha,
      degron = summary(model)$coefficients[3,4] < alpha
    )
  )
}


#this function is the same as in the example, just with effect sizes for both predictors
monte_carlo_power <- function(n, n_sims = 1000, effect_promoter = 0.1, effect_degron = 0.1, alpha = 0.05) {
  results <- replicate(n_sims, simulate_binomial_glm(n, effect_promoter, effect_degron, alpha), simplify = "array")
  return(
    list(
      promoter = mean(as.numeric(results[1,])),
      degron = mean(as.numeric(results[2,]))
    )
    )
}

#Power curve----

#run in sequence :>

sample_sizes <- seq(10, 300, by = 10)

#create an empty tibble and add a row for each sample size
power_results <- tibble(n = numeric(), promoter = numeric(), degron = numeric())
for (i in sample_sizes){
  x <- monte_carlo_power(i, effect_promoter = 0.1, effect_degron = 0.1)
  
  power_results <- power_results |>
    add_row(n = i, promoter = x$promoter, degron = x$degron)
}

#make longer for ggplot
power_results <- power_results |>
  pivot_longer(
    cols = !n,
    names_to = "variable",
    values_to = "power"
  )

#plot power curve
power_results |>
  ggplot(aes(x = n, y = power, colour = variable)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  scale_x_continuous(breaks = seq(0,300, by = 20))
