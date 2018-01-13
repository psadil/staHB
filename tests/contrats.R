

library(brms)
library(rstanarm)

d2 <- d %>%
  mutate(condition = factor(condition, ordered=FALSE))

f <- brms::brmsformula(evidence_x ~ mo(condition), family = gaussian())

fit_arm <- stan_glm(evidence_x ~ condition, data=d2, QR = TRUE)
