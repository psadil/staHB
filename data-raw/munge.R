cutoff <- 1
type <- c("full","mon","nmon")

raw_dir <- file.path(devtools::package_file(), "data-raw")

targets <- readr::read_csv(file.path(raw_dir, 'objectNames.csv'), col_names = FALSE)


d <- readr::read_csv(file.path(raw_dir, 'data.csv')) %>%
  dplyr::mutate(targets = purrr::map(object, ~targets[.x,])) %>%
  dplyr::mutate(firstTarget = purrr::map_chr(targets, ~magrittr::extract2(.x,1))) %>%
  dplyr::mutate(dl = purrr::map2(targets, responses, ~stringdist::stringdist(purrr::as_vector(.x), .y, method="dl"))) %>%
  dplyr::mutate(minDist = purrr::map_dbl(dl, ~min(.x, na.rm=TRUE))) %>%
  dplyr::mutate(y2 = dplyr::if_else(minDist < cutoff, 1L, 0L)) %>%
  dplyr::rename(y1 = afc,
         item = object) %>%
  dplyr::select(-pas, -pas2,-objectPair,-responses,-targets,-dl, -minDist, -firstTarget,-list) %>%
  tidyr::crossing(., type = type) %>%
  dplyr::mutate(expt = dplyr::case_when(type == "mon" ~ 1,
                                         type == "nmon" ~ 2,
                                         type == "full" ~ 3)) %>%
  dplyr::group_by(expt) %>%
  tidyr::nest() %>%
  dplyr::mutate(stan_data = purrr::map(data, ~gen_stan_data(.x)))


usethis::use_data(d, overwrite = TRUE)

expt1_post <- d %>%
  dplyr::mutate(post = purrr::map(stan_data, ~  run_stan(.x, iter=500, warmup=1000, chains=6, model=2)))

usethis::use_data(expt1_post, overwrite = TRUE)


expt1_post_rhos <- d %>%
  dplyr::mutate(post = purrr::map(stan_data, ~  run_stan(.x, iter = 1000, warmup=1000, chains=6, model=3)))

usethis::use_data(expt1_post_rhos, overwrite = TRUE)


post1 <- expt1_post %>%
  dplyr::mutate(log_lik = purrr::map(post, ~loo::extract_log_lik(.x, merge_chains = FALSE)),
                r_eff = purrr::map(log_lik, ~loo::relative_eff(.x, cores=4))) %>%
  dplyr::mutate(waic = purrr::map(log_lik, ~loo::waic(.x))) %>%
  dplyr::mutate(loo = purrr::map2(log_lik, r_eff, ~ loo::loo(.x, r_eff = .y)))


post <- expt1_post_rhos %>%
  dplyr::mutate(log_lik = purrr::map(post, ~loo::extract_log_lik(.x,merge_chains = FALSE)),
                r_eff = purrr::map(log_lik, ~loo::relative_eff(.x, cores=4))) %>%
  dplyr::mutate(waic = purrr::map(log_lik, ~loo::waic(.x))) %>%
  dplyr::mutate(loo = purrr::map2(log_lik, r_eff, ~ loo::loo(.x, r_eff = .y)))

