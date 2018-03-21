# boiler plate code for creation of directory groups

root <- devtools::package_file()

dir.output <- file.path(root,"output")
dir.stan <- file.path(root,'stan')
dir.report <- file.path(root, 'reports')

dir.stan_samples <- file.path(dir.output, 'stan_samples')
dir.stan_data <- file.path(dir.output, 'stan_data')
dir.figures <- file.path(dir.output, 'figures')
dir.report_output <- file.path(dir.output, 'reports')
dir.report_params <- file.path(dir.output, 'report_params')


if (!dir.exists(dir.stan)){
  dir.create(dir.stan)
}
if (!dir.exists(dir.figures)){
  dir.create(dir.figures)
}
if (!dir.exists(dir.output)){
  dir.create(dir.output)
}
if (!dir.exists(dir.stan_data)){
  dir.create(dir.stan_data)
}
if (!dir.exists(dir.stan_samples)){
  dir.create(dir.stan_samples)
}
if (!dir.exists(dir.report_output)){
  dir.create(dir.report_output)
}
if (!dir.exists(dir.report_params)){
  dir.create(dir.report_params)
}
