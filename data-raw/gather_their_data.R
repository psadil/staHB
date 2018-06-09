library(magrittrm )

delay <- readr::read_table2("https://raw.githubusercontent.com/michaelkalish/STA/master/STACMR-R/Data%20files/delay.dat",
                            col_names = c("participant", "group", "dv", "condition_1", "condition_2", "condition_3", "condition_4")) %>%
  as.data.frame()
dfie <- readr::read_table2("https://raw.githubusercontent.com/michaelkalish/STA/master/STACMR-R/Data%20files/dfie.dat",
                           col_names = c("participant", "condition", "dv", "successes", "failures"))  %>%
  as.data.frame()

# imports file called data
load(file.path(devtools::package_file(), "data-raw", "Supplementary", "WMdata.RData"))
wmdata <- tibble::as_tibble(data)

usethis::use_data(delay, dfie, wmdata, overwrite = TRUE)
