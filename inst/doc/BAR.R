## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, warning = FALSE---------------------------------------------------
library(BAR)

## -----------------------------------------------------------------------------
get_oc_BAR(success_prob = c(.1, .5, .8), n_burn_in = 10, tot_num = 150, block_size = 1, reptime = 100, seed = 100)

## -----------------------------------------------------------------------------
get_oc_BAR(success_prob = c(.1, .5, .8), n_burn_in = 10, tot_num = 150, block_size = 1, reptime = 1, output = "raw", seed = 100)

## -----------------------------------------------------------------------------
get_oc_BAR(success_prob = c(.1, .5, .8), n_burn_in = 10, tot_num = 150, block_size = 3, reptime = 100, seed = 100)

## -----------------------------------------------------------------------------
get_oc_BAR(success_prob = c(.1, .5, .8), n_burn_in = 10, tot_num = 150, block_size = 1, reptime = 100, control_arm = "fixed", seed = 100)

## -----------------------------------------------------------------------------
next_allocation_rate_BAR(n = c(20, 25, 30), success_count = c(2, 9, 15), tot_num = 150, power_c = "n/2N")

## -----------------------------------------------------------------------------
next_allocation_rate_BAR(n = c(20, 25, 30), success_count = c(2, 10, 15), tot_num = 150, power_c = .5)

