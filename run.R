#!/usr/local/bin/Rscript

task <- dyncli::main()

library(dplyr, warn.conflicts = FALSE)
library(purrr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(tibble, warn.conflicts = FALSE)
library(dynwrap, warn.conflicts = FALSE)

library(mfa, warn.conflicts = FALSE)

#   ____________________________________________________________________________
#   Load data                                                               ####

expression <- as.matrix(task$expression)
end_n <- task$priors$end_n
parameters <- task$parameters

# make sure we have at least one end state
if (end_n == 0) {
  end_n <- 1
}

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# perform MFA
m <- mfa::mfa(
  y = expression,
  b = end_n,
  iter = parameters$iter,
  thin = parameters$thin,
  zero_inflation = parameters$zero_inflation,
  pc_initialise = parameters$pc_initialise,
  prop_collapse = parameters$prop_collapse,
  scale_input = parameters$scale_input
)

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

# obtain results
ms <- summary(m) %>%
  mutate(cell_id = rownames(expression)) %>%
  select(cell_id, everything()) %>%
  group_by(cell_id) %>%
  mutate(
    branch = paste0("M", branch),
    branch_certainty = branch_certainty / sum(branch_certainty)
  ) %>%
  ungroup()

end_state_probabilities <- ms %>%
  select(cell_id, branch, branch_certainty) %>%
  spread(branch, branch_certainty, fill = 0)

pseudotime <-
  ms %>%
  group_by(cell_id) %>%
  summarise(pseudotime = sum(branch_certainty * pseudotime)) %>%
  deframe()

#   ____________________________________________________________________________
#   Save output                                                             ####

output <- 
  dynwrap::wrap_data(cell_ids = names(pseudotime)) %>%
  dynwrap::add_end_state_probabilities(
    end_state_probabilities = end_state_probabilities,
    pseudotime = pseudotime
  ) %>%
  dynwrap::add_timings(checkpoints)

dyncli::write_output(output, task$output)
