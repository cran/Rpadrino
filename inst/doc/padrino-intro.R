## ----message = FALSE----------------------------------------------------------

library(Rpadrino)

pdb <- pdb_download(save = FALSE)



## -----------------------------------------------------------------------------

pdb


## -----------------------------------------------------------------------------
sub_ind <- setdiff(pdb$Metadata$ipm_id, 
                   c("aaaa15", "aaaa16", "aaaa21", "aaaa54", "aaaa59"))

det_pdb <- pdb_subset(pdb, ipm_ids = sub_ind)


## -----------------------------------------------------------------------------

aster_ind <- det_pdb$Metadata$ipm_id[det_pdb$Metadata$tax_family == "Asteraceae"]

aster_pdb <- pdb_subset(pdb, ipm_ids = aster_ind)


## -----------------------------------------------------------------------------

proto_list <- pdb_make_proto_ipm(aster_pdb)


## -----------------------------------------------------------------------------

proto_list 


## -----------------------------------------------------------------------------

cirsiums <- pdb_make_ipm(proto_list)

lambdas  <- lambda(cirsiums)

cirsiums

## -----------------------------------------------------------------------------

# This creates a table of the number of state variables per ipm_id. Since
# simple IPMs, by defintion, have only 1 state variable, we can use the names
# of the vector that it returns to choose only those. 

n_state_vars <- table(pdb$StateVariables$ipm_id) 
simple_mod_ind <- names(n_state_vars)[n_state_vars == 1]


# Next, we want to capture any models that have stochastic dynamics, multiple
# values per parameter, or density dependence. 

rm_mod_ind <- unique(c(pdb$EnvironmentalVariables$ipm_id, # Stochastic models
                       pdb$ParSetIndices$ipm_id,          # Multiple parameter values
                       pdb$Metadata$ipm_id[pdb$Metadata$has_dd], # Density dependent
                       pdb$Metadata$ipm_id[pdb$Metadata$continent != "n_america"]))

# Remove these IDs from our simple_mod_ind vector

simple_mod_ind <- simple_mod_ind[!simple_mod_ind %in% rm_mod_ind]

# We're also going to remove monocarpic perennials, as their survival/growth
# kernels are slightly trickier to work with (note that there is an example
# of working with these in the manuscript/si/case_study_1.pdf file in PADRINO
# GitHub repository).

monocarps <- pdb$Metadata$ipm_id[pdb$Metadata$organism_type == "biennial"]

simple_mod_ind <- simple_mod_ind[!simple_mod_ind %in% monocarps]

# Finally, we'll subset the database and build the IPM objects!

simple_pdb <- pdb_subset(pdb, simple_mod_ind) %>%
  pdb_make_proto_ipm() %>%
  pdb_make_ipm()



## -----------------------------------------------------------------------------

l_a <- function(ipm, a) {
    
  P <- ipm$sub_kernels$P
  
  # %^% is a function from ipmr that raises matrices to a power, rather than
  # a pointwise power that ^ does.
  
  P_a <- P %^% a
  
  colSums(P_a)
}

f_a <- function(ipm, a) {
    
  P <- ipm$sub_kernels$P
  F <- ipm$sub_kernels$F
  
  l_age <- l_a(ipm, a)
  
  P_a <- P %^% a
  
  colSums(F %*% P_a) / l_age

}


## ---- fig.height = 8, fig.width = 8-------------------------------------------

l_as <- lapply(simple_pdb,
               function(x, a) l_a(x, a), 
               a = 5)

f_as <- lapply(simple_pdb,
               function(x, a) f_a(x, a), 
               a = 5)

# This only plots the figures for the first two species.
# Remove the [1:2] to see all of them.
# Uncomment the par(mfrow = c(...)) line to get an arrangement you like

# par(mfrow = c(2, 2))
for(i in seq_along(l_as)[1:2]) {
  
  nm <- pdb$Metadata$species_accepted[pdb$Metadata$ipm_id == names(l_as)[i]]
  
  plot(l_as[[i]], type = "l",
       # ylim = c(0, 1),
       main = paste0(nm,": Probability of survival to age 5"),
       xlab = expression(paste("Initial size z"[0])),
       ylab = "Pr(s)")
  
  plot(f_as[[i]], type = "l",
       # ylim = c(0, 1),
       main = paste0(nm,": Expected Fecundity at age 5 (given survival)"),
       xlab = expression(paste("Initial size z"[0])),
       ylab = "E[f]")
  
}


## ---- fig.height = 8, fig.width = 8-------------------------------------------

keep_ind <- pdb_subset(pdb, simple_mod_ind) %>%
  .$Metadata %>%
  .[!duplicated(.$species_accepted), "ipm_id"]

use_ipms <- simple_pdb[keep_ind]


n_yrs <- 10

init_pops <- right_ev(use_ipms)

# As above, remove the [1:2] to see all plots and use par() to control their
# arrangement
# par(mfrow = c(3, 2))

for(i in seq_along(init_pops)[1:2]) {
  
  f_age <- l_age <- numeric(n_yrs) 
  
  P_a <- diag(nrow(use_ipms[[i]]$sub_kernels[[1]]))
  
  
  for(j in seq_len(n_yrs)) {
    
    P_now <- use_ipms[[i]]$sub_kernels$P
    F_now <- use_ipms[[i]]$sub_kernels$F

    l_age[j] <- sum(colSums(P_a) * init_pops[[i]][[1]]) 
    f_age[j] <- sum(colSums(F_now %*% P_a) * init_pops[[i]][[1]])
    
    P_a <- P_now %*% P_a
  }
  
  f_age <- f_age / l_age
  
  nm <- pdb$Metadata$species_accepted[pdb$Metadata$ipm_id == names(init_pops)[i]]
  
  plot(l_age, type = "l",
       ylim = c(0, 1),
       main = paste0(nm, ": Probability of survival"),
       xlab = "Age",
       ylab = "Pr(s)")
  
  plot(f_age, type = "l",
       # ylim = c(0, 1),
       main = paste0(nm, ": Average Fecundity"),
       xlab = "Age",
       ylab = "E[f]")
  
}


## -----------------------------------------------------------------------------

make_N <- function(ipm) {
  
  P <- ipm$sub_kernel$P
  I <- diag(nrow(P))
  N <- solve(I - P)
  
  return(N)
}
eta_bar_z0 <- function(ipm) {
  
  N <- make_N(ipm)
  return(colSums(N))
  
}

mean_lifespan <- lapply(use_ipms, eta_bar_z0) 


## -----------------------------------------------------------------------------

sigma_eta_z0 <- function(ipm) {
  
  N <- make_N(ipm)
  
  out <- colSums(2 * N %^% 2 - N) - (colSums(N)) ^ 2
  
  return(out)
  
}

var_lifespan <- lapply(use_ipms, sigma_eta_z0)

sd_lifespan  <- lapply(var_lifespan, sqrt)


## -----------------------------------------------------------------------------

vapply(mean_lifespan, range, numeric(2L))

vapply(var_lifespan, range, numeric(2L))


## ---- fig.height = 10, fig.width = 10-----------------------------------------

mean_lifespan <- mean_lifespan[!names(mean_lifespan) %in% "aaa341"]
sd_lifespan  <- sd_lifespan[!names(sd_lifespan) %in% "aaa341"]

library(ggplot2)

all_data <- data.frame(
  id      = NA,
  species = NA,
  mean_ls = NA,
  upper   = NA,
  lower   = NA,
  z_0     = NA
)

for(i in seq_along(mean_lifespan)) {
  
  temp <- data.frame(
    id = names(mean_lifespan)[i],
    species = pdb$Metadata$species_accepted[pdb$Metadata$ipm_id == names(mean_lifespan)[i]],
    mean_ls = mean_lifespan[[i]],
    upper   = mean_lifespan[[i]] + 1.96 * sd_lifespan[[i]],
    lower   = mean_lifespan[[i]] - 1.96 * sd_lifespan[[i]],
    z_0     = seq(1, length(mean_lifespan[[i]]), 1)
  )
  
  all_data <- rbind(all_data, temp)
  
}

# Remove the NA dummy row, and restrict the lower CI to >= 0 (can't have negative
# lifespan)

all_data <- all_data[-1, ]
all_data$lower <- ifelse(all_data$lower < 0, 0, all_data$lower)


# Now, ggplot using facet wrap and geom_ribbon to get the confidence interval

ggplot(all_data, aes(x = z_0, y = mean_ls)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower,
                  ymax = upper),
              fill = "grey50",
              alpha = 0.5) +
  facet_wrap( ~ species,
              scales = "free") +
  theme_bw()



## -----------------------------------------------------------------------------

lonicera_proto <- pdb_make_proto_ipm(pdb, "aaa341")

# Inspect the vital rate expressions
vital_rate_exprs(lonicera_proto)


## ----eval = FALSE-------------------------------------------------------------
#  
#  vital_rate_exprs(proto_ipms) <- pdb_new_fun_form(
#    list(
#      <ipm_id_1> = list(
#        <vital_rate_name_1> = <expression_1>
#        <vital_rate_name_2> = <expression_2>
#      ),
#      <ipm_id_2> = list(
#        <vital_rate_name_3> = <expression_3>
#      )
#    )
#  )
#  

## -----------------------------------------------------------------------------

vital_rate_exprs(lonicera_proto) <-pdb_new_fun_form(
  list(
    aaa341 = list(
      s = pmin(0.98, 1 / (1 + exp(-(si + ss1 * size_1 + ss2 * size_1 ^ 2))))
    )
  )
)

vital_rate_exprs(lonicera_proto)

## ----fig.height = 10, fig.width = 10------------------------------------------

# Rebuild the IPM, then use our functions for mean, variance, and SD on it.
lonicera_ipm <- pdb_make_ipm(lonicera_proto)

lonicera_mu_ls  <- eta_bar_z0(lonicera_ipm$aaa341)
lonicera_var_ls <- sigma_eta_z0(lonicera_ipm$aaa341)
lonicera_sd_ls  <- sqrt(lonicera_var_ls)

temp <- data.frame(
    id = "aaa341",
    species = "Lonicera_maackii",
    mean_ls = lonicera_mu_ls,
    upper   = lonicera_mu_ls + 1.96 * lonicera_sd_ls,
    lower   = lonicera_mu_ls - 1.96 * lonicera_sd_ls,
    z_0     = seq(1, length(lonicera_mu_ls), 1)
  )


all_data <- rbind(all_data, temp)

# Again, restrict the lower CI to have minimum of 0
all_data$lower <- ifelse(all_data$lower < 0, 0, all_data$lower)


# Rebuild our plot!
ggplot(all_data, aes(x = z_0, y = mean_ls)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower,
                  ymax = upper),
              fill = "grey50",
              alpha = 0.5) +
  facet_wrap( ~ species,
              scales = "free") +
  theme_bw()

