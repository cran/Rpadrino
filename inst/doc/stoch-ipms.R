## ---- message = FALSE, warning = FALSE----------------------------------------
library(Rpadrino)

pdb <- pdb_download(FALSE)

unique(pdb$EnvironmentalVariables$ipm_id)

## -----------------------------------------------------------------------------
stoch_ids <- c("aaaa15", "aaaa16")


## -----------------------------------------------------------------------------
# Step 1 - create proto_ipms, and then figure out which parameters in the
# vital rate expressions we need to work with.

protos <- pdb_make_proto_ipm(pdb, stoch_ids)

protos[[1]]



## ----eval = FALSE-------------------------------------------------------------
#  pdb$EnvironmentalVariables[pdb$EnvironmentalVariables$ipm_id %in% stoch_ids, ]

## ----echo = FALSE-------------------------------------------------------------

library(knitr)
kable(pdb$EnvironmentalVariables[pdb$EnvironmentalVariables$ipm_id %in% stoch_ids, ],
             caption = "EnvironmentalVariables Table, subsetted using 'stoch_ids'")


## -----------------------------------------------------------------------------

use_j <- sample(1:5, size = 50, replace = TRUE)
use_a <- runif(50, 5, 8)

# The names of this data.frame need to be the same as the names of the variables
# we're substituting in the IPM. Otherwise, the function we write in the next
# block won't work!

use_env_data <- data.frame(j = use_j, A_max = use_a)


## -----------------------------------------------------------------------------
sample_env <- function(env_data, index) {

  as.list(env_data[index, ]) %>%
    setNames(names(env_data))

}


## -----------------------------------------------------------------------------

# For this demonstration, we'll create copy in the 'protos' object and overwrite
# that.

protos$aaaa15_2 <- protos$aaaa15

protos$aaaa15_2$env_state <- lapply(
  protos$aaaa15_2$env_state, 
  function(x) return(NA)
)

# In practice, you may just want to to overwrite the proto_ipm object in place.
# We can still recover the object by re-building the proto_ipm list from the
# intact version of the 'pdb' object.

# protos$aaaa15$env_state <- lapply(
#   protos$aaaa15$env_state, 
#   function(x) return(NA)
# )

# If you want to do get rid of both at the same time, use a double-lapply:

# protos <- lapply(protos, 
#                  function(x) {
#                    x$env_state <- lapply(x$env_state, 
#                                          function(y) return(NA))
#                    
# })



## -----------------------------------------------------------------------------

# Replace the environemntal set up with the one we just created using ipmrs'
# define_env_state.

protos$aaaa15_2 <- define_env_state(
  protos$aaaa15_2,
  env_vars = sample_env(use_env_data, index = t),
  data_list = list(use_env_data = use_env_data,
                   sample_env   = sample_env)
)


## ---- message = FALSE---------------------------------------------------------

ipms <- pdb_make_ipm(protos)

lambda(ipms)


## ----echo = FALSE-------------------------------------------------------------
cap_1 <- "Environmental covariate sequence from unaltered PADRINO data\n(first 10 iterations)"
cap_2 <- "Environmental covariate sequence from custom function\n(first 10 iterations)"

kable(head(ipms$aaaa15$env_seq, n = 10), 
             caption = cap_1,
             format = "html",
             table.attr = "style='width:70%;'")


## ----echo = FALSE-------------------------------------------------------------
kable(head(ipms$aaaa15_2$env_seq, n = 10),
             caption = cap_2,
             format = "html",
             table.attr = "style='width:70%;'")


## -----------------------------------------------------------------------------

all.equal(ipms$aaaa15_2$env_seq$j,     use_env_data$j) 
all.equal(ipms$aaaa15_2$env_seq$A_max, use_env_data$A_max)


## ----echo = FALSE-------------------------------------------------------------

knitr::kable(pdb$EnvironmentalVariables[pdb$EnvironmentalVariables$ipm_id %in% stoch_ids, ],
             caption = "EnvironmentalVariables table, subsetted using 'stoch_ids'")

## -----------------------------------------------------------------------------

pdb_2 <- pdb

inds <- which(pdb_2$EnvironmentalVariables$vr_expr_name %in% paste0("a_max_",
                                                                    c("low", 
                                                                      "high")) &
                pdb_2$EnvironmentalVariables$ipm_id == "aaaa15")

pdb_2$EnvironmentalVariables$env_range[inds[1]] <- 3 # Bump min down two steps
pdb_2$EnvironmentalVariables$env_range[inds[2]] <- 9 # Bump max up two steps

# The rest is the usual building workflow: make a proto_ipm, then convert to
# ipm.

new_protos <- pdb_make_proto_ipm(pdb_2, "aaaa15")

new_ipms   <- pdb_make_ipm(new_protos)

lambda(new_ipms)


## ----message = FALSE----------------------------------------------------------

library(rlang)

sub_env_fun <- function(db,          # pdb object
                        parameter,   # parameter created by the distribution
                        new_fun,     # The new distribution
                        new_pars,    # The parameters for the new distribution
                        id = NULL) { # ipm_id to operate on
  
  if(!is.null(id)) {
    
    if(length(id) > 1L) stop("'id' should only have length 1!")
    
    db <- pdb_subset(db, id)
  } 
  
  ev_tab <- db$EnvironmentalVariables %>%
    sub_rng_fun(parameter, new_fun, new_pars) %>%
    sub_rng_pars(parameter, new_pars)
  
  db$EnvironmentalVariables <- ev_tab
  
  return(db)
}


## -----------------------------------------------------------------------------

sub_rng_fun <- function(ev_tab, parameter, new_fun, new_pars) {

  # Get the old function form
  old_fun_ind <- which(ev_tab$vr_expr_name == parameter)

  # Create a new function call, and insert new arguments
  arg_nms     <- syms(names(new_pars))
  new_fun     <- call2(new_fun, !!! syms(names(new_pars)))
  
  # Convert back to a string and insert into the table
  ev_tab$env_function[old_fun_ind] <- expr_text(new_fun)
  
  return(ev_tab)
}


## -----------------------------------------------------------------------------

sub_rng_pars <- function(ev_tab, parameter, new_pars) { 
  
  # Remove the existing parameters associated with the distribution
  # that we've replaced.
  
  if(sum(ev_tab$env_variable == parameter) > 1) {
    rm_ind <- ev_tab$env_variable == parameter & ev_tab$model_type == "Parameter"
    ev_tab <- ev_tab[!rm_ind, ]
  }
  
  
  new_rows <- data.frame(
    ipm_id       = rep(unique(ev_tab$ipm_id), length(new_pars)),
    env_variable = rep(parameter, length(new_pars)),
    vr_expr_name = names(new_pars),
    env_range    = unlist(new_pars),
    env_function = rep("NULL", length(new_pars)),
    model_type   = "Parameter"
  )
  rownames(new_rows) <- NULL
  
  out <- rbind(ev_tab, new_rows)
  return(out)
}

small_pdb <- pdb %>%
  pdb_subset("aaaa15")

new_pdb <- sub_env_fun(db = small_pdb, 
                       id = NULL, 
                       "A_max",
                       "Gamma", 
                       list(a_max_shape = 14, a_max_rate = 3.5))


## ----echo = FALSE-------------------------------------------------------------

kable(new_pdb$EnvironmentalVariables,
      caption = "New 'EnvironmentalVariables' table")

## ---- message = FALSE---------------------------------------------------------

new_ipm <- new_pdb %>% 
  pdb_make_proto_ipm() %>%
  pdb_make_ipm()

lambda(new_ipm)


