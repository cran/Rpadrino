## ----message = FALSE----------------------------------------------------------
library(Rpadrino)

data(pdb)

problem_ipm <- pdb_make_proto_ipm(pdb, "aaa341") %>%
  pdb_make_ipm()

P <- problem_ipm$aaa341$sub_kernels$P
N <- solve(diag(nrow(P)) - P)

range(colSums(N))



## -----------------------------------------------------------------------------

kernel_formulae(problem_ipm)
vital_rate_exprs(problem_ipm)



## -----------------------------------------------------------------------------

problem_ipm <- pdb_make_proto_ipm(pdb, "aaa341") %>%
  pdb_make_ipm(addl_args = list(aaa341 = list(return_all_envs = TRUE)))


vr_funs <- vital_rate_funs(problem_ipm)

vr_funs

## -----------------------------------------------------------------------------

problem_proto <- pdb_make_proto_ipm(pdb, "aaa341")

vital_rate_exprs(problem_proto) <- pdb_new_fun_form(
  list(
    aaa341 = list(
      s = pmin(0.98, plogis(si + ss1 * size_1 + ss2 * size_1 ^ 2))
    )
  )
)

good_ipm <- pdb_make_ipm(problem_proto,
                         addl_args = list(aaa341 = list(return_all_envs = TRUE)))

vital_rate_funs(good_ipm)


P <- good_ipm$aaa341$sub_kernels$P
N <- solve(diag(nrow(P)) - P)

range(colSums(N))


