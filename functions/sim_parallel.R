library(foreach)
library(doParallel)

sim_parallel <- function(iter,
                         data_gen_fun,
                         data_gen_fun_args,
                         method_funs,
                         methods_funs_args,
                         packages_vec = NULL,
                         export_vec = NULL,
                         res_filename) {
  no_methods <- length(method_funs)
  init_list <- vector("list", no_methods)
  cl <- makeCluster(detectCores(logical = FALSE))
  registerDoParallel(cl)
  sim_res <- foreach(
    i = 1:iter,
    .combine = "comb",
    .packages = packages_vec,
    .export = export_vec,
    .init = init_list
  ) %dopar% {
    data_i <- data_gen_fun(i, data_gen_fun_args)
    res_list <- vector("list", no_methods)
    for (meth in seq_along(method_funs)) {
      res_list[[meth]] <- method_funs[[meth]](data_i, i, methods_funs_args[[meth]])
    }
    res_list
  }
  stopCluster(cl)
  save(sim_res, file = res_filename)
  return(sim_res)
}

# List1 is list from all previous iterations (.init in foreach at first)
# List2 is from this iteration
comb <- function(list1, list2) {
  comb_list <- vector("list", length(list2))
  for (i in seq_along(list2)) {
    comb_list[[i]] <- rbind(list1[[i]], list2[[i]])
  }
  return(comb_list)
}
