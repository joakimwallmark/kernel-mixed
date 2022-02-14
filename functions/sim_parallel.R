library(foreach)
library(doParallel)

sim_parallel <- function(iter,
                   data_gen_FUN,
                   data_gen_FUN_args,
                   method_FUNs,
                   method_FUNs_args,
                   packages_vec = NULL,
                   export_vec = NULL,
                   res_filename) {
  no_methods <- length(method_FUNs)
  init_list <- vector("list", no_methods)
  cl <- makeCluster(detectCores(logical = F))
  registerDoParallel(cl)
  sim_res <- foreach(i = 1:iter,
                    .combine = "comb",
                    .packages = packages_vec,
                    .export = export_vec,
                    .init = init_list) %dopar% {
    data_i <- data_gen_FUN(i, data_gen_FUN_args)
    res_list <- vector("list", no_methods)
    for (meth in 1:length(method_FUNs)) {
      res_list[[meth]] <- method_FUNs[[meth]](data_i, i, method_FUNs_args[[meth]])
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
  for (i in 1:length(list2)) {
    comb_list[[i]] <- rbind(list1[[i]], list2[[i]])
  }
  return(comb_list)
}
