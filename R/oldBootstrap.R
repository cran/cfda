# 
# 
# bootstrapEncoding <- function(data, basisobj, computeCI = TRUE, nBootstrap = 50, propBootstrap = 0.5, nCores = max(1, ceiling(detectCores()/2)), verbose = TRUE, ...)
# {
#   t1 <- proc.time()
#   ## check parameters
#   checkData(data)
#   checkDataBeginTime(data)
#   checkDataEndTmax(data)
#   if(!is.basis(basisobj))
#     stop("basisobj is not a basis object.")
#   if(any(is.na(nCores)) || (length(nCores) != 1) || !is.whole.number(nCores) || (nCores < 1))
#     stop("nCores must be an integer > 0.")
#   if(any(is.na(nBootstrap)) || (length(nBootstrap) != 1) || !is.whole.number(nBootstrap) || (nBootstrap < 1))
#     stop("nBootstrap must be an integer > 0.")
#   if(any(is.na(propBootstrap)) || !is.numeric(propBootstrap) || (length(propBootstrap) != 1) || (propBootstrap > 1) || (propBootstrap < 0))
#     stop("nCores must be an integer > 0.")
#   ## end check
# 
#   if(verbose)
#     cat("######### Compute encoding #########\n")
# 
#   # change state as integer
#   out <- stateToInteger(data$state)
#   data$state = out$state
#   label <- out$label
#   rm(out)
# 
#   # refactor labels as 1:nbId
#   uniqueId <- unique(data$id)
#   nId <- length(uniqueId)
#   id2 <- refactorCategorical(data$id, uniqueId, seq_along(uniqueId))
#   uniqueId2 <- unique(id2)
# 
#   nCores <- min(max(1, nCores), detectCores()-1)
# 
#   Tmax <- max(data$time)
#   K <- length(label$label)
#   nBasis <- basisobj$nbasis
# 
#   if(verbose)
#   {
#     cat(paste0("Number of individuals: ", nId, "\n"))
#     cat(paste0("Number of states: ", K, "\n"))
#     cat(paste0("Basis type: ", basisobj$type, "\n"))
#     cat(paste0("Number of basis functions: ", nBasis, "\n"))
#     cat(paste0("Number of cores: ", nCores, "\n"))
#   }
# 
# 
#   V <- computeVmatrix(data, uniqueId2, id2, basisobj, K, nCores, verbose, ...)
# 
#   Uval <- computeUmatrix(data, uniqueId2, id2, basisobj, K, nCores, verbose, ...)
# 
# 
#   fullEncoding <- computeEncoding(Uval, V, K, nBasis, uniqueId, label, verbose)
# 
# 
# 
#   # outEnc <- list()
#   # t3 <- proc.time()
#   # if(verbose)
#   #   cat("---- Compute Bootstrap Encoding:\n")
#   #
#   # for(i in 1:nBootstrap)
#   # {
#   #   cat("*")
#   #   idToKeep <- sample(nId, sizeBootstrap, replace = TRUE)
#   #
#   #   outEnc[[i]] = computeEncoding(Uval[idToKeep, ], V[idToKeep, ], K, nBasis, idToKeep, label, verbose = FALSE)
#   #
#   #   outEnc[[i]] = c(outEnc[[i]] , list(basisobj = basisobj))
#   #   class(outEnc[[i]]) = "fmca"
#   # }
#   #
#   # # reorder alpha, pc such that representation have the same sign for each bootstrap sample
#   # outEnc = unifySign(outEnc)
#   #
#   #
#   # alpha <- lapply(seq_along(outEnc[[1]]$alpha), function(harm) Reduce("+", lapply(outEnc, function(x) x$alpha[[harm]]))/nBootstrap)
#   # invF05vec <- computeMeanInvF05vec(outEnc)
#   #
#   # t4 <- proc.time()
#   # if(verbose)
#   #   cat(paste0("\nDONE in ", round((t4-t3)[3], 2), "s\n"))
# 
#   # out <- list(out = outEnc, alpha = alpha, pc = V%*%invF05vec, V = V, invF05vec = invF05vec, basisobj = basisobj, nBootstrap = nBootstrap)
# 
#   if(computeCI)
#   {
#     bootEncoding <- computeBootStrapEncoding(Uval, V, K, nBasis, label, nId, propBootstrap, nBootstrap, verbose)
#     varAlpha <- computeVarianceAlpha(bootEncoding)
#   }
# 
#   out <- c(fullEncoding, list(V = V, basisobj = basisobj, bootstrap = bootEncoding))
#   class(out) = c("fmca", "fmcaBootstrap")
# 
#   t2 <- proc.time()
# 
#   if(verbose)
#     cat(paste0("Run Time: ", round((t2-t1)[3], 2), "s\n"))
# 
#   return(out)
# }
# 
# 
# # compute the variance of encoding
# #
# # a = sum_i alpha_i * phi_i 
# # Var(a) = sum_i var(alpha_i) * phi_i^2 + sum_{i<j}  2 * phi_i * phi_j * cov(alpha_i, alpha_j)
# # computeVarianceEncoding <- function(bootEncoding, harm = 1)
# # {
# #   nBasis <- bootEncoding[[1]]$basisobj$nbasis
# #   phi <- fd(diag(nBasis), bootEncoding[[1]]$basisobj)
# #   nState <- ncol(bootEncoding[[1]]$alpha[[harm]])
# #   
# #   funcVar <- list()
# #   for(iState in 1:nState)
# #   {
# #     varAlpha <- var(t(sapply(bootEncoding, function(x) x$alpha[[harm]][, iState])))
# #     funcVar[[iState]] <- sum(diag(varAlpha) * phi^2)
# #     for(i in 1:(nBasis-1))
# #     {
# #       for(j in i:nBasis)
# #       {
# #         funcVar[[iState]] = funcVar[[iState]] + 2 * phi[i] * phi[j] * varAlpha[i, j]
# #       }
# #     }
# #   }
# #   
# #   return(funcVar)
# # }
# # 
# # 
# # 
# # computeVarianceEncodingb <- function(bootEncoding, harm = 1, nx = 200)
# # {
# #   nBasis <- bootEncoding[[1]]$basisobj$nbasis
# #   phi <- fd(diag(nBasis), bootEncoding[[1]]$basisobj)
# #   nState <- ncol(bootEncoding[[1]]$alpha[[harm]])
# #   
# #   timeVal <- seq(bootEncoding[[1]]$basisobj$rangeval[1],
# #                bootEncoding[[1]]$basisobj$rangeval[2],
# #                length = nx)
# #   
# #   
# #   Phi <- matrix(nrow = length(timeVal), ncol = nBasis)
# #   for(i in 1:nBasis)
# #     Phi[, i] = eval.fd(timeVal, phi[i])
# #   
# #   funcVar <- list()
# #   for(iState in 1:nState)
# #   {
# #     varAlpha <- var(t(sapply(bootEncoding, function(x) x$alpha[[harm]][, iState])))
# #     funcVar[[iState]] = rep(NA, nx)
# #     for(j in seq_along(timeVal))
# #     {
# #       # print(dim(Phi[,j, drop = FALSE]))
# #       # print(dim(varAlpha))
# #       funcVar[[iState]][j] <- Phi[j,, drop = FALSE] %*% varAlpha %*% t(Phi[j, , drop = FALSE])
# #       
# #     }
# #     
# #   }
# #   
# #   return(funcVar)
# # }
# 
# 
# 
# 
# # 
# # 
# # plot.fmcaBootstrap <- function(x, harm = 1, coeff = 2)
# # {
# #   encod <- get_encoding(x, harm = harm, fdObject = FALSE, nx = 128)
# #   matplot(encod$x, encod$y, type = "l", lty = 1)
# #   
# #   
# #   for(i in 1:ncol(encod$y))
# #   {
# #     lines(encod$x, encod$y[,i] + as.vector(sqrt(eval.fd(encod$x, x$varAlpha[[harm]][[i]])) * coeff), lwd = 2, col = i, lty = 2)
# #     lines(encod$x, encod$y[,i] - as.vector(sqrt(eval.fd(encod$x, x$varAlpha[[harm]][[i]])) * coeff), lwd = 2, col = i, lty = 2)
# #   }
# # }
# 
# 
# 
# computeVarianceEncoding <- function(varAlpha, basisobj, harm = 1, nx = 200)
# {
#   nBasis <- basisobj$nbasis
#   phi <- fd(diag(nBasis), basisobj)
#   nState <- length(varAlpha[[harm]])
#   
#   timeVal <- seq(basisobj$rangeval[1], basisobj$rangeval[2], length = nx)
#   
#   Phi <- matrix(nrow = length(timeVal), ncol = nBasis)
#   for(i in 1:nBasis)
#     Phi[, i] = eval.fd(timeVal, phi[i])
#   
#   funcVar <- list()
#   for(iState in 1:nState)
#   {
#     funcVar[[iState]] = rep(NA, nx)
#     for(j in seq_along(timeVal))
#     {
#       funcVar[[iState]][j] <- Phi[j,, drop = FALSE] %*% varAlpha[[harm]][[iState]] %*% t(Phi[j, , drop = FALSE])
#     }
#     
#   }
#   
#   return(funcVar)
# }
# 
# 
# plot2 <- function(x, harm = 1, coeff = 2, nx = 128)
# {
#   encod <- get_encoding(x, harm = harm, fdObject = FALSE, nx = nx)
#   matplot(encod$x, encod$y, type = "l", lty = 1)
#   
#   variance <- computeVarianceEncoding(x$varAlpha, x$basisobj, harm = harm, nx = nx)
#   
#   for(i in 1:ncol(encod$y))
#   {
#     lines(encod$x, encod$y[,i] + sqrt(variance[[harm]][[i]]) * coeff, lwd = 2, col = i, lty = 2)
#     lines(encod$x, encod$y[,i] - sqrt(variance[[harm]][[i]]) * coeff, lwd = 2, col = i, lty = 2)
#   }
# }
# 
# 
# 
# 
