# cfd <- function(data, Tmax = NULL, time = NULL, idLabel = NULL, byrow = TRUE) {
#   out <- data
#   if (is.list(data) && (c("id", "time", "state") %in% names(data))) {
#     out <- dataframeToCfd(data)
#   }
#
#   if (is.matrix(data)) {
#     out <- matrixToCfd(data, time, idLabel, byrow)
#   }
#
#   out <- sortTime(out)
#   out <- remove_duplicated_states(out, keep.last = TRUE)
#
#   if (!is.null(Tmax)) {
#     out <- cut_cfd(out, Tmax)
#   }
#
#   return(out)
# }
#
#
# dataframeToCfd <- function(data) {
#   requiredColumns <- c("id", "time", "state")
#   isMissing <- !(requiredColumns %in% colnames(data))
#   if (any(isMissing)) {
#     stop(paste("Missing column(s):", paste(requiredColumns[isMissing], collapse = ", ")))
#   }
#
#   out <- data.frame(id = data$id, time = data$time, state = data$state)
#   class(out) <- c(class(out), "cfd")
#
#   out <- sortTime(out)
#   out <- remove_duplicated_states(out, keep.last = TRUE)
#
#   return(out)
# }
#
#
# matrixToCfd <- function(X, time = NULL, idLabel = NULL, byrow = TRUE) {
#   nInd <- ifelse(byrow, nrow(X), ncol(X))
#   nTime <- ifelse(byrow, ncol(X), nrow(X))
#
#   ## manage id labels
#   if (!is.null(idLabel)) {
#     if (length(idLabel) != nInd) {
#       stop(paste("idLabel must have the same length as the number of individuals:", nInd))
#     }
#
#     name <- idLabel
#   } else {
#     rowName <- rownames(X)
#     if (!is.null(idLabel)) {
#       name <- rowName
#     } else {
#       name <- seq_len(nInd)
#     }
#   }
#   id <- rep(seq_len(nInd), each = nTime)
#
#
#   ## manage time values
#   if (is.null(time)) {
#     time <- seq_len(nTime)
#   } else {
#     if (length(time) != nTime) {
#       stop(paste("time must have the same length as the number of time values:", nTime))
#     }
#   }
#
#   time2 <- rep(time, nInd)
#
#
#   ## manage state
#   state <- c()
#   if (byrow) {
#     state <- as.vector(t(as.matrix(X)))
#   } else {
#     state <- as.vector(as.matrix(X))
#   }
#
#   out <- data.frame(id = id, time = time2, state = state)
#   class(out) <- c(class(out), "cfd")
#
#   out <- sortTime(out)
#   out <- remove_duplicated_states(out, keep.last = TRUE)
#
#   return(out)
# }
#
# sortTime <- function(data) {
#   data <- do.call(rbind, lapply(split(data, data$id), function(x) {
#     ord <- order(x$time)
#     x[ord, ]
#   }))
#   rownames(data) <- NULL
#
#   return(data)
# }
