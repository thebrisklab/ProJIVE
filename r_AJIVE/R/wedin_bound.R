#'
#' Estimates the wedin bound for a data matrix with the resampling procedure.
#'
#' returns min(max(||E tilde{V}||, ||E^T tilde{U}||) / sigma_min(widetilde{A}), 1) from equation (7)
#'
#' @param X Matrix. The data matrix.
#' @param SVD List. The SVD decomposition of X.
#' @param signal_rank Integer. The estimated signal rank of X.
#'
#' @return The the wedin bound samples.
get_wedin_bound_samples <- function(X, SVD, signal_rank, num_samples=1000){

    # resample for U and V
    U_perp <- matrix(SVD[['u']][ , -(1:signal_rank)], nrow = nrow(SVD[['u']]))
    U_sampled_norms <- wedin_bound_resampling(X=X,
                                              perp_basis=U_perp,
                                              right_vectors=FALSE,
                                              num_samples=num_samples,
                                              signal_rank)

    V_perp <- matrix(SVD[['v']][ , -(1:signal_rank)], nrow = nrow(SVD[['v']]))
    V_sampled_norms <- wedin_bound_resampling(X=X,
                                              perp_basis=V_perp,
                                              right_vectors=TRUE,
                                              num_samples=num_samples,
                                              signal_rank)

    sigma_min <- SVD[['d']][signal_rank]
    wedin_bound_samples <- mapply(function(u, v)  asin(min(max(u, v)/sigma_min, 1)), U_sampled_norms, V_sampled_norms)

    wedin_bound_samples
}


#' Resampling procedure for the wedin bound
#'
#' @param X Matrix. The data matrix.
#' @param perp_basis Matrix. Either U_perp or V_perp: the remaining left/right singluar vectors of X after estimating the signal rank.
#' @param right_vectors Boolean. Right multiplication or left multiplication.
#' @param num_samples Integer. Number of vectors selected for resampling procedure.
#' @param signal_rank Integer. The estimated signal rank of X.
wedin_bound_resampling <- function(X, perp_basis, right_vectors, num_samples=1000, signal_rank){

    rank <- signal_rank
    resampled_norms <- rep(0, num_samples)

    for(s in 1:num_samples){

        sampled_col_index <- sample.int(n=dim(perp_basis)[2],
                                        size=rank,
                                        replace=TRUE)


        perp_resampled <- perp_basis[ , sampled_col_index]

        if(right_vectors){
            resampled_projection <- X %*% perp_resampled
        } else{
            resampled_projection <- t(perp_resampled) %*% X
        }

        # operator L2 norm
        resampled_norms[s] <- svd(resampled_projection)[['d']][1]
    }

    resampled_norms
}

