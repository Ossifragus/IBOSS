#' The IBOSS method
#'
#' This function implements the IBOSS method for the input covariate Z and response vector Y. It returns a list with elements: beta, the least squares estimate based on the subdata; se, the standard errors; sigma, variance estimate for the error term, index, index of the subdata.
#' @useDynLib IBOSS
#' @param Z the input covariate matrix or covariate vector
#' @param Y the response vector
#' @param k the subdata size
#' @param int.adj whether to calculate the adjusted estimate of the intercept. It is \code{TRUE} by default. 
## ' @keywords iboss.od 
#' @export
#' @examples
#' library(IBOSS)
#' library(mvtnorm)
#' beta.true  <- rep(1, 51)
#' d <- length(beta.true) - 1
#' corr  <- 0.5
#' sigmax <- matrix(corr, d, d) + diag(1-corr, d)
#' n <- 5000
#' k <- 100
#' set.seed(0)
#' X  <- rmvt(n, sigmax, 2)
#' mu  <- beta.true[1] + c(X %*% beta.true[-1])
#' Y  <- mu + rnorm(n, 0, 3)
#' fit <- iboss.od(X, Y, k)
#' beta.od <- fit$beta
#' beta.odia <- fit$beta0.adj


iboss.od <- function(Z, Y, k, int.adj="TRUE") {
    call <- match.call()
    ## method <- match.arg(method)
    nd <- dim(Z)
    if (is.null(nd)) {
        n <- length(Z)
        d <- 1
        if (k / 2 != k %/% 2)
            warning("k/d/2 is not an integer; its foor is used as r.")
        r <- max(k %/% 2, 1)
        idx.od <- .Call("getIdx", as.integer(r), as.numeric(Z),
                        PACKAGE="IBOSS")
        x.od <- cbind(1, Z[idx.od])
        if (int.adj)
            Zbar <- mean(Z)
    }
    else if (is.matrix(Z)) {
        n <- nd[1]
        d <- nd[2]
        if (k / d / 2 != k %/% d %/% 2)
            warning("k/d/2 is not an integer; its floor is used as r.")
        r <- k %/% d %/% 2
        idx.od <- .Call("getIdx", as.integer(r), as.numeric(Z[,1]),
                        PACKAGE="IBOSS")
        for(j in 2:d) {
            tmp <- .Call("getIdxR", as.integer(r), as.numeric(Z[,j]),
                         as.integer(idx.od), PACKAGE="IBOSS")
            idx.od <- c(idx.od, tmp)
        }
        x.od <- cbind(1, Z[idx.od,])
        if (int.adj)
            Zbar <- colMeans(Z)
    }
    y.od <- Y[idx.od]
    eig.od <- eigen(t(x.od) %*% x.od)
    iI.od <- eig.od$vectors %*% (t(eig.od$vectors)/eig.od$values)
    beta.od <- c(iI.od %*% t(x.od) %*% y.od)
    if (int.adj) {
        Ybar <- mean(Y)
        beta.odia <- Ybar - sum(Zbar * beta.od[-1])
    }
    res.od <- y.od - c(x.od %*% beta.od)
    sigma.od <- sum(res.od^2) / (n - d)
    se.od <- diag(iI.od) * sigma.od
    result <- list(beta=beta.od, se=se.od, sigma=sigma.od,
                   index=idx.od, beta0.adj=beta.odia)
    return(result)
}
