#' Firth's penalized-likelihood logistic regression with more decimal places of
#' p-value than \code{logistf} function in the R package \sQuote{logistf}
#'
#' Adapted from \code{logistf} in the R package \sQuote{logistf}, this is 
#' the same as \code{logistf} except that it provides more decimal places 
#' of p-value that would be useful for Genome-Wide Association Study (GWAS) 
#' or Phenome Wide Association Study (PheWAS).
#'
#' @param formula a formula object, with the response on the left of the
#' operator, and the model terms on the right. The response must be a vector
#' with 0 and 1 or FALSE and TRUE for the  outcome, where the higher value (1 or
#' TRUE) is modeled. It is possible to include contrasts, interactions, nested
#' effects, cubic or polynomial splines and all S features as well, e.g.
#' \code{Y ~ X1*X2 + ns(X3, df=4)}. From version 1.10, you may also include
#' offset() terms.
#' @param data a data.frame where the variables named in the formula can be
#' found, i. e. the variables containing the binary response and the covariates.
#' @param pl specifies if confidence intervals and tests should be based on the
#' profile penalized log likelihood (\code{pl=TRUE}, the default) or on the Wald
#' method (\code{pl=FALSE}).
#' @param alpha the significance level (1-\eqn{\alpha} the confidence level,
#' 0.05 as default).
#' @param control Controls Newton-Raphson iteration. Default is \cr
#' \code{control=logistf.control(maxstep, maxit, maxhs, lconv, gconv, xconv})
#' @param plcontrol Controls Newton-Raphson iteration for the estimation of the
#' profile likelihood confidence intervals. Default is \cr
#' \code{plcontrol=logistpl.control(maxstep, maxit,}
#' \code{maxhs, lconv, xconv, ortho, pr)}
#' @param firth use of Firth's penalized maximum likelihood (\code{firth=TRUE},
#' default) or the standard maximum likelihood method (\code{firth=FALSE}) for
#' the logistic regression. Note that by specifying \code{pl=TRUE} and
#' \code{firth=FALSE} (and probably a lower number of iterations) one obtains
#' profile likelihood confidence intervals for maximum likelihood logistic
#' regression parameters.
#' @param init specifies the initial values of the coefficients for the fitting
#' algorithm.
#' @param weights specifies case weights. Each line of the input data set is
#' multiplied by the corresponding element of \code{weights}.
#' @param plconf specifies the variables (as vector of their indices) for which
#' profile likelihood confidence intervals should be computed. Default is to
#' compute for all variables.
#' @param flic If TRUE, intercept is altered such that the predicted
#' probabilities become unbiased while keeping all other coefficients constant
#' @param model If TRUE the corresponding components of the fit are returned.
#' @param \dots Further arguments to be passed to logistf.
#'
#' @return
#' same as \code{logistf} except for providing more decimal places of p-value.
#'
#' @templateVar author choibeck
#' @template auth
#'
#' @references same as those provided in the R package \sQuote{logistf}.
#'
#' @examples
#' data(dataPheWAS)
#' fit <- Logistf(X264.3 ~ exposure + age + race + gender, data=dd)
#' summary(fit)
#' @export

`Logistf` <- 
function (formula, data, pl = TRUE, alpha = 0.05, control, plcontrol, 
    firth = TRUE, init, weights, plconf = NULL, flic = FALSE, 
    model = TRUE, ...) 
{
    if(!requireNamespace("logistf", quietly = TRUE)) {
      stop("Logistf requires the logistf package, please install it.",
        call. = FALSE)
    }
    ns <- loadNamespace('logistf')
    logistpl <- getFromNamespace('logistpl', ns)
    logistf.fit <- getFromNamespace('logistf.fit', ns)
    logistf.control <- getFromNamespace('logistf.control', ns)
    logistpl.control <- getFromNamespace('logistpl.control', ns)
    isspecnum <- getFromNamespace('isspecnum', ns)

    call <- match.call()
    extras <- list(...)
    call_out <- match.call()
    if (missing(control)) 
        control <- logistf.control()
    if (pl == TRUE & missing(plcontrol)) 
        plcontrol <- logistpl.control()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "weights", "na.action", "offset"), 
        names(mf), 0L)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.response(mf)
    n <- length(y)
    x <- model.matrix(mt, mf)
    k <- ncol(x)
    cov.name <- labels(x)[[2]]
    weight <- as.vector(model.weights(mf))
    offset <- as.vector(model.offset(mf))
    if (is.null(offset)) 
        offset <- rep(0, n)
    if (is.null(weight)) 
        weight <- rep(1, n)
    if (missing(init)) 
        init <- rep(0, k)
    if (is.null(plconf)) {
        if (isspecnum(cov.name, "(Intercept)")) {
            plconf <- NULL
            rest_plconf <- 1
        }
        else {
            plconf <- 1:k
            rest_plconf <- NULL
        }
    }
    else {
        rest_plconf <- (1:k)[-plconf]
    }
    if (cov.name[1] == "(Intercept)") {
        int <- 1
        coltotest <- 2:k
    }
    else {
        int <- 0
        coltotest <- 1:k
    }
    if (!is.null(extras$terms.fit)) {
        colfit <- eval(extras$terms.fit)
        matched <- match(colfit, cov.name)
        if (any(is.na(matched))) 
            stop(paste0("term(s): ", paste(colfit[is.na(matched)], 
                collapse = ", "), " not found in formula.\n"))
        nterms <- length(matched)
        colfit <- (1:k)[matched]
        call_out$terms.fit <- extras$terms.fit
        rest_plconf <- (1:k)[-matched]
    }
    else {
        colfit <- 1:k
        nterms <- k
    }
    fit.full <- logistf.fit(x = x, y = y, weight = weight, offset = offset, 
        firth, col.fit = colfit, init, control = control)
    fit.null <- logistf.fit(x = x, y = y, weight = weight, offset = offset, 
        firth, col.fit = int, init, control = control)
    if (fit.full$iter >= control$maxit) {
        warning(paste("Maximum number of iterations exceeded. Try to increase the number of iterations or alter step size by passing 'logistf.control(maxit=..., maxstep=...)' to parameter control"))
    }
    fit <- list(coefficients = fit.full$beta, alpha = alpha, 
        terms = colnames(x), var = fit.full$var, df = nterms - 
            int, loglik = c(fit.null$loglik, fit.full$loglik), 
        iter = fit.full$iter, n = sum(weight), y = y, formula = formula(formula), 
        call = call_out, conv = fit.full$conv)
    names(fit$conv) <- c("LL change", "max abs score", "beta change")
    beta <- fit.full$beta
    covs <- fit.full$var
    pi <- fit.full$pi
    fit$firth <- firth
    fit$linear.predictors <- as.vector(x %*% beta + offset)
    fit$predict <- fit.full$pi
    fit$hat.diag <- fit.full$Hdiag
    if (firth) 
        fit$method <- "Penalized ML"
    else fit$method <- "Standard ML"
    if (!is.null(extras$terms.fit)) {
        loc <- match(extras$terms.fit, fit$terms)
        var.red <- fit$var[loc, loc]
        vars <- diag(var.red)
        waldprob <- pchisq((beta[loc]^2/vars), 1, lower.tail = FALSE)
        wald_ci.lower <- as.vector(beta[loc] + qnorm(alpha/2) * 
            vars^0.5)
        wald_ci.upper <- as.vector(beta[loc] + qnorm(1 - alpha/2) * 
            vars^0.5)
    }
    else {
        vars <- diag(covs)
        waldprob <- pchisq((beta^2/vars), 1, lower.tail = FALSE)
        wald_ci.lower <- as.vector(beta + qnorm(alpha/2) * vars^0.5)
        wald_ci.upper <- as.vector(beta + qnorm(1 - alpha/2) * 
            vars^0.5)
    }
    fit$alpha <- alpha
    fit$conflev <- 1 - alpha
    if (pl) {
        betahist.lo <- vector(length(plconf), mode = "list")
        betahist.up <- vector(length(plconf), mode = "list")
        pl.conv <- matrix(0, length(plconf), 4)
        dimnames(pl.conv)[[1]] <- as.list(plconf)
        dimnames(pl.conv)[[2]] <- as.list(c("lower, loglik", 
            "lower, beta", "upper, loglik", "upper, beta"))
        LL.0 <- fit.full$loglik - qchisq(1 - alpha, 1)/2
        pl.iter <- matrix(0, k, 2)
        icount <- 0
        iters <- vector()
        for (i in plconf) {
            icount <- icount + 1
            inter <- logistpl(x, y, beta, i, LL.0, firth, -1, 
                offset, weight, plcontrol)
            fit$ci.lower[i] <- inter$beta
            pl.iter[i, 1] <- inter$iter
            betahist.lo[[icount]] <- inter$betahist
            pl.conv.lower <- t(inter$conv)
            inter <- logistpl(x, y, beta, i, LL.0, firth, 1, 
                offset, weight, plcontrol)
            fit$ci.upper[i] <- inter$beta
            pl.iter[i, 2] <- inter$iter
            betahist.up[[icount]] <- inter$betahist
            pl.conv.upper <- t(inter$conv)
            pl.conv[icount, ] <- cbind(pl.conv.lower, pl.conv.upper)
            fit.i <- logistf.fit(x, y, weight = weight, offset = offset, 
                firth, col.fit = (1:k)[-i], control = control)
            iters <- c(iters, fit.i$iter)
            fit$prob[i] <- pchisq(2 * (fit.full$loglik - fit.i$loglik), 1, lower.tail = FALSE)
            fit$method.ci[i] <- "Profile Likelihood"
        }
        fit$pl.iter <- pl.iter
        if (sum(fit$pl.iter >= plcontrol$maxit)) {
            notconv <- cov.name[apply(fit$pl.iter >= plcontrol$maxit, 
                1, sum) > 0]
            warning(paste("Nonconverged PL confidence limits: maximum number of iterations for variables:", 
                paste0(notconv, collapse = ", ")), " exceeded. Try to increase the number of iterations by passing 'logistpl.control(maxit=...)' to parameter plcontrol")
        }
        fit$betahist <- list(lower = betahist.lo, upper = betahist.up)
        fit$pl.conv <- pl.conv
        if (sum(iters >= control$maxit) > 0) {
            notconv <- cov.name[iters >= control$maxit]
            warning(paste("Maximum number of iterations for variables:", 
                paste0(notconv, collapse = ", ")), " exceeded. P-value may be incorrect. Try to increase the number of iterations by passing 'logistf.control(maxit=...)' to parameter control")
        }
    }
    else {
        fit$prob <- waldprob
        fit$method.ci <- rep("Wald", k)
        fit$ci.lower <- wald_ci.lower
        fit$ci.upper <- wald_ci.upper
    }
    for (i in rest_plconf) {
        fit$ci.lower[i] <- 0
        fit$ci.upper[i] <- 0
        fit$prob[i] <- 0
        fit$method.ci[i] <- "-"
    }
    names(fit$prob) <- names(fit$ci.upper) <- names(fit$ci.lower) <- names(fit$coefficients) <- dimnames(x)[[2]]
    if (flic) {
        if(!requireNamespace("formula.tools", quietly = TRUE)) {
          stop("using `flic` requires the formula.tools package, please install it.",
            call. = FALSE)
        }
        fit$flic <- TRUE
        lp_flic <- fit$linear.predictors - fit$coef[1]
        response <- formula.tools::lhs.vars(formula)
        fit_flic <- glm(as.formula(paste(response, paste("1"), 
            sep = " ~ ")), family = binomial(link = logit), data = data, 
            offset = lp_flic)
        W <- diag(fit_flic$fitted.values * (1 - fit_flic$fitted.values))
        tmp.var <- solve(t(x) %*% W %*% x)
        beta0.se <- sqrt(tmp.var[1, 1])
        fit$coefficients <- c(fit_flic$coef, fit$coef[-1])
        fit$ci.lower <- c(fit_flic$coef - beta0.se * 1.96, fit$ci.lower[-1])
        fit$ci.upper <- c(fit_flic$coef + beta0.se * 1.96, fit$ci.upper[-1])
        fit$linear.predictors <- fit_flic$linear
        fit$predict <- fit_flic$fitted
        fit$var <- tmp.var
        fit$method.ci[1] <- "Wald"
    }
    else fit$flic <- FALSE
    fit$control <- control
    if (model) 
        fit$model <- mf
    attr(fit, "class") <- c("logistf")
    fit
}
