\name{TBEinf.glm.fit}
\alias{TBEinf.glm.fit}
\title{TBE infinity applied to GLM}
\description{\eqn{\text{TBE}_\infty} is applied to the GLM for estimation.}
\usage{
TBEinf.glm.fit(y, X, beta0 = NULL, tau = 1e-8, max.iter = 100)
}
\arguments{
\item{y}{Vector of responses \eqn{\mbox{\bf y}\in \mathbb{R}^{n}}.}
\item{X}{The design matrix \eqn{\mbox{\bf X}\in \mathbb{R}^{n\times p}}.}
\item{beta0}{Initial parameter vector \eqn{\boldsymbol \beta_0}. Default is \eqn{\boldsymbol \beta_0}=NULL
in which \eqn{\boldsymbol \beta_0} is given by \eqn{\boldsymbol \beta_0=(-0.01,-0.01,\cdots,-0.01)^\top}.}
\item{tau}{The threshold \eqn{\tau} for judging convergence. Default is 1e-8.}
\item{max.iter}{The maximum number of outer-loop iterations allowed. Default is 100.}
}
\value{
A list.
\itemize{
    \item beta - Estimate of coefficients of covariates.
    \item iter - Iterations.
    \item m - Estimate of mean obtained through link function.
    \item phi - Estimate of dispersion.
}
}
\details{
    The \eqn{\text{TBE}_\infty} GLM assumes a systematic component where 
    the linear predictor \eqn{\eta =\beta _{0}+\sum_{j=1}^{p}\beta _{j}x_{j}} 
    is linked to the mean \eqn{m} through a link function \eqn{\eta =-e^{-m}}.
    The estimation of \eqn{\boldsymbol \beta} is obtained through Iteratively 
    Reweighted Least Squares (IRLS). The estimate of \eqn{\varphi} is obtained 
    through mean deviance.
}
\references{
	Dunn, P. K. and Smyth, G. K. (2018).
   \emph{Generalized Linear Models with Examples in R}.
   Springer: Berlin/Heidelberg, Germany.
}

\examples{
set.seed(123) # For reproducibility
link_inverse <- function(eta) -log(- eta) # Inverse link function

n <- 200
p <- 9
pp <- p+1
beta_0 <- -0.1*(1:pp)
phi_0 <- 0.01

x <- matrix(runif(n * p), nrow = n, ncol = p)
X <- cbind(1, x)
eta_0 <- X %*% beta_0
m_0 <- link_inverse(eta_0)
y <- rnorm(n, mean = m_0, sd = phi_0*exp(m_0))

initial <- rep(-0.01, pp)
beta <- TBEinf.glm.fit(y, X, initial)$beta
iter <- TBEinf.glm.fit(y, X, initial)$iter
m <- TBEinf.glm.fit(y, X)$m
phi <- TBEinf.glm.fit(y, X)$phi
}
\seealso{\code{\link{TBEinf.Dev}} }
