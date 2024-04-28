\name{TBEinf.glm.pre}
\alias{TBEinf.glm.pre}
\title{TBE infinity applied to GLM}
\description{\eqn{\text{TBE}_\infty} is applied to the GLM for prediction.}
\usage{
TBEinf.glm.pre(X, beta)
}
\arguments{
\item{X}{new data \eqn{x}}
\item{beta}{Parameter vector \eqn{\boldsymbol \beta}. Generally estimated values.}
}
\value{
The prediction of response \eqn{y}.
}
\details{
    \code{TBEinf.glm.pre} predicts \eqn{y} by \eqn{\hat{y}=-\log(- \mbox{\bf x}^\top \beta)}.
}
\examples{
set.seed(123) # For reproducibility
link_inverse <- function(eta) -log(- eta) # Inverse link function

n <- 800
p <- 9
pp <- p+1
beta_0 <- -0.1*(1:pp)
phi_0 <- 0.01

x <- matrix(runif(n * p), nrow = n, ncol = p)
X <- cbind(1, x)
eta_0 <- X %*% beta_0
m_0 <- link_inverse(eta_0)
y <- rnorm(n, mean = m_0, sd = phi_0*exp(m_0))

beta <- TBEinf.glm.fit(y, X)$beta
phi <- TBEinf.glm.fit(y, X)$phi

x_new <- runif(pp)
y_pre <- TBEinf.glm.pre(x_new, beta)
}