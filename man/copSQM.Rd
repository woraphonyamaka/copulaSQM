\name{copSQM}
\alias{copSQM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Copula based Stochastic Frontier Quantile Model
}
\description{
Estimating Copula based Stochastic Frontier Quantile Model
}
\usage{
copSQM(Y=Y,X=X,family=1,tau=0.5,RHO=0.5,LB=-0.99,UB=0.99)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{Y}{ vector of dependent variable}
\item{X}{ matrix of independent variable}
\item{family}{ Copula function eg. Gaussain=1, Student-t=2 (see, Vinecopula package)}
\item{tau}{ The qauntile level range between 0-1}
\item{RHO}{ The initail value of the copula parameter}
\item{LB}{  The lower bound of the copula parameter}
\item{UB}{ The upper bound of the copula parameter}
}
\details{
The Copula based Stochastic Frontier model of (Pipitpojanakarnet al.(2016) is the new class the technical efficiency measure. Methodologically, the model concern about dependency between two-sided error term and one-sided inefficiency. This model becomes more flexible to the outlier and it can measure the relationship between output and input levels across efficiency quantiles. In addition, this model also provides the different slopes of parameters describing the production of Asian countries rather than average value.
}
\value{

\item{result }{The result contain the estimated parameters, standard errors, t-stat, and p-value}
\item{AIC }{Akaiki Information Criteria}
\item{BIC }{Bayesian Information Criteria}
\item{Loglikelihood }{Maximum Log-likelihood function}

}
\references{
Pipitpojanakarn, V., Maneejuk, P., Yamaka, W., & Sriboonchitta, S. (2016, November). Analysis of agricultural production in Asia and measurement of technical efficiency using copula-based stochastic frontier quantile model. In International Symposium on Integrated Uncertainty in Knowledge Modelling and Decision Making (pp. 701-714). Springer, Cham.

Pipitpojanakarn, V., Yamaka, W., Sriboonchitta, S., & Maneejuk, P. (2017). Frontier Quantile Model Using a Generalized Class of Skewed Distributions. Advanced Science Letters, 23(11), 10737-10742.
}
\author{
Woraphon Yamaka
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(truncnorm)
library(mvtnorm)
library("VineCopula")
library("frontier")
library(ald)
library("LaplacesDemon")

# example included in FRONTIER 4.1 (cross-section data)
data(front41Data)
attach(front41Data)
# Cobb-Douglas production frontier
cobbDouglas <- sfa( log(output)~log(capital)+log(labour),data=front41Data)
summary(cobbDouglas)


# Select familty  copula upper and lower bouubd ( look at Vinecopula package)
# family=1   # 1 is Gaussian, 2 is Student-t, 3 is Clayton and so on....

#Gaussian (-.99, .99)
#Student t (-.99, .99)
#Clayton (0.1, Inf)
Y=log(output)
X=cbind(log(capital),log(labour))
model=copSQM(Y=Y,X=X,family=1,tau=0.5,RHO=0.5,LB=-0.99,UB=0.99)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
