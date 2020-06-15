\name{TE}
\alias{TE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Technical efficiency measure.
}
\description{
Computing and plotting the technical efficiency.
}
\usage{
TE(theta,Y,X,family=1,tau)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{theta}{ The estimated parameters form the model}
 \item{Y}{ Vector of dependent variable}
\item{X}{ Matrix of independent variable}
\item{family}{ Copula function eg. Gaussain=1, Student-t=2 (see, Vinecopula package)}
\item{tau}{ The qauntile level range between 0-1}
}
\details{
Computing and plotting the technical efficiency.
}
\value{

item{Output }{Technical efficiency series}
item{plot }{Plot of Technical efficiency}

}
\references{
Pipitpojanakarn, V., Maneejuk, P., Yamaka, W., & Sriboonchitta, S. (2016, November). Analysis of agricultural production in Asia and measurement of technical efficiency using copula-based stochastic frontier quantile model. In International Symposium on Integrated Uncertainty in Knowledge Modelling and Decision Making (pp. 701-714). Springer, Cham.
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

#EX: Plot the technical efficiency
tau=0.5
te1=TE(model$result[,1],Y=Y,X=X,family=1,tau)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
