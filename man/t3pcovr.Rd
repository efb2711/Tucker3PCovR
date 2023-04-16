\name{t3pcovr}
\alias{t3pcovr}
\title{Estimation of Tucker3 Principal Covariates Regression parameters}
\description{Tucker3-PCovR model reduces the predictors (the two-way matrix) to a few components and predict the criterion (the three-way array) by using these components as predictors at the same time that the three way data is fitted through the Tucker3 model.}
\usage{
t3pcovr(X, Y3D,...)}
\arguments{
  \item{X}{(n x q): 2D-data matrix}
  \item{Y3D}{(n x m x p): 3D-data matrix (represented as a 3D-array)}
  \item{n}{number of elements of first mode of 3D/2D (the common mode: rows)}
  \item{m}{number of elements of second mode of 3D (columns 3D)}
  \item{p}{number of elements of third mode of 3D (slabs)}
  \item{q}{number of elements of second mode of 2D (columns 2D)}
  \item{r1}{Number of extracted components for the A-mode}
  \item{r2}{Number of extracted components for the B-mode}
  \item{r3}{Number of extracted components for the C-mode}
  \item{conv}{value for convergence (tolerance value)}
  \item{OriginalAlfa}{(0-1): importance that degree reduction and prediction have in the analysis}
  \item{AlternativeLossF}{Using the alternative loss function? 0 = no (use original loss function: weighted SSQ;               weighted met alfa) 1 = yes (use weighted loss function with scaled SSQ: scaled by the SSQ in X and y )}
  \item{nRuns}{Number of runs}
  \item{StartSeed}{Seed for the analysis}
}
\value{
  \item{A}{Component matrix for the A-mode)}
  \item{B1}{Component matrix for the B-mode}
  \item{C}{Component matrix for the C-mode}
  \item{H}{Matrized core array (frontal slices)}
  \item{B2}{Loading matrix of components (components x predictors)}
  \item{...}{Further arguments}
}

\references{\cite{De Jong, S., & Kiers, H. A. (1992). Principal covariates regression: Part I. Theory. Chemometrics and Intelligent Laboratory Systems , 155-164.}

\cite{Marlies Vervloet, Henk A. Kiers, Wim Van den Noortgate, Eva Ceulemans (2015). PCovR: An R Package for Principal Covariates Regression. Journal of Statistical Software, 65(8), 1-14. URL http://www.jstatsoft.org/v65/i08/.}

\cite{Smilde, A. K., Bro, R., & Geladi, P. (2004). Multi-way analysis with applications in the chemical sciences. Chichester, UK: Wiley.}}

\author{Elisa Frutos-Bernal (\email{efb@usal.es})}
\examples{
  data(X3way)
  data(X2way)
  t3pcovr (X2way, X3way, n = 8, m= 7, p=6, q=10, r1=2, r2=2, r3 =2, conv=1e-6, OriginalAlfa=0.5,AlternativeLossF=1,nRuns = 100, StartSeed = 10)
}
