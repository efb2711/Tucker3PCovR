\name{plot.t3pcovr}
\alias{plot.t3pcovr}
\title{Biplots}
\description{Program for producing joint biplots and interactive biplots.}
\usage{
plot.t3pcovr(sol,type, fixmode = NULL, fixunit = NULL, laba, labb, labc, labbc = NULL))}
\arguments{
  \item{sol}{A t3pcovr.out object}
  \item{type}{Type of biplot: 1- joint biplot, 2- interactive biplot. Mandatory}
  \item{fixmode}{Mode for which one unit is to be chosen (1 for   A-mode, 2 for B-mode, 3 for C-mode). Mandatory for joint        biplots}
  \item{fixunit}{Number of component for which joint plot is      desired. Mandatory for joint biplots}
  \item{laba}{Vector of length n containing the labels of the A   -mode entities. Mandatory.}
  \item{labb}{Vector of length m containing the labels of the B   -mode entities. Mandatory for joint biplot.}
  \item{labc}{Vector of length p containing the labels of the C   -mode entities. Mandatory for joint biplot.}
  \item{labbc}{Vector of length mxp containing the labels of the   BxC-mode entities. Mandatory for interactive biplot.}
}

\references{\cite{Gabriel, K. R. (1971). The biplot graphic display of matrices with application to principal component analysis. Biometrika, 58 , 453-467.}

\cite{Carlier, A., & Kroonenberg, P. (1996). Decompositions and biplots in three-way correspondence analysis. Psychometrika, 61 , 355-373.}}

\author{Elisa Frutos-Bernal (\email{efb@usal.es})}
\examples{
  data(X3way)
  data(X2way)
  laba = c("Person1","Person2","Person3","Person4","Person5","Person6","Person7","Person8")

  labb = c("Other anger","Shame","Love","Sorrow","Fear","Guilt","Self(anger)")

  labc = c("Quarrelling with someone","Partner leaves you","Someone is telling lies about you","Giving a bad speech","Failing a test","Writing a bad paper")

  labbc = c("S1R1","S1R2","S1R3","S1R4","S1R5","S1R6","S1R7","S2R1","S2R2","S2R3","S2R4","S2R5","S2R6","S2R7","S3R1","S3R2","S3R3","S3R4","S3R5","S3R6","S3R7","S4R1","S4R2","S4R3","S4R4","S4R5","S4R6","S4R7","S5R1","S5R2","S5R3","S5R4","S5R5","S5R6","S5R7","S6R1","S6R2","S6R3","S6R4","S6R5","S6R6","S6R7")

  sol<- t3pcovr (X2way, X3way, n = 8, m= 7, p=6, q=10, r1=2, r2   =2, r3 =2, conv=1e-6, OriginalAlfa=0.5,AlternativeLossF=1       ,nRuns = 100, StartSeed = 10)
  #joint biplot for the third mode, first component
  plot(sol,1,3,1,laba,labb,labc)
  #interactive biplot
  plot(sol,2,NULL,NULL,laba,NULL,NULL,labbc)
}
