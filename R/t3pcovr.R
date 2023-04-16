t3pcovr <- function( X , Y3D , n , m , p , q , r1 , r2 , r3 , conv , OriginalAlfa , AlternativeLossF , nRuns , StartSeed)
{
  set.seed( StartSeed )
  #
  checkinput = 1

  if( r1 > min(n,m*p,q) )
  {
    cat(" ",fill=TRUE)
    cat("rank1 should be an integer between 1 and " , min(i,l) , fill=TRUE )
    cat(" ",fill=TRUE)
    checkinput=0
  }

  if( r2 > min(m,n*p) )
  {
    cat(" ",fill=TRUE)
    cat("rank2 should be an integer between 1 and " , min(j,i*k) , fill=TRUE)
    cat(" ",fill=TRUE)
    checkinput=0
  }

  if( r3 > min(p,n*m) )
  {
    cat(" ",fill=TRUE)
    cat("rank3 should be an integer between 1 and " , min(k,i*j) , fill=TRUE)
    cat(" ",fill=TRUE)
    checkinput=0
  }
  if ( (r1 > r2*r3) | (r2 > r1*r3) | (r3 > r1*r2) )
  {
    cat(" ",fill=TRUE)
    cat("None of the ranks can be larger than the products of the other two (e.g., rank1 > rank2*rank3 is not allowed)",fill=TRUE)
    cat(paste(r1,r2,r3))
    cat(" ",fill=TRUE)
    checkinput=0
  }

  if ( (OriginalAlfa < 0) || (OriginalAlfa >= 1) )
  {
    cat(" ",fill=TRUE)
    cat("OriginalAlfa should be between 0 and 1 (but not 0 or 1)",fill=TRUE)
    cat(" ",fill=TRUE)
    checkinput=0
  }

  if ( (AlternativeLossF !=0) && (AlternativeLossF != 1) )
  {
    cat(" ",fill=TRUE)
    cat("AlternativeLossF should be 0 or 1",fill=TRUE)
    cat(" ",fill=TRUE)
    checkinput=0
  }

  if ( checkinput == 1 )
  {

    X = as.matrix(X)
    Y = as.matrix(Y3D)
    Y = matrix(Y , n , m * p) # P slices are concatenated horizontally
    ssq3D = sum(Y ^ 2)
    ssq2D = sum(X ^ 2)

    if (AlternativeLossF == 1)
    {
      Alfa = (OriginalAlfa * ssq3D) / ((OriginalAlfa * ssq3D) + ((1 - OriginalAlfa) * ssq2D))
    } else
    {
      Alfa = OriginalAlfa
    }
    BestA = matrix(0 , n , r1)
    BestB1 = matrix(0 , m , r2)
    BestC = matrix(0 , p , r3)
    BestH = array(0 , cbind(r1, r2, r3)) # CoreArray
    BestB2 = matrix(0 , q , r1)
    BestIter = -9999
    BestLoss = 999999999999
    cputime = system.time({
      FitValues = matrix(0 , 1 , nRuns + 1)
      nIterValues = matrix(0 , 1 , nRuns + 1)
      for (run in 1:nRuns + 1)
      {
        # initialize A, B1 and C (orthonormal)
        if (run  == 1)
        {
          # rational starts via eigendecomposition (gives orthonormal starting values)
          EIG = eigen(cbind(X, Y) %*% t(cbind(X, Y)))
          A = EIG$vectors[, 1:r1]
          rm(EIG)

          Z = permnew(Y , n , m , p)		# yields m x p x n array
          EIG = eigen(Z %*% t(Z))
          B1 = EIG$vectors[, 1:r2]
          rm(EIG)

          Z = permnew(Z , m , p , n)		# yields p x n x mrray
          EIG = eigen(Z %*% t(Z))
          C = EIG$vectors[, 1:r3]
          rm(EIG, Z)
        }
        else
        {
          # random start (orthonormal)
          A = orth(matrix(rnorm(n * r1 , 0 , 1) , n , r1))
          B1 = orth(matrix(rnorm(m * r2 , 0 , 1) , m , r2))
          C = orth(matrix(rnorm(p * r3 , 0 , 1) , p , r3))
        }
        # Calculate initial Core
        Z = permnew(t(A) %*% Y , r1 , m , p)
        Z = permnew(t(B1) %*% Z , r2 , p , r1)
        H = permnew(t(C) %*% Z , r3 , r1 , r2) # H (r1 x r2r3)
        rm(Z)

        # Calcule initial B2 (is in general not orthogonal !!! )
        B2 = t(solve(t(A) %*% A) %*% t(A) %*% X)

        # Evaluate f
        Model3D = A %*% H %*% kronecker(t(C) , t(B1))
        Model2D = A %*% t(B2)
        f = (1 - Alfa) * sum((Y - Model3D) ^ 2)  +  Alfa * sum((X - Model2D) ^ 2)
        iter = 0
        fold = f + (2 * conv * f)
        V = cbind((sqrt(Alfa) * X) , (sqrt(1 - Alfa) * Y))
        while ((fold - f) > (f * conv))
        {
          iter = iter + 1
          fold = f

          U = t(cbind((sqrt(Alfa) * t(B2)), (
            sqrt(1 - Alfa) * H %*% kronecker(t(C) , t(B1))
          )))

          # update B1 (orthonormal)
          Z = permnew(Y , n , m , p)
          Z = permnew(Z , m , p , n)
          Z = permnew(t(C) %*% Z , r3 , n , m)
          Z = permnew(t(A) %*% Z , r1 , m , r3)			 # yields m x r3 x r1 array
          B1 = qr.Q(qr(Z %*% (t(Z) %*% B1)) , complete = FALSE)
          rm(Z)

          # update C (orthonormal)
          Z = permnew(t(A) %*% Y , r1 , m , p)
          Z = permnew(t(B1) %*% Z , r2 , p , r1)			 # yields p x r1 x r2 array
          C = qr.Q(qr(Z %*% (t(Z) %*% C)) , complete = FALSE)
          rm(Z)

          # Update H (Core)
          Z = permnew(t(A) %*% Y , r1 , m , p)
          Z = permnew(t(B1) %*% Z , r2 , p , r1)
          H = permnew(t(C) %*% Z , r3 , r1 , r2)
          rm(Z)

          A = t(ginv(t(U) %*% U) %*% t(U) %*% t(V))

          # Update B2 (not necessarily orthogonal !!)
          B2 = t(solve(t(A) %*% A) %*% t(A) %*% X)


          # Evaluate f
          Model3D = A %*% H %*% kronecker(t(C) , t(B1))
          Model2D = A %*% t(B2)
          f = (1 - Alfa) * sum((Y - Model3D) ^ 2)  +  Alfa * sum((X - Model2D) ^ 2)
        }   #end of while-loop (alternating part of the algorithm)
        FitValues[run] = f
        nIterValues[run] = iter

        if (f < BestLoss)
        {
          BestA = A
          BestB1 = B1
          BestC = C
          BestB2 = B2
          BestH = H
          BestLoss = f
          BestIter = iter
        }
      } # end of for-loop (nRuns)
    }) # end of system.time
    # compute "intrinsic eigenvalues"
    La = BestH %*% t(BestH)
    J = permnew(BestH , r1 , r2 , r3)
    Lb = J %*% t(J)
    J = permnew(J , r2 , r3 , r1)
    Lc = J %*% t(J)

    # Compute BOF
    Model3D = BestA %*% BestH %*% kronecker(t(BestC) , t(BestB1))
    Model2D = BestA %*% t(BestB2)
    BOF3D = sum((Model3D - Y) ^ 2)
    BOF2D = sum((Model2D - X) ^ 2)

    FitPercentage = (((1 - Alfa) * (sum(BestH ^ 2) / ssq3D)) + (Alfa  * (sum(Model2D ^
                                                                               2) / ssq2D))) * 100
    FitPercentage3D = 100 * ((ssq3D - BOF3D) / ssq3D)
    FitPercentage2D = 100 * ((ssq2D - BOF2D) / ssq2D)

    out = list()
    out$Info = list()
    out$Info$nRows = n
    out$Info$nColumns3D = m
    out$Info$nSlices = p
    out$Info$nColumns2D = q
    out$Info$RankVector = cbind(r1 , r2 , r3)
    out$Info$TolPercentage = conv
    out$Info$OriginalAlfa = OriginalAlfa
    out$Info$Alfa = Alfa
    out$Info$AlternativeLossF = AlternativeLossF
    out$Info$nRuns = nRuns

    out$A = BestA
    out$B1 = BestB1
    out$C = BestC
    out$H = BestH
    out$B2 = BestB2
    out$LossWeighted = BestLoss
    out$LossUnweighted = BOF3D + BOF2D
    out$FitPercentage = FitPercentage
    out$FitPercentage3D = FitPercentage3D
    out$FitPercentage2D = FitPercentage2D
    out$nIter = BestIter
    out$FitValues = FitValues
    out$nIterValues = nIterValues
    out$La = La
    out$Lb = Lb
    out$Lc = Lc
    out$Fit3D = BOF3D
    out$Fit2D = BOF2D
    out$Fit = (Alfa * BOF3D) + ((1 - Alfa) * BOF2D)
    out$Fit3Dsize = BOF3D / (n * m * p)
    out$Fit2Dsize = BOF2D / (n * q)
    out$Fitsize = out$Fit / ((n * m * p) + (n * q))
    out$CpuTime = cputime[1]
    out$TimeSeconds = round(cputime[1] , 2)

    class(out) = "t3pcovr"
    return(out)
  }
}


