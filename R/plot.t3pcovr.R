plot.t3pcovr <- function  (sol, type, fixmode = NULL, fixunit = NULL, laba, labb, labc, labbc = NULL)
{
  r1 = ncol(sol$A)
  n = nrow(sol$A)
  r2 = ncol(sol$B1)
  m = nrow(sol$B1)
  r3 = ncol(sol$C)
  p = nrow(sol$C)
  checkinput=1
  if (type !=1 & type != 2)
  {
    cat(" ",fill=TRUE)
    cat("Type of biplot should be 1 (jointbiplot) or 2 (interactive biplot)" , fill=TRUE )
    cat(" ",fill=TRUE)
    checkinput=0
  }
  if (type ==1 & (fixmode != 1 && fixmode != 2 && fixmode != 3))
  {
    cat(" ",fill=TRUE)
    cat("Fix mode should be 1, 2 or 3" , fill=TRUE )
    cat(" ",fill=TRUE)
    checkinput=0
  }
  if (checkinput == 1)
  {
  if (type == 1) ## Jointbiplot
  {
    if (fixmode == 1)
    {
        Gmat = matrix(t(sol$H[fixunit, ]), r2, r3)
        Smat = sol$B1 %*% Gmat %*% t(sol$C)
        SVD = svd(Smat)
        fit = SVD$d[1]^2 + SVD$d[2]^2
        fit = fit/SUM(Smat)$ssq * 100
        Bmat = (m/p)^0.25 * SVD$u[, 1:2] %*% diag(SVD$d[1:2])^0.5
        Cmat = (p/m)^0.25 * SVD$v[, 1:2] %*% diag(SVD$d[1:2])^0.5
        xmax = max(Bmat[,1], Cmat[,1]) + 0.1
        xmin = min(Bmat[,1], Cmat[,1]) - 0.1
        ymax = max(Bmat[,2], Cmat[,2]) + 0.1
        ymin = min(Bmat[,2], Cmat[,2]) - 0.1
        datos = data.frame(rbind(data.frame(Bmat, type = "B",label = labb), data.frame(Cmat, type = "C", label = labc)))
        colnames(datos)[1:2] <- c("d1", "d2")
        datos$type = factor(datos$type)
        jointbiplot = ggplot(data = datos, aes(x = d1, y = d2, group = "type"))+theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
          labs(title='',x='PC1',y='PC2')+ geom_hline(yintercept = 0,linetype="dashed") + geom_vline(xintercept = 0,linetype="dashed")
        jointbiplot <- jointbiplot + scale_x_continuous(limits = c(xmin,xmax)) + scale_y_continuous(limits = c(ymin,ymax))
        jointbiplot <- jointbiplot + geom_segment(xend = 0, yend = 0,data = subset(datos, type == "B"))
        jointbiplot <- jointbiplot + geom_text(size=3, data = subset(datos, type == "C"),
                                               aes(label = label), show.legend = FALSE)
        jointbiplot <- jointbiplot + geom_text(size=3, data = subset(datos, type == "B"),
                                               aes(label = label), show.legend = FALSE)+ coord_equal()

    }
    if (fixmode == 2)
    {
      K = permnew(sol$H, r1, r2, r3)
      Gmat = matrix(t(K[fixunit, ]), r3, r1)
      Smat = sol$C %*% Gmat %*% t(sol$A)
      SVD = svd(Smat)
      fit = SVD$d[1]^2 + SVD$d[2]^2
      fit = fit/SUM(Smat)$ssq * 100
      Cmat = (p/n)^0.25 * SVD$u[, 1:2] %*% diag(SVD$d[1:2])^0.5
      Amat = (n/p)^0.25 * SVD$v[, 1:2] %*% diag(SVD$d[1:2])^0.5
      xmax = max(Amat[,1], Cmat[,1]) + 0.1
      xmin = min(Amat[,1], Cmat[,1]) - 0.1
      ymax = max(Amat[,2], Cmat[,2]) + 0.1
      ymin = min(Amat[,2], Cmat[,2]) - 0.1
      datos = data.frame(rbind(data.frame(Amat, type = "A",label = laba), data.frame(Cmat, type = "C", label = labc)))
      colnames(datos)[1:2] <- c("d1", "d2")
      datos$type = factor(datos$type)
      jointbiplot = ggplot(data = datos, aes(x = d1, y = d2, group = "type"))+theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
        labs(title='',x='PC1',y='PC2')+ geom_hline(yintercept = 0,linetype="dashed") + geom_vline(xintercept = 0,linetype="dashed")
      jointbiplot <- jointbiplot + scale_x_continuous(limits = c(xmin,xmax)) + scale_y_continuous(limits = c(ymin,ymax))
      jointbiplot <- jointbiplot + geom_segment(xend = 0, yend = 0,data = subset(datos, type == "C"))
      jointbiplot <- jointbiplot + geom_text(size=3, data = subset(datos, type == "A"),
                                             aes(label = label), show.legend = FALSE)
      jointbiplot <- jointbiplot + geom_text(size=3, data = subset(datos, type == "C"),
                                             aes(label = label), show.legend = FALSE)+ coord_equal()

    }
    if (fixmode == 3)
    {
      K = permnew(sol$H, r1, r2, r3)
      K = permnew(K, r2, r3, r1)
      Gmat = matrix(t(K[fixunit, ]), r1, r2)
      Smat = sol$A %*% Gmat %*% t(sol$B1)
      SVD = svd(Smat)
      fit = SVD$d[1]^2 + SVD$d[2]^2
      fit = fit/SUM(Smat)$ssq * 100
      Amat = (p/n)^0.25 * SVD$u[, 1:2] %*% diag(SVD$d[1:2])^0.5
      Bmat = (n/p)^0.25 * SVD$v[, 1:2] %*% diag(SVD$d[1:2])^0.5
      xmax = max(Amat[,1], Bmat[,1]) + 0.1
      xmin = min(Amat[,1], Bmat[,1]) - 0.1
      ymax = max(Amat[,2], Bmat[,2]) + 0.1
      ymin = min(Amat[,2], Bmat[,2]) - 0.1
      datos = data.frame(rbind(data.frame(Amat, type = "A",label = laba), data.frame(Bmat, type = "B",label = labb)))
      colnames(datos)[1:2] <- c("d1", "d2")
      datos$type = factor(datos$type)
      jointbiplot = ggplot(data = datos, aes(x = d1, y = d2, group = "type"))+theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
        labs(title='',x='PC1',y='PC2')+ geom_hline(yintercept = 0,linetype="dashed") + geom_vline(xintercept = 0,linetype="dashed")
      jointbiplot <- jointbiplot + scale_x_continuous(limits = c(xmin,xmax)) + scale_y_continuous(limits = c(ymin,ymax))
      jointbiplot <- jointbiplot + geom_segment(xend = 0, yend = 0,data = subset(datos, type == "B"))
      jointbiplot <- jointbiplot + geom_text(size=3, data = subset(datos, type == "A"),
                                             aes(label = label), show.legend = FALSE)
      jointbiplot <- jointbiplot + geom_text(size=3, data = subset(datos, type == "B"),
                                             aes(label = label), show.legend = FALSE)+ coord_equal()
    }
    print(jointbiplot)
  }
  if (type ==2)
  {
    Matriz_BC <- t(sol$H %*% kronecker( t(sol$C) , t(sol$B1) ))
    xmax = max(Matriz_BC[,1], sol$A[,1],sol$B2[,1]) + 0.1
    xmin = min(Matriz_BC[,1], sol$A[,1],sol$B2[,1]) - 0.1
    ymax = max(Matriz_BC[,2], sol$A[,2],sol$B2[,2]) + 0.1
    ymin = min(Matriz_BC[,2], sol$A[,2],sol$B2[,2]) - 0.1


    datos = data.frame(rbind
                       (data.frame(sol$A, type = "persons", label = laba),
                         data.frame(Matriz_BC, type = "reactions-situations", label = labbc)))
    colnames(datos)[1:2] <- c("d1", "d2")
    datos$type = factor(datos$type)
    interactivebiplot = ggplot(data = datos, aes(x = d1, y = d2, group = "type"))+theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
      labs(title='',x='PC1',y='PC2')+ geom_hline(yintercept = 0,linetype="dashed") + geom_vline(xintercept = 0,linetype="dashed")
    interactivebiplot <- interactivebiplot + scale_x_continuous(limits = c(xmin,xmax)) + scale_y_continuous(limits = c(ymin,ymax))
    interactivebiplot <- interactivebiplot + geom_segment(xend = 0, yend = 0,data = subset(datos, type == "reactions-situations"))+ coord_equal()
    interactivebiplot <- interactivebiplot + geom_text(size=3, data = subset(datos, type == "persons"),
                                                       aes(label = label), show.legend = FALSE)
    interactivebiplot <- interactivebiplot + geom_text(size=3, data = subset(datos, type == "reactions-situations"),
                                                       aes(label = label), show.legend = FALSE)
    print(interactivebiplot)
  }
  }
}
