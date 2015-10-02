DGLM_norm <- function(m.form, d.form, indata, maxiter = 20, conv = 1e-6)
{
    X.mean <- model.matrix(m.form, data = indata)
    X.disp <- model.matrix(d.form, data = indata)
    y.name <- all.vars(m.form)[1]
    y <- indata[,y.name]
    w <- rep(1, nrow(indata))
    convergence <- 1
    iter <- 0
    while (convergence > conv & iter < maxiter) {
      iter <- iter + 1
      w.old <- w
      glm1 <- lm(y ~ . - 1, weights = w, data = data.frame(X.mean))
      res <- resid(glm1)
      q <- hatvalues(glm1)
      y2 <- (res ^ 2)/(1 - q)
      glm2 <- glm(y2 ~ . - 1, family = Gamma(link = log), weights = (1 - q)/2, data = data.frame(X.disp))
      w <- 1/fitted(glm2)
      convergence <- (max(abs(w.old - w)) + (summary(glm1)$sigma - 1) )
    }
    return(list(mean = glm1, disp = glm2, iter = iter))
}
