TARV_marginal_effect <- function(pheno, geno, form.adj, family)
{
  form <- as.formula(paste(as.character(form.adj), "+ geno[i, ]"))
  Z <- rep(NA, nrow(geno))
  for(i in 1:nrow(geno)) {
    model <- glm(form, data = pheno, family = family)
    Z[i] <- tail(summary(model)$coefficient, n = 1)[3]
  }
  return(Z)
}
