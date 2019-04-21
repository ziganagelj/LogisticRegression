

GetSeperateP <- function(model) {
  as.numeric(summary(model)$coefficients[2, 4])
}

GetStadijP <- function(data) {
  # data - genmerirani podatki
  # data <- stadij.data.H0
  
  p <- array(0,
             dim = c(m, 2, 4),
             dimnames = list(NULL,
                             c("U", "M"),
                             c('S1', 'S2', 'S3', 'S4')))
  
  for (i in 1:m) {
    vzorec <- rbind(data[y == 0][sample(nrow(data[y == 0]), 150),],
                    data[y == 1][sample(nrow(data[y == 1]), 150),])
    
    # M, multivariate 
    model <- glm(y ~ stadij, family=binomial(link='logit'), data = vzorec)
    p[i, 2, ] <- anova(model, test = 'LRT')$`Pr(>Chi)`[-1]
    
    # U vsak stadij posebaj
    model.S1 <- glm(y ~ S1, family=binomial(link='logit'), data = vzorec)
    model.S2 <- glm(y ~ S2, family=binomial(link='logit'), data = vzorec)
    model.S3 <- glm(y ~ S3, family=binomial(link='logit'), data = vzorec)
    model.S4 <- glm(y ~ S4, family=binomial(link='logit'), data = vzorec)
    
    p[i, 1, ] <- sapply(list(model.S1,
                             model.S2, 
                             model.S3,
                             model.S4),
                        GetSeperateP)
  }
  p
}

GetAlpha <- function(data) {
  # data - matrika k x 2 p vrednosti ki jih primerjamo
  # data <- dt.H0.p[, c(1, 2), 1]
  # ZAVRNEM U 
  zavrnjene_U <- mean(data[, 1] < 0.05)
  
  # ZAVRNEM M
  zavrnjene_M <- mean(data[, 2] < 0.05)
  
  # ZAVRNE U in ZAVRNE M
  zavrnjene_U_M <- mean(data[, 1] < 0.05 & data[, 2] < 0.05)
  
  # ZAVRNE U in NE ZAVRNE M
  zavrnjene_U_neM <- mean(data[, 1] < 0.05 & data[, 2] > 0.05)
  
  # Ne ZAVRNE U in ZAVRNE  M
  zavrnjene_neU_M <- mean(data[, 1] > 0.05 & data[, 2] < 0.05)
  
  c('P(zavrneU)' = zavrnjene_U,
    'P(zavrneM)' = zavrnjene_M,
    'P(zavrneM | zavrenU)' = round(zavrnjene_U_M/zavrnjene_U, 3), 
    'P(zavrneU | zavrenM)' = round(zavrnjene_U_M/zavrnjene_M, 3))
}


GetPlotUnivariate <- function(data, name) {
  for (v in unlist(dimnames(data)[3])) {
    plt <- ggplot(data.frame(data[, , v]), aes(x = U)) +
      geom_histogram(aes(y = ..count../m), breaks = seq(0, 1, 0.1), fill = '#71a0c5') +
      xlab("p-vrednost") + 
      ylab("Relativna frekvenca") +
      geom_vline(aes(xintercept = 0.05, color=("Stopnja značilnosti")), size=0.1 ,show.legend = FALSE) +
      theme_minimal()
    ggsave(paste0('./Images/', name, '_U_', v, '.pdf'), plt, width = 6, height = 6)
  }
}

GetPlotMultivariate <- function(data, name) {
  plt <- ggplot(data.frame(data[, , 1]), aes(x = M)) +
    geom_histogram(aes(y = ..count../m), breaks = seq(0, 1, 0.1), fill = '#71a0c5') +
    xlab("p-vrednost") + 
    ylab("Relativna frekvenca") +
    geom_vline(aes(xintercept = 0.05, color=("Stopnja značilnosti")), size=0.1 ,show.legend = FALSE) +
    theme_minimal()
  ggsave(paste0('./Images/', name, '_M.pdf'), plt, width = 6, height = 6)
}


GetResults <- function(data, name) {
  # UNIVARIATE
  GetPlotUnivariate(data = data, name = name)
  # MULTIVARIATE
  GetPlotMultivariate(data = data, name = name)
  
  apply(data, 3, GetAlpha)
  
}

GetZapletP <- function(data) {
  # data - genmerirani podatki
  # data <- zaplet.data.H0
  
  p <- array(0,
             dim = c(m, 2, 10),
             dimnames = list(NULL,
                             c("U", "M"),
                             c('Z1', 'Z2', 'Z3', 'Z4', 'Z5', 'Z6', 'Z7', 'Z8', 'Z9', 'Z10')))
  
  for (i in 1:m) {
    vzorec <- rbind(data[y == 0][sample(nrow(data[y == 0]), 150),],
                    data[y == 1][sample(nrow(data[y == 1]), 150),])
    
    # celoten model
    model <- glm(y ~ ., family=binomial(link='logit'), data = vzorec)
    
    # F, posamezne za full model
    p[i, 2, ] <-  as.numeric(summary(model)$coefficients[, 4])[-1] # brez intercepta
    
    # vsak stadij posebaj
    model.Z1 <- glm(y ~ Z1, family=binomial(link='logit'), data = vzorec)
    model.Z2 <- glm(y ~ Z2, family=binomial(link='logit'), data = vzorec)
    model.Z3 <- glm(y ~ Z3, family=binomial(link='logit'), data = vzorec)
    model.Z4 <- glm(y ~ Z4, family=binomial(link='logit'), data = vzorec)
    model.Z5 <- glm(y ~ Z5, family=binomial(link='logit'), data = vzorec)
    model.Z6 <- glm(y ~ Z6, family=binomial(link='logit'), data = vzorec)
    model.Z7 <- glm(y ~ Z7, family=binomial(link='logit'), data = vzorec)
    model.Z8 <- glm(y ~ Z8, family=binomial(link='logit'), data = vzorec)
    model.Z9 <- glm(y ~ Z9, family=binomial(link='logit'), data = vzorec)
    model.Z10 <- glm(y ~ Z10, family=binomial(link='logit'), data = vzorec)
    
    # S, za vsak model poseben 
    p[i, 1, ] <- sapply(list(model.Z1, 
                             model.Z2,
                             model.Z3, 
                             model.Z4, 
                             model.Z5,
                             model.Z6, 
                             model.Z7, 
                             model.Z8, 
                             model.Z9,
                             model.Z10), 
                        GetSeperateP)
  }
  p
}
