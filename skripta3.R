library(data.table)
library(ggplot2)

# load('env.RData')
set.seed(24)
k <- 100000
m <- 10000
# najrepj ti k + 2000 zato da povem o hkrathnem testirnaju hipotez, poltem damo vzorec na 10000, naredi loop za k, drugac eje komentar da ker subsamplamo iz velike populacije ni tok velika tezava hkratnega testiranja hipotez, drugace bi pa tudi to bil problem
# I STADIJ ----
# generirali bomo vecjo populacijo torej 5000 vzorcev, da potem subsamplamo in zagotovimo da je 150 operacij vsakega tipa

# POPULACIJSKI DELEŽI STADIJEV
stadij.verjetnosti  <- c(0.60, 0.25, 0.10, 0.05)
stadij <- sample(c(1, 2, 3, 4), size = k, prob = stadij.verjetnosti, replace = TRUE)
S <- as.matrix(data.frame(ifelse(stadij == 1, 1, 0), 
                          ifelse(stadij == 2, 1, 0), 
                          ifelse(stadij == 3, 1, 0), 
                          ifelse(stadij == 4, 1, 0)))
colnames(S) <- c('S1', 'S2', 'S3', 'S4')

# 1. DATA ----
# NICELNA DOMNEVA, torej regresijski koficient je 0

# linearna kobinacija z biasom
z.H0 = 0 + c(0, 0, 0, 0) %*% t(S)    

# verjetnosti za vsak vzorec glede na linearno kombinacijo
pr.H0 = 1/(1+exp(-z.H0))

# odvisna spremenljivka generirana iz binomske glede na verjetnost
y.H0 = rbinom(k, 1, pr.H0)

stadij.data.H0 <- data.table(y = y.H0, S, stadij = factor(stadij))
table(stadij.data.H0$y)

# ALTERNATVINA DOMENVA, torej gregresijski koficient je različen od 0, v našem, primeru 2
z.HA = 0 + c(0.1, 0.3, 0.6, -1.2) %*% t(S)    
pr.HA = 1/(1+exp(-z.HA))
y.HA = rbinom(k, 1, pr.HA)

stadij.data.HA <- data.table(y = y.HA, S, stadij = factor(stadij))
table(stadij.data.HA$y)

# 2. MODEL ----
# U - 4x univariat
# M - multivariate

stadij.p.H0 <- GetStadijP(stadij.data.H0)
stadij.p.HA <- GetStadijP(stadij.data.HA)

# 3. ANALIZA ----
# pooling stadijev  vs celotna informacija

# 3.a H0 ----
table.stadij_H0 <- GetResults(data = stadij.p.H0, name = 'stadij_H0')

# 3.b HA ----
table.stadij_HA <- GetResults(data = stadij.p.HA, name = 'stadij_HA')


# II ZAPLET ----
# generirali bomo vecjo populacijo torej 5000 vzorcev, da potem subsamplamo in zagotovimo da je 150 operacij vsakega tipa

# POPULACIJSKI DELEŽI STADIJEV
zaplet.verjetnost  <- c(0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.3, 0.4, 0.5, 0.6)
Z <- t(replicate(k, rbinom(10, 1, zaplet.verjetnost)))
colnames(Z) <- c('Z1', 'Z2', 'Z3', 'Z4', 'Z5', 'Z6', 'Z7', 'Z8', 'Z9', 'Z10')

# 1. DATA ----
# NICELNA DOMNEVA, torej regresijski koficient je 0
# linearna kobinacija z biasom
z.H0 = 1 + c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0) %*% t(Z)  

# verjetnosti za vsak vzorec glede na linearno kombinacijo
pr.H0 = 1/(1+exp(-z.H0))

# odvisna spremenljivka generirana iz binomske glede na verjetnost
y.H0 = rbinom(k, 1, pr.H0)

zaplet.data.H0 <- data.table(y = y.H0, Z)
table(zaplet.data.H0$y)

# ALTERNATVINA DOMENVA, torej gregresijski koficient je različen od 0, v našem, primeru 2
z.HA = 0 + c(-3.5, 3, -2.5, 2, -1.5, 1, -0.8, 0.6, -0.4, 0.2) %*% t(Z)     
pr.HA = 1/(1+exp(-z.HA))
y.HA = rbinom(k, 1, pr.HA)

zaplet.data.HA <- data.table(y = y.HA, Z)
table(zaplet.data.HA$y)

# 2. MODEL ----
# U - 10x univariat
# M - multivariate
zaplet.p.H0 <- GetZapletP(zaplet.data.H0)
zaplet.p.HA <- GetZapletP(zaplet.data.HA)

# 3. ANALIZA ----
# pooling stadijev  vs celotna informacija

# 3.a H0 ----
table.zaplet_H0 <- GetResults(data = zaplet.p.H0, name = 'zaplet_H0')

# 3.b HA ----
table.zaplet_HA <- GetResults(data = zaplet.p.HA, name = 'zaplet_HA')

# III TABLES ----
table.list <- list(stadij_H0 = table.stadij_H0,
                   stadij_HA = table.stadij_HA,
                   zaplet_H0 = table.zaplet_H0,
                   zaplet_HA = table.zaplet_HA)

params.pop.stadij <- data.frame(stadij.verjetnosti)
rownames(params.pop.stadij) <- paste('Stadij', 1:4)

params.pop.zaplet <- data.frame(zaplet.verjetnost)
rownames(params.pop.zaplet) <- paste('Zaplet', 1:10)

colnames(params.pop.stadij) <- colnames(params.pop.zaplet) <- 'verjetnost'

params.pop <- list(stadij = params.pop.stadij,
                   zaplet = params.pop.zaplet)

params.hypothesis <- data.frame(stadij = c(0, 0.1, 0.3, 0.6, -1.2, NA, NA, NA, NA, NA, NA),
                                zaplet = c(0, -3.5, 3, -2.5, 2, -1.5, 1, -0.8, 0.6, -0.4, 0.2))

rownames(params.hypothesis) <- paste0('beta',0:10)

params.hypothesis <- params.hypothesis



saveRDS(table.list, "./table_list.rds")
saveRDS(params.pop, "./params_pop.rds")
saveRDS(params.hypothesis, "./params_hypothesis.rds")
