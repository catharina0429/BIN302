Sys.setlocale("LC_TIME", "C")
options(scipen = 10000)
setwd("C:/Users/catha/Desktop/Jisoo/NMBU_DataScience/2021_Autumn/BIN302/Assignment/Finalpaper")

suppressMessages({
  library(readr)      # read file
  library(tidyverse)  # data manipulation
  library(skimr)      # data summary
  library(data.table) # data handling
  library(ggplot2);theme_set(theme_minimal())
  library(ggpubr)     # multiple plots in one figure
  library(pavo)       # spectral data
  library(pls)        # Partial Least Square Regression
})

Omega <- readr::read_csv("Omega_Assignment.csv")
NIR <- readr::read_csv("NIR_Assignment.csv")
# RAMAN <- readr::read_csv("Raman_Assignment.csv")

# Omega3 + NIR (Same as Practical) =======
ON <- left_join(Omega, NIR, by = "n") %>% 
  select(-c(3:4)) %>% 
  rename("ID" = "n")

df_NIR <- ON[, -c(1:2)]
df_NIRt <- t(df_NIR) # transpose
df_NIRt <- cbind(wl = row.names(df_NIRt),df_NIRt)
df_NIRt <- apply(df_NIRt, 2, as.numeric)
table(sapply(df_NIR, class), useNA = "always") # all is numeric
df_NIRt <- as.rspec(df_NIRt)
plot(df_NIRt, ylab = "Absorbance", 
     col = rainbow(3), main = "NIR Spectra with 200 obs")

df_NIRm = as.matrix(ON[, -c(1:2)])
df_NIRm <- data.frame(NIR = I(df_NIRm),
                      Omega3 = ON$Omega3 ,
                      ID = ON$ID)
class(df_NIRm$NIR) <- "matrix"
str(df_NIRm)

set.seed(1032)
idx <- caret::createDataPartition(df_NIRm$ID, p = 0.7)
tr_df_NIRm <- df_NIRm[idx$Resample1,]
te_df_NIRm <- df_NIRm[-idx$Resample1,]

df_NIRm_plsr <- plsr(Omega3 ~ NIR, data = tr_df_NIRm, 
                 ncomp = 10, validation = "CV", 
                 #jackknife = TRUE
                 )
summary(df_NIRm_plsr)
c <- pls::RMSEP(df_NIRm_plsr)
c <- as.data.frame(c$val) 
c <- as.data.frame(t(c)) %>% 
  slice(-1) %>% 
  mutate(Index = seq(1, 10, 1) %>% as.factor())

ggplot(c, aes(x = Index, y = CV)) + 
  geom_bar(aes(fill = Index),
           stat = "identity", show.legend = FALSE) + 
  geom_point(size = 3) + 
  labs(y = "RMSEcv",
       title = "RMSEcv by the number of Components") + 
  scale_x_discrete(labels=c(paste(seq(1, 10, 1), "comps"))) + 
  scale_fill_viridis_d(direction = -1) + 
  theme(axis.title.x = element_blank())

plot(RMSEP(df_NIRm_plsr), legendpos = "topright")
plot(df_NIRm_plsr, "loadings", comps = 1:3, legendpos = "topleft",labels = "numbers", xlab = "nm",ylim=c(-0.10,0.10))
abline(h = 0)

head(te_df_NIRm)

pred <- te_df_NIRm %>% 
  select(ID, Omega3) %>% 
  mutate(pred_NIR_comp2 = predict(df_NIRm_plsr, ncomp = 2,
                                  newdata = te_df_NIRm),
         pred_NIR_comp3 = predict(df_NIRm_plsr, ncomp = 3,
                                  newdata = te_df_NIRm))

ggpubr::ggscatter(data = pred, x = "Omega3", 
                  y = "pred_NIR_comp2", add = "reg.line") +
  stat_cor(label.y = 30, 
           aes(label = paste(..rr.label.., 
                             ..p.label.., sep = "~`,`~"))) +
  stat_regline_equation(label.y = 28)

ggpubr::ggscatter(data = pred, x = "Omega3", 
                  y = "pred_NIR_comp3", add = "reg.line") +
  stat_cor(label.y = 30, 
           aes(label = paste(..rr.label.., 
                             ..p.label.., sep = "~`,`~"))) +
  stat_regline_equation(label.y = 28)

# Omega3 + NIR (PCR) ====
# 그냥 pca를 했을 경우와
NIR_pca <- prcomp(NIR[, -1], 
                  center = TRUE, scale. = TRUE, rank. = 10)
summary(NIR_pca)
# 지수평활법을 써서 calibration을 한 다음 pca를 적용한 경우

