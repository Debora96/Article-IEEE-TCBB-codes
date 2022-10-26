---
title: "R Notebook - Supplementary Figures: Collinearity, Correlation and Variable Ranking"
output: html_notebook
---


```{r}

if(!require("tidyverse")){install.packages("tidyverse")}
if(!require("DALEX")){install.packages("DALEX")}
if(!require("randomForest")){install.packages("randomForest")}
if(!require("ggfortify")){install.packages("ggfortify")}
if(!require("survival")){install.packages("survival")}
if(!require("cowplot")){install.packages("cowplot")}
if(!require("pals")){install.packages("pals")}
if(!require("colorspace")){install.packages("colorspace")}
if(!require("ggsci")){install.packages("ggsci")}

#install https://indrajeetpatil.github.io/ggstatsplot/index.html
# remotes::install_github("IndrajeetPatil/ggstatsplot")
# sudo apt-get install libgmp3-dev
# sudo apt-get install libmpfr-dev
install.packages("ggstatsplot")
library("ggstatsplot")

```

# Load Data ---- 

```{r}

GENESLUSCSELEC <- read_csv("../data/GENESLUSCSELEC.csv")
LUSC_clin <- read_csv("../data/LUSC_clin.csv")

datatime <- subset(LUSC_clin, select= c(submitter_id,time))
datagenes <- full_join(GENESLUSCSELEC, datatime, by = "submitter_id")
datagenes <- datagenes[complete.cases(datagenes[ , c('time')]), ] 

datagenes <- na.omit(datagenes)

datagenes[datagenes==-1] <- 0

sign<- c("PKHD1","CNTNAP","HCN1","ERICH3","RELN","DNAH8","PKHD1L1","DNAH5","KMT2D","FAM135B","SYNE1","SI","CDH10","PAPPA2","DAMTS12","RYR3","MUC17","PCDH15","PCLO","COL11A1","NAV3","SPTA1","FLG","XIRP2","ZFHX4","LRP1B","RYR2","CSMD3","MUC16","TP53","TTN")
sign <- intersect(sign, colnames(datagenes))

```

# Aalen's additive regression model for censored data

The additive model provides an alternative to Cox regression. It can provide detailed information about temporal influence of each of the covariates not available in Cox's model.


```{r}
library(survival)
  
aa_fit <-aareg(as.formula(paste0("Surv(time, vital_status) ~ " , paste(sign, collapse= " + "))),
               data = datagenes)
aa_fit

```


```{r}

vars<- c("PKHD1","CNTNAP","HCN1","ERICH3","RELN","DNAH8","PKHD1L1","DNAH5","KMT2D","FAM135B","SYNE1","SI","CDH10","PAPPA2","DAMTS12","RYR3","MUC17","PCDH15","PCLO","COL11A1","NAV3","SPTA1","FLG","XIRP2","ZFHX4","LRP1B","RYR2","CSMD3","MUC16","TP53","TTN")

variables <- rev(factor(vars, levels=vars))

#-- color decrease of lightness
#cols1 <- as.vector(glasbey(15))
#cols1 <- readhex(file = textConnection(paste(cols1, #collapse = "\n")),
#                 class = "RGB")
#transform to hue/lightness/saturation colorspace
#cols1 <- as(cols1, "HLS")
#additive decrease of lightness
#cols1@coords[, "L"] <- pmax(0, cols1@coords[, "L"] + 0.05)
#cols1 <- as(cols1, "RGB")
#cols1 <- hex(cols1)


p1 <- ggcoefstats(
  x = aa_fit,
  title = "Aalen's additive regression model",
  subtitle = "(for censored data)",
  only.significant = F,
  #point.args = list(color = "green", shape = 9),
  package = "pals",
  stats.label.args = list(
    max.time = 3,
    direction = "y",
    point.padding = 0.2,
    nudge_x = .15,
    nudge_y = .5,
    segment.curvature = -0.1,
    segment.angle = 10,
    segment.size  = 0.2,
    segment.linetype = 2
    # nudge_x = .15,
    # box.padding = 0.5,
    # nudge_y = 1,
    # segment.curvature = -0.1,
    # segment.ncp = 3,
    # segment.angle = 20
    ),
  palette = "polychrome",
  sort = "none", 
  k = 3
) 

p2 <- ggplot2::autoplot(aa_fit) +
  theme(legend.position="none")

p2$layers[[1]]$datagenes$variable <- factor(p2$layers[[1]]$datagenes$variable,
                                       levels= variables)
p2$layers[[2]]$datagenes$variable <- factor(p2$layers[[2]]$datagenes$variable,
                                       levels= variables)
p2$datagenes$variable <- factor(p2$datagenes$variable,
                                       levels= variables)


#p2 <- p2 +
#  scale_fill_manual(values=rev(cols1))

#Figure 5

pdf("fig5.pdf", width = 16, height = 9)
  cowplot::plot_grid(p1, p2, labels = "auto", nrow=1, rel_widths = c(0.7,1))
dev.off()

#rm(list = setdiff(ls(),c("sign", "data")))
```
