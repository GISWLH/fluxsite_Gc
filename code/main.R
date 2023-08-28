library(tidyverse)
library(data.table)

flux_cs <- read.csv("D:/acad3_211001/csv1/flux_sample.csv", header = TRUE)
IGBP_pars_new <- 
  read.csv("D:/acad3_211001/PML_R/PMLV2_parameters_TERRA_rm_pro_smooth_GEE.csv", header = TRUE) %>% t

IGBP_pars_new <- IGBP_pars_new[-1, ]
colnames(IGBP_pars_new) <- c('Alpha','Thelta', 'm', 'Am_25', 'D0', 'kQ', 'kA',
                             'S_sls', 'fER0', 'VPDmin', 'VPDmax', 'gsx', 'hc', 'ID')
IGBP_pars_new <- data.frame(IGBP_pars_new)
IGBP_pars_new$ID = 0:17
IGBP_pars_new$IGBPname <- rownames(IGBP_pars_new)
colnames(IGBP_pars_new)[14] = 'IGBPcode'
IGBP_pars_new[,1:13] <- as.numeric(unlist(IGBP_pars_new[,1:13]))
fulldata <- flux_cs %>% left_join(IGBP_pars_new)

output <- PML_site(#pars = pars, 
                   input = fulldata,
                   is_PMLV2 = TRUE)