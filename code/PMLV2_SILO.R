# references
# 1. Zhang, Y., Kong, D., Gan, R., Chiew, F. H. S., McVicar, T. R., Zhang, Q., & Yang, Y. (2019). Coupled estimation 
#    of 500 and 8-day resolution global evapotranspiration and gross primary production in 2002â€?2017. Remote Sensing 
#    of Environment, 222, 165-182. https://doi.org/10.1016/j.rse.2018.12.031
# 2. Gan, R., Zhang, Y., Shi, H., Yang, Y., Eamus, D., Cheng, L., et al. (2018). Use of satellite leaf area index 
#    estimating evapotranspiration and gross assimilation for Australian ecosystems. Ecohydrology, 11(5), e1974. 
#    https://doi.org/10.1002/eco.1974


# dependecy
library(data.table)
library(tidyverse)
# default parameters for 17-type IGBP landcover 
# while some are not specifid in RSE paper (no fluxnet site)
# be catious that parameters change for AU EBF PFT

IGBP_pars_new <- 
  read.csv("D:/acad3_211001/PML_R/PMLV2_parameters_TERRA_rm_pro_smooth_GEE.csv") %>% t
IGBP_pars_new <- IGBP_pars_new[-1, ]
colnames(IGBP_pars_new) <- c('Alpha','Thelta', 'm', 'Am_25', 'D0', 'kQ', 'kA',
                               'S_sls', 'fER0', 'VPDmin', 'VPDmax', 'gsx', 'hc', 'ID')
IGBP_pars_new <- data.frame(IGBP_pars_new)
IGBP_pars_new$ID = 0:17
IGBP_pars_new$IGBPname <- rownames(IGBP_pars_new)
colnames(IGBP_pars_new)[14] = 'IGBPcode'

#=====================================================
# PML GLOBAL PARAMETERS
# define for all circumstances  
#=====================================================
Gsc         = 0.0820  # solar constant in unit MJ m-2 min-1
as          = 0.25    # parameter Rs/Ra=as+bs*n/N calibration from our solar radiation measurement
bs          = 0.50    # parameter Rs/Ra=as+bs*n/N
alfa        = 0.23    # surface albedo of grass
alfa_forest = 0.22    # surface albedo of forest
alfa_crop   = 0.14    # surface albedo of crop
kmar        = 0.40    # von Karman's constant 0.40 
Zob         = 15      # m making sure higher than hc
Cp          = 1.0164  # 4.2 * 0.242 specific heat at constant pressure 1.013  [J g-1 0C-1
epsl        = 0.622   # ratio molecular weight of water vapour/dry air

# PML_v1 parameters for Gc 
kQ          = 0.4488  # extinction coefficient
kA          = 0.7     # the attenuation of net all-wave irradicance typically about 0.6-0.8 (Denmend 1976 Kelliher FM et al. (1995))
Q50         = 30      # the value of absorbed PAR when gs=gsx/2 W/m2
D0          = 0.7     # the value of VPD when stomtal conductance is reduced  kpa 

#====================================================================================
# MAIN FUNCTION
# require parameters and inputs (data.frame or data.table, equal-size list)
# 1. diff with GEE javascript version in input,
#    use actual vapor pressure (ea), instead of q (specific humidity)
# 2. an R version, requires 
#====================================================================================
#' PML modelling at site scale
#'
#' @param pars a vector that has 13 parameters attached for the PFT (IGBP code), must in sequence 
#'  Alpha, 
#'  Thelta, 
#'  m, 
#'  Am_25,
#'  D0, 
#'  kQ, 
#'  kA, 
#'  S_sls, 
#'  fERO, 
#'  VPD_min
#'  VPDmax, 
#'  gsx, 
#'  hc, 
#' @param input  a data.frame or data.table, equal-size list that includes all following parameters,
#' be cautious about the unit of each input
#'  co2, umol mol-1, need to be converted when it is ppm
#'  ea, kPa, actual vapor pressure (or q, kg/kg, if ea not directly available) 
#'  Pa, kPa, air pressure 
#'  U2,  m/s, wind speed at 2m height 
#'  Tmax, degC
#'  Tmin, degC
#'  Tavg, degC
#'  Rln, W/m2/s, not MJ/m2/d, long-wave radiation 
#'  Rs,  W/m2/s, solar radiation (shortwave radiation / direct)
#'  Prcp, mm/d, precipitation
#'  Albedo, 0-1, surface albedo
#'  Emiss, edmissivity, ~0.97
#'  LAI,  0-7, m^2 m^-2
#'  Landcover, IGBP landcover ID, 0-17
#'  
#' @param is_PMLV2 bool, default using PML_V2, otherwise V1
#'
#' @return
#' @export
#'
#' @examples
#' 
PML_site <- function(input, is_PMLV2 = TRUE){
  #==============================================================
  # PML pars and input
  #==============================================================
  
  # parameters, see RSE paper 
  
  Alpha <- input[['Alpha']]
  Thelta <- input[['Thelta']]
  m <- input[['m']]
  Am_25 <- input[['Am_25']]
  D0 <- input[['D0']]
  kQ <- input[['kQ']]
  kA <- input[['kA']]
  S_sls <- input[['S_sls']]
  fER0 <- input[['fER0']]
  VPDmin <- input[['VPDmin']]
  VPDmax <- input[['VPDmax']]
  gsx <- input[['gsx']]
  hc <- input[['hc']]
  LAIref <- 5
  
  
  # PML INPUTS
  
  Ca <-  input[['co2']]   # umol mol-1, need to be converted when it is ppm
  # ea, kPa, actual vapure pressure 
  if(!is.null(input[['ea']])){
    ea <- input[['ea']]     
  }
  else{
    if(!is.null(input[['q']])){
      q <- input[['q']]
      # actual vapure pressure converted from q (specific humidity)
      # https://www.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html, Eq-17
      ea <- q * p / (0.622 + 0.378 * q);
    }else{
      LE1 <- input[['LE']]
      #error('there is no ea or q in the input dataframe')
    }
  }
  p <-  input[['Pa']]    # kPa, air pressure 
  u2 <-  input[['U2']]    # m/s, wind speed at 2m height 
#  Tmax <-  input[['Tmax']]  # degC
#  Tmin <-  input[['Tmin']]  # degC
  Tavg <-  input[['Tavg']]  # degC
  Rln <-  input[['Rln']]   # W/m2/s, not MJ/m2/d, long-wave radiation 
  Rs <-  input[['Rs']]    # W/m2/s, solar radiation (shortwave radiation / direct)
  P <- input[['Prcp']] # mm/d, precipitation
  
  # modis
  albedo <-  input[['Albedo']]# 0-1, surface albedo
  emiss <-  input[['Emiss']] # 0-1, edmissivity, ~0.97
  LAI <-  input[['LAI']]   # 0-7, m^2 m^-2
  land <- input[['Landcover']] # IGBP landcover ID, 0-17
  
  
  #==============================================================
  # PML model 
  #==============================================================
  lambda <-  2500; # latent heat of vaporization, 2500 [J g-1]  at 25 degC
  lambda <- Tavg * -2.2 + lambda # xu.edit 2.4
  
#  es_tmax <- vapor_pressure(Tmax)
#  es_tmin <- vapor_pressure(Tmin)
   es_tavg <- vapor_pressure(Tavg)
#  es <- (es_tmax + es_tmin)/2
  
#  VPD <- es - ea
#  VPD <- rmax(VPD, 0.001)
  VPD <- input[['VPD']]
  rou_a <- 3846 * p / (Tavg + 273.15)
  
  gama <- Cp * p / (0.622 * lambda)
  
  slop <- 4098 * es_tavg / (Tavg + 237.3) ^ 2
  
  Stefan <- 4.903e-9
  
  Rns <- (1 - albedo) * Rs
  RLout <- (emiss * Stefan * (Tavg + 273.15) ^ 4) / 0.0864
  
  Rnl <- Rln - RLout #Rnl is net longwave radiation (input)
  
  Rn <- (Rns + Rnl) %>% rmax(0)
  
  # units convert: http://www.egc.com/useful_info_lighting.php
  PAR <- (Rs * 0.45) %>% rmax(0)
  
  fvpd_gc <- 1/(1 + VPD / D0)
  
  fvpd <- (VPDmax - VPD)/(VPDmax - VPDmin)
  fvpd <- fvpd %>% rmin(1) %>% rmax(0)
  
  
  # CANOPY CONDUCTANCE, Gc
  # GPP
  # for full descriptions, see supplement of Gan et al., 2018
  if(is_PMLV2){
    
    
    PAR_mol <- PAR * 4.57 # from [W m-2] to [umol m-2 s-1]
    
    
    fT2 <- exp(0.031 * (Tavg - 25)) / (1 + exp(0.115 * (Tavg - 41)))
    fT2 <- fT2 %>% rmin(1)
    
    Am <- Am_25 * fT2
    
    
    P1 <- Am * Alpha * Thelta * PAR_mol
    P2 <- Am * Alpha * PAR_mol
    P3 <- Am * Thelta * Ca
    P4 <- Alpha * Thelta * PAR_mol * Ca 
    
    Ags <- Ca*P1/(P2*kQ + P4*kQ) * (kQ*LAI + log((P2+P3+P4)/(P2+P3*exp(kQ*LAI) + P4))) # umol cm-2 s-1
    
    # 1.0368 = 60.*60./10^6.*24*12;  [umol m-2 s-1] to [g C m-2 d-1]
    GPP <- Ags * 1.0368 * fvpd
    
    # 1.6 = conductance of water / conductance of CO2 (mol m-2 s-1)
    Gc <- m/Ca*Ags*1.6*fvpd_gc   
    
    # unit convert to m s-1
    Gc <- Gc * 1e-2/(0.446 * (273 / (273+Tavg)) *(p/101.3)) 
    
  }else{ 
    
    Gc <- gsx/kQ*log((PAR+Q50)/(PAR*exp(-kQ*LAI)+Q50))*fvpd_gc
    
  }
  
  Gc <- Gc %>% rmax(1e-6)
  
  
  # known bug: bare, ice & snow, unc, all zero parameters will lead to p1, p2, p3, p4 = 0,
  #            GPP = 0/0(masked), and Ec = masked.
  
  # AERODYNAMIC CONDUCTANCE 
  d <- hc * 0.64 # xu.edit 0.667
  zom <- hc * 0.13  # xu.edit 0.123
  zoh <- zom * 0.1
  uz <- log(67.8*Zob - 5.42)/4.87 * u2  #u2 wind speed at 2m height  zob height when meansure
  Ga <- uz * kmar * kmar / (log((Zob - d)/zom) * log((Zob - d)/zoh))
  
  
  # PET (ET_pot), and ET_water
  Eeq <- (slop / (slop + gama) * Rn 
          / lambda * 86.4) %>% # convert W/m2/s into mm
    rmax(0.0001)
  
  Evp <- (gama/(slop + gama))*((6430 * (1 + 0.536*u2) * VPD)/lambda) %>%
    rmax(0)
  
  ET_pot <- Eeq + Evp
  ET_water <- ET_pot
  
  
  # Conductance and ET component
  Tou <- exp(-kA*LAI)
  # Transpiration from plant cause by radiation water transfer
  LEcr <- slop/gama*Rn *(1 - Tou)/(slop/gama + 1 + Ga/Gc)
  # Transpiration from plant cause by aerodynamic water transfer
  LEca <- (rou_a * Cp * Ga * VPD / gama)/(slop/gama + 1 + Ga/Gc)
  
  # making sure vegetation transpiration is negaligable, this is very important for very dry Sahara
  # Should take it seriously. LAI = 0, will lead to a extremely large value. 
  # Update 24 Aug'2017, kongdd
  LAI0_mask <- LAI
  LAI0_mask[LAI0_mask >= 0] <- 1
  LEcr <- LEcr * LAI0_mask
  LEca <- LEca * LAI0_mask
  
  LEc <- LEcr + LEca
  
  # Soil evaporation at equilibrium
  LEs_eq <- (slop/gama)* Rn *Tou/(slop/gama + 1)
  
  
  # cover unit to mm
  coef_MJ2mm <- lambda/86.4 # ET./lambda*86400*10^-3
  Es_eq <- LEs_eq/coef_MJ2mm
  Ecr <- LEcr/coef_MJ2mm
  Eca <- LEca/coef_MJ2mm
  Ec <- LEc/coef_MJ2mm
  
  
  # Interception Precipitation Evaporation: prcp_real = prcp - Ei 
  # @references 
  # Van Dijk, A.I.J.M. and Warren, G., 2010. The Australian water resources assessment system. Version 0.5, 3(5). P39
  fveg <- 1 - exp(-LAI/LAIref)
  Sveg <- S_sls * LAI
  
  fER <- fveg * fER0
  Pwet <- -log(1 - fER0) / fER0 * Sveg / fveg
  
  # evaporation from canopy interception 
  Ei <- (P < Pwet) * fveg * P + (P >= Pwet) * ( fveg*Pwet + fER*(P - Pwet) )
  
  
  # throughfall 
  Pi <- P - Ei
  
  # calculate actual Es from Es_eq and fval_soil
  Es_eq_mean <- movmean(Es_eq, 3)
  Pi_mean <- movmean(Pi, 3)
  # Es_eq_mean <- mean(Es_eq)
  # Pi_mean <- mean(Pi)
  fval_soil <- (Pi_mean / Es_eq_mean) %>% rmax(0) %>% rmin(1)
  Es <- Es_eq * fval_soil
  
  
  # set up mask and calculate ET
  if(land[1] != 0 & land[1] != 15){
    ET_water <- NA
  }
  else{
    Es <- NA
    Ec <- NA
    Ei <- NA
  }
  ET <- Es + Ec + Ei
  
  LEc1 <- LE1 - LEs_eq
  #LEc1 <- Ec1 * coef_MJ2mm
  # <- Ga / (((slop/gama) * Rn * (1-Tou) + (rou_a * Cp * Ga * VPD / gama)) / LEc1 - 1 - slop / gama )
  Gc1 <- (Ga * LEc1) / ((slop/gama) * Rn * (1-Tou) + rou_a * Cp * Ga * VPD / gama - LEc1 - LEc1*(slop/gama))
  if(is_PMLV2){
  output <- data.table(Date = input$Date, site = input$site, IGBPname = input$IGBPname, Prcp = input$Prcp, Tavg = input$Tavg, LAI = input$LAI,
                       ID = input$ID, LE = input$LE, GPP, ET, Es, Ec, Ei, ET_pot, Gc, Gc1, hw <- input$hw) 
  }else{
    output <- data.table(Date = input$Date, site = input$site, IGBPname = input$IGBPname, Prcp = input$Prcp, Tavg = input$Tavg, LAI = input$LAI,
                         ID = input$ID, LE = input$LE, LEc, ET, Es, Ec, Ei, ET_water, ET_pot, GPP, Gc, Gc1, hw <- input$hw)
  }
  
  output
  
}# PML_site



#==============================================================
# FUNCTIONS
#==============================================================
vapor_pressure = function(t) {
  return(0.6108 * exp(17.27 * t / (t + 237.3)))
}

get_lambda <- function(img) {
  # Cp = 4.2 * 0.242;   
  # specific heat at constant pressure, 1.013[kJ kg - 1 0C - 1]
  # lambda = 2500 - 2.2 * Tavg; // J g - 1
  # return lambda;
  return(2500 - 2.2 * img)
}


#' convert W m-2 to mm
#' @param LE LE with the unit of 
#' @param Tavg mean air temperature, in degC
W2mm = function(LE, Tavg) {
  lambda = get_lambda(Tavg)
  return(LE * 86400 * 1e-3 / lambda)
}

#================================================================================
# raster related 
#================================================================================
rmax <- function(r, number){
  r[r < number] <- number
  return(r)
}
rmin <- function(r, number){
  r[r > number] <- number
  return(r)
}

#================================================================================
# time
#================================================================================
dn <- function(Date, n){
  DOY <- yday(Date)
  return(ceiling(DOY/8))
}

# moving mean
movmean <- function(x, win_back){
  start <- sapply(seq_along(x) - win_back, function(x) max(x, 1))
  end <- seq_along(x)
  
  imovmean <- function(iend, start, x){
    val <- mean(x[iend:start[iend]])
  }
  movmean <- sapply(end, start, x, FUN = imovmean)
}

