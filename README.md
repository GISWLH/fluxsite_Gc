# fluxsite_Gc
## variable
* Date date of specific year
* pecp precipitation
* Gc canopy conductance
* hw whether heatwave day
* indensity heatwave indensity
* Tavg daily mean temperature

## Principle

The main inverse calculation is based on these formulas:
$$
G_c=\frac{\lambda E_cG_a}{\varepsilon Q_{A,c}-\lambda E_c(1+\varepsilon)+(\frac{\rho c_p }{\gamma})D_a G_a}
$$

$$
\lambda E_c =\lambda E-\lambda E_s
$$

$$
\lambda E_s=\frac{f \varepsilon Q_{A,s}}{\varepsilon +1}
$$



It does not involve PML_v2

The GPP of v2 only compares the Gc and backcalculation results of PML in the code.

The input latent heat flux is LE1

line132: `LE1 <- input[['LE']]`

Subtracting soil evaporation

line318: `LEc1 <- LE1 - LEs_eq`

Obtain the inverse calculation result Gc1

line321: `Gc1 <- (Ga * LEc1) / ((slop/gama) * Rn * (1-Tou) + rou_a * Cp * Ga * VPD / gama - LEc1 - LEc1*(slop/gama))`

## Input data

![image-20221027151942593](https://imagecollection.oss-cn-beijing.aliyuncs.com/img/image-20221027151942593.png)

co2、LE、Rln、Prcp、Pa、Rs、Tavg、VPD、U2

It's data from the flux station

LAI, Albedo, and Emiss are time series that I extracted from GEE and need to be interpolated

PMLV2_ Parameters_ TERRA_ Rm_ Pro_ Smooth_ GEE is a parameter table corresponding to different land cover types

