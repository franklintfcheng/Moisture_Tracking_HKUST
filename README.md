# Lagrangian moisture tracking diagnostics

**Code Developer:** Dr. Franklin Tat Fan Cheng (franklin.cheng@ust.hk)

**Affiliation:** Hong Kong University of Science and Technology (HKUST)

---

**Description**  
This repository consists of the source codes for the backward **`WaterSip`** and forward **`WaterDrip`** diagnostics based on the particle output position (`partpositXXXX`) from [FLEXPART](https://www.flexpart.eu/) simulations. 

**Usage**  
- To perform moisture tracking, go to `moisture_tracking.m`. Adapt paths/parameters to your needs.

- To plot the backtracked moisture source, go to `flexpart_tutorial.m`. Adapt paths/parameters to your needs.
  
- The programs were all written in Matlab (ver 2023a).


**Updates** 

*-- July 5, 2024 --*
1. Added AR-related moisture contribution in WaterSip diagnostic (`Pm_AR`)
2. Minimized the hard-coded portions from `moisture_tracking.m` and `flexpart_tutorial.m`


**References** 

- Cheng, T. F., and M. Lu, 2023: Global Lagrangian Tracking of Continental Precipitation Recycling, Footprints, and Cascades. J. Climate, 36, 1923–1941, https://doi.org/10.1175/JCLI-D-22-0185.1.

- Sodemann, H., C. Schwierz, and H. Wernli (2008), Interannual variability of Greenland winter precipitation sources: Lagrangian moisture diagnostic and North Atlantic Oscillation influence, J. Geophys. Res., 113, D03107, https://doi.org/10.1029/2007JD008503.
