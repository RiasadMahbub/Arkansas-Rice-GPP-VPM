# Arkansas Rice GPP: Spatial Estimation and Yield Analysis via VPM (2008–2020)

This repository contains the full analytical workflow for:

> **Mahbub, R. B., Reba, M. L., & Runkle, B. R. (2025).** Magnitude, drivers, and patterns of gross primary productivity of rice in Arkansas using a calibrated vegetation photosynthesis model. *Agricultural and Forest Meteorology, 369*, 110583.

---

## Overview

We estimate daily and cumulative Gross Primary Productivity (GPP) of rice across Arkansas from 2008 to 2020 using the Vegetation Photosynthesis Model (VPM). The model is first calibrated and validated against eddy covariance (EC) tower measurements at 10 rice fields (16 site-seasons, 2015–2018), then applied statewide using MODIS satellite reflectance and NCEP CFS gridded meteorology.

**Site-calibrated parameters for Arkansas rice:**
- ε₀ (LUEmax) = 0.060 mol CO₂ mol⁻¹ PPFD
- T_opt = 30.02 °C

**Key findings:**
- Mean statewide cumulative GPP (2008–2020): 1784 ± 156 g C m⁻² year⁻¹
- Grand Prairie exhibited the highest GPP among six rice production zones (1950 ± 295 g C m⁻² year⁻¹)
- EVI explained 68% of spatial variation in cumulative GPP
- Rice yield increased significantly over time (Mann-Kendall p = 0.04), but statewide GPP trend was not significant — suggesting yield gains were driven by agronomic rather than photosynthetic improvements

---

## Study Sites

Ten EC tower rice fields in Arkansas, covering 2015–2018 (16 site-seasons):

![Site map](Figures/Figure.png)

Sites include US-HRA and US-HRC (central Arkansas, Lonoke County) and US-BDA, US-BDC, and US-OF1 through US-OF6 (northeastern Arkansas, Mississippi County).

---

## Repository Structure

```
GPPPVPM/
│
├── Scripts/
│   └── Workflow/                         # All R scripts, in execution order
│       ├── Function.R
│       ├── SiteScaleAnalysis_DataReadingMerging.R
│       ├── AllDriverinOneScript.R
│       ├── AllDriverinOneScript-Projecting500.R
│       ├── ModeledPVPMVPM_satelliteProcessing.R
│       ├── VPMYieldGPP_CountyMaps_Yearwise.R
│       ├── VPMMeanRasterImageAnalysis2008_2020.R
│       ├── VPM_InterannualGraph.R
│       ├── Drivers_EVI_T_Precipitation_LSWI.R
│       └── Distributioncheck50percentricepixelversusfullricepixel.R
│
├── Figures/
│   ├── Figure.png                        # EC tower site map
│   ├── GPPriceproductionregion.png       # GPP by rice production zone
│   └── VPMcumulativearranged.png         # Statewide cumulative GPP map (2008–2020)
│
└── Data/                                 # Input data (see Data Sources below)
```

---

## Workflow

### Step 1 — Site-Scale Data Assembly and Calibration

```
Scripts/Workflow/Function.R
Scripts/Workflow/SiteScaleAnalysis_DataReadingMerging.R
```

Reads and merges EC tower GPP with site meteorological data (temperature, PAR) and MODIS reflectance (EVI, LSWI) at 8-day intervals. Calibrates ε₀ and T_opt using 78-iteration k-fold cross-validation (11 training / 2 testing / 3 validation site-seasons). The calibrated values (ε₀ = 0.060, T_opt = 30.02 °C) reduced RMSE and MAE by 23% and 25% respectively versus biome-default parameters.

### Step 2 — Driver Raster Preparation

```
Scripts/Workflow/AllDriverinOneScript.R
Scripts/Workflow/AllDriverinOneScript-Projecting500.R
```

Reprojects and harmonizes EVI, LSWI, temperature, and PAR rasters to a common 500 m resolution to match the MODIS spatial grid. Applies a 0.9 bias correction factor to satellite DSWR and converts to PAR using a factor of 2.02 μmol m⁻² s⁻¹ per W m⁻².

### Step 3 — Statewide GPP Estimation (2008–2020)

```
Scripts/Workflow/ModeledPVPMVPM_satelliteProcessing.R
```

Runs VPM across all Arkansas MODIS rice pixels from 2008 to 2020. Two pixel filters are compared:
- All rice pixels (`VPMspatial`)
- Pixels with ≥50% rice land cover (`VPMspatial_50percentpixel`) — used as the primary model

Rice pixel masks are derived from the USDA Cropland Data Layer (CDL) via the AgKit4EE toolkit in GEE.

**Cumulative GPP map (2008–2020 mean):**

![Cumulative GPP](Figures/VPMcumulativearranged.png)

### Step 4 — GPP by Rice Production Zone

```
Scripts/Workflow/VPMMeanRasterImageAnalysis2008_2020.R
```

Computes and compares mean cumulative GPP across six rice production ecological zones: Grand Prairie, White River, West of Crawley's Ridge, North Delta, Middle Delta, and South Delta.

![GPP by production region](Figures/GPPriceproductionregion.png)

Grand Prairie (highest GPP: 1950 g C m⁻² yr⁻¹) and Middle Delta (lowest: 1697 g C m⁻² yr⁻¹) differ significantly in EVI (p = 7.5×10⁻⁶) and LSWI (p = 4.9×10⁻⁶) but not in temperature or PAR, pointing to agronomic and soil drivers of GPP variability.

### Step 5 — GPP–Yield Relationship (County Scale)

```
Scripts/Workflow/VPMYieldGPP_CountyMaps_Yearwise.R
```

Compares annual cumulative GPP averaged to county scale against USDA-reported rice yield for 26 Arkansas counties (2008–2020). County yield data sourced from USDA NASS Quick Stats.

Results:
- County-scale R² = 0.16 between GPP_cum and yield (p < 0.001)
- Relationship improves with larger planted area per county (R² up to 0.47 for counties with 40,000–59,000 ha)
- Site-calibrated parameters + 50% pixel filter: R² = 0.16, RMSE = 58.5 g rice m⁻² season⁻¹
- Biome-default parameters, no pixel filter: R² = 0.0008

### Step 6 — Interannual Trends

```
Scripts/Workflow/VPM_InterannualGraph.R
```

Mann-Kendall trend tests on statewide GPP and yield time series. Yield trend is significant (p = 0.04, slope = 7.67 g rice m⁻² yr⁻¹); GPP trend is not (p = 0.25), suggesting yield improvements are driven by agronomic/genetic factors rather than increased photosynthesis.

### Step 7 — Driver Analysis

```
Scripts/Workflow/Drivers_EVI_T_Precipitation_LSWI.R
```

Linear regression of GPP against EVI, LSWI, temperature, and PAR across space and time. EVI dominates (R² = 0.68 spatial); temperature is secondary (R² = 0.025); PAR and LSWI contribute minimally.

### Step 8 — Pixel Filter Sensitivity

```
Scripts/Workflow/Distributioncheck50percentricepixelversusfullricepixel.R
```

Compares the GPP distribution from all rice pixels versus the ≥50% rice coverage filter to verify that the 50% threshold does not introduce systematic bias.

---

## R Script Reference

| Script | Purpose |
|---|---|
| `Function.R` | Shared helper functions (VPM equations, scalars, fAPAR) |
| `SiteScaleAnalysis_DataReadingMerging.R` | Site data assembly and VPM calibration |
| `AllDriverinOneScript.R` | Load and merge all driver rasters |
| `AllDriverinOneScript-Projecting500.R` | Reproject drivers to 500 m MODIS grid |
| `ModeledPVPMVPM_satelliteProcessing.R` | Statewide VPM GPP estimation and site validation |
| `VPMYieldGPP_CountyMaps_Yearwise.R` | County-scale GPP–yield analysis and maps |
| `VPMMeanRasterImageAnalysis2008_2020.R` | Year-by-year raster processing and zone summaries |
| `VPM_InterannualGraph.R` | Interannual GPP and yield trend plots |
| `Drivers_EVI_T_Precipitation_LSWI.R` | Driver vs. GPP regression analysis |
| `Distributioncheck50percentricepixelversusfullricepixel.R` | Pixel filter sensitivity check |

---

## Key Dependencies

```r
install.packages(c(
  "raster", "terra", "sf", "sp",
  "ggplot2", "ggpubr", "cowplot", "tidyverse",
  "viridis", "trend", "Metrics", "ggspatial",
  "ggrepel", "ggtext", "gridExtra"
))
```

---

## Data Sources

| Dataset | Variable(s) | Source |
|---|---|---|
| Eddy covariance GPP | GPPEC (16 site-seasons) | Leavitt, Reba, Massey et al. |
| MODIS MOD09A1 v6 | EVI, LSWI (500 m, 8-day) | NASA LP DAAC |
| NCEP CFS v2 | T_mean, T_max, DSWR | Saha et al. (2010, 2014) |
| USDA CDL | Rice pixel masks (2008–2020) | USDA NASS |
| USDA NASS Quick Stats | County rice yield | quickstats.nass.usda.gov |
| PRISM | MAT, MAP (site table) | Daly et al. (2008, 2015) |

---

## Citation

Please cite the published paper when using this code or data:

> Mahbub, R. B., Reba, M. L., & Runkle, B. R. (2025). Magnitude, drivers, and patterns of gross primary productivity of rice in Arkansas using a calibrated vegetation photosynthesis model. *Agricultural and Forest Meteorology, 369*, 110583.

---

## Contact

**Riasad Bin Mahbub**
University of Arkansas
rbmahbub@uark.edu
