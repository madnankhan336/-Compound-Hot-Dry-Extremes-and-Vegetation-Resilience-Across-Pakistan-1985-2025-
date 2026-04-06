# Compound Hot–Dry Extremes and Vegetation Resilience Across Pakistan (1985–2025)

> An event-based Landsat NDVI anomaly assessment of vegetation resistance, recovery, and elevation-dependent vulnerability across Pakistan.

## Overview

This repository contains the analytical workflow, scripts, figures, tables, and supporting material for a national-scale study of compound temperature–moisture extremes and vegetation response across Pakistan from **1985 to 2025**. The study develops a **month-resolved, event-based framework** using **Landsat-derived NDVI anomalies** and complementary hydroclimatic variables to quantify where vegetation is more resistant to compound hot–dry stress and where recovery is delayed.

The analysis was designed to address three linked questions:

1. How are compound hot–dry and hot–wet extremes distributed across Pakistan?
2. How does vegetation resistance and recovery vary across space and elevation?
3. Which climatic, topographic, and antecedent ecological conditions explain spatial variation in vegetation response?

The main modelling and attribution workflow focuses on **compound hot–dry (HD)** events because these were widespread across Pakistan and supported stable national inference, whereas **hot–wet (HW)** events were comparatively sparse under the adopted month-of-year percentile framework.

## Study area and temporal scope

- **Region:** Pakistan  
- **Period:** January 1985 to December 2025  
- **Temporal resolution:** Monthly  
- **Vegetation indicator:** Landsat NDVI  
- **Climate framework:** Compound heat–moisture extremes  
- **Topographic synthesis:** Five elevation belts  
  - `<500 m`
  - `500–1500 m`
  - `1500–2500 m`
  - `2500–3500 m`
  - `≥3500 m`

## Data sources

### Vegetation
- **Landsat Collection 2, Level 2 Surface Reflectance**
  - Landsat 5 TM
  - Landsat 7 ETM+
  - Landsat 8/9 OLI/OLI-2

Monthly NDVI composites were generated from clear-sky observations using a **P75 compositing strategy**, with conservative short-gap interpolation applied only where justified.

### Climate
- **TerraClimate**
  - Maximum temperature (`tmmx`)
  - Precipitation (`pr`)
  - Vapour pressure deficit (`VPD`)
  - Soil moisture
  - Potential evapotranspiration (`PET`)
  - Actual evapotranspiration (`AET`)
  - Climatic water balance (`pr - PET`)

### Topography
- Elevation was used for both **stratified analysis** and **predictive modelling**.

## Methodological summary

### 1. NDVI anomaly construction
Monthly Landsat NDVI composites were transformed into **calendar-month standardised anomalies (z-scores)** at pixel level. This allowed vegetation response to be compared across ecosystems with very different seasonal behaviour and mean productivity.

### 2. Compound event detection
Compound extremes were identified using **pixel-wise, calendar-month percentile thresholds**, ensuring that “extreme” conditions were defined relative to the local seasonal baseline rather than by a single fixed annual threshold.

The baseline hot–dry definition used:
- **Hot:** maximum temperature > 95th percentile
- **Dry:** precipitation < 10th percentile

Flagged months were converted into **discrete compound events** using run-length segmentation.

### 3. Vegetation response metrics
For each event, two core ecological response metrics were derived:

- **Resistance:** the minimum NDVI anomaly during the event  
- **Recovery:** the number of months required for NDVI anomaly to return to baseline after event termination

### 4. Elevation-aware synthesis
Exposure, resistance, recovery, and amplification were summarised across five elevation zones to assess whether mountain environments differ systematically from lowland and mid-elevation systems.

### 5. Interpretable modelling
Spatial variation in hot–dry resistance was modelled using:
- **Random Forest**
- **Quantile Regression Forest**
- **Grouped cross-validation by pixel ID**
- **Permutation importance**
- **SHAP-based attribution**

This modelling framework was used to explain resistance patterns and generate uncertainty-aware spatial prediction surfaces under typical hot–dry conditions.

### 6. Robustness analysis
Sensitivity tests examined whether the main results were stable under:
- alternative hot thresholds
- stricter dry thresholds
- alternative dryness proxies
- climatic water balance
- combined VPD–soil moisture limitations

## Main findings

The analysis shows that **compound hot–dry stress is the dominant compound climatic hazard for vegetation across Pakistan**, whereas compound hot–wet events are relatively rare under the adopted monthly framework.

Several broad patterns emerged:

- Hot–dry events occurred across most of Pakistan and across the full elevation range.
- Vegetation resistance was generally weakly positive in low to mid elevations, but shifted towards stronger suppression in the highest terrain.
- Recovery was usually rapid in lower elevations but became markedly slower at high elevations.
- Elevation strongly shaped vegetation vulnerability, with the **highest recovery delays observed in the ≥3500 m belt**.
- Antecedent vegetation condition, climatic water balance, temperature, and elevation emerged as the leading controls on hot–dry resistance.
- Sensitivity tests showed that exposure estimates changed with threshold choice, but the broader ecological interpretation remained stable.

## Repository structure

```text
.
├── manuscript/              # Main manuscript files
├── supplementary/           # Supplementary tables, figures, and notes
├── scripts/                 # Data processing, analysis, modelling, and plotting scripts
├── data/
│   ├── raw/                 # Raw input data metadata or download instructions
│   └── processed/           # Derived analysis-ready datasets
├── figures/                 # Main manuscript figures
├── tables/                  # Main manuscript tables
├── results/                 # Model outputs, summaries, and intermediate products
├── docs/                    # Additional project notes and documentation
├── README.md
├── CITATION.cff
└── .gitignore
