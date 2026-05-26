# LBiST: Look–Bin–Space–Time Bayesian Retracking for Inland SAR Altimetry

This repository contains the MATLAB implementation of the **LBiST** framework, a Bayesian waveform decontamination and retracking approach for inland water monitoring using **Sentinel-3 SRAL SAR Level-1BS** data.

The method operates on **look-bin maps** before the standard multilooking stage and exploits information in four domains:

- **Look**: Doppler looks in the SAR stack  
- **Bin**: range-delay samples  
- **Space**: along-track neighboring measurements  
- **Time**: repeat-cycle information across overpasses  

By combining spatial likelihoods and temporal priors, the method identifies a **Nadir water-only echo zone** in the look-bin map, regenerates cleaner waveforms, and improves retracking performance in challenging inland environments such as narrow rivers, fragmented reservoirs, and complex shorelines.

---

## Overview

Conventional SAR altimetry products provide multilooked Level-2 waveforms, where land and water contributions may already be mixed. This repository implements a preprocessing framework that works directly on **Sentinel-3 Level-1BS delay–Doppler data** to reduce contamination before waveform formation.

Main steps:

1. Load Sentinel-3 Level-1BS look-bin data
2. Align look-bin maps in the bin domain
3. Build spatial and temporal information stacks
4. Apply Bayesian inference to isolate the water-only echo zone
5. Regenerate decontaminated waveforms
6. Retrack the regenerated waveforms
7. Compare against in-situ water levels and baseline retrackers

---

## Features

- Processing of **Sentinel-3 SRAL SAR Level-1BS** data
- Construction of **4D Look–Bin–Space–Time stacks**
- Bayesian filtering of look-bin maps
- Regeneration of cleaner waveforms before retracking
- Comparison with standard retrackers such as:
  - OCOG
  - OCEAN / SAMOSA-based ocean range from Level-2 products
  - ICE retracker
- Evaluation using:
  - Correlation
  - KGE
  - NSE
  - RMSE
  - Outlier counts

---

## Repository Structure

```text
.
├── main_script.m              # Main processing script
├── Lib_Initialization.m       # Initialization of paths / constants
├── Lib_Load_Datasets.m        # Loads case-study datasets and references
├── Lib_Bayesian_PCA_FFT.m         # Spatial-temporal Bayesian / manifold processing
├── Lib_Cycle_Analysis2.m      # Cycle-wise analysis
├── Lib_samosa_retracker.m     # SAMOSA retracker interface / implementation
├── *.mat                      # Intermediate saved data products
└── README.md
