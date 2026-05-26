# LBiST: Look–Bin–Space–Time Bayesian Retracking for Inland SAR Altimetry

This repository contains the MATLAB implementation of the **LBiST** framework, a Bayesian waveform decontamination and retracking approach for inland water monitoring using **Sentinel-3 SRAL SAR Level-1BS** data.

The method operates on **look-bin maps** before the standard multilooking stage and exploits information in four domains:

- **Look**: Doppler looks in the SAR stack  
- **Bin**: range-delay samples  
- **Space**: along-track neighboring measurements  
- **Time**: repeat-cycle information across overpasses  

By combining spatial likelihoods and temporal priors, the method identifies the **nadir-water echo zone** in the look-bin map, regenerates cleaner waveforms, and improves retracking performance in challenging inland environments such as narrow rivers, fragmented reservoirs, and complex shorelines.

---

## Overview

Conventional SAR altimetry products provide multilooked Level-2 waveforms, where land and water contributions may already be mixed. This repository implements a preprocessing and retracking framework that works directly on **Sentinel-3 Level-1BS delay–Doppler data** to reduce contamination before waveform formation.

Main steps:

1. Load Sentinel-3 Level-1BS look-bin data
2. Align look-bin maps in the bin domain
3. Build spatial and temporal information stacks
4. Apply Bayesian inference to isolate the nadir-water echo zone
5. Regenerate decontaminated waveforms
6. Retrack the regenerated waveforms
7. Compare against in-situ water levels and baseline retrackers

---

## Features

- Processing of **Sentinel-3 SRAL SAR Level-1BS** data
- Construction of **4D Look–Bin–Space–Time stacks**
- PCA-based extraction of spatial echo structure
- FFT-based extraction of temporal/seasonal behavior
- Bayesian fusion of spatial and temporal information
- Regeneration of cleaner waveforms before retracking
- Comparison with standard retrackers such as:
  - OCOG
  - OCEAN / SAMOSA-based ocean range from Level-2 products
  - ICE retracker
  - Iterative Threshold retracker
  - optional SAMOSA+ interface
- Evaluation using:
  - Correlation
  - KGE
  - NSE
  - RMSE
  - Outlier counts
  - cycle-median analysis

