# Resolution limit of the eye

This repository contains data and MATLAB scripts for the project: 'Resolution limit of the eye: how many pixels can we see?'

---

## System Requirements

- **Software:** MATLAB R2021b or later
- **Dependencies:** No external toolboxes required
- **OS tested on:** Windows 11, Ubuntu 22.04
- **Hardware:** No special hardware required; runs on standard desktop/laptop

---

## Installation Guide

No installation is required.

1. Clone or download the repository as a `.zip` file.
2. Open MATLAB and add the `matlab` folder (and `utils/` subfolder) to the MATLAB path:
   ```matlab
   addpath(genpath('path_to_repository/matlab'))

---


## Data Files

* `data/resolution_limit_data.csv` - Observer-level results from the main resolution limit experiment. Each row corresponds to a single threshold measurement. Observer identifiers have been anonymized. The description of the variables is as follows:

| Variable               | Description                                                                 |
|--------------------------|-----------------------------------------------------------------------------|
| `observer`               | Anonymized observer ID (e.g., Obs001, Obs002, ...)                          |
| `luminance_cd_m2`        | Stimulus luminance in candela per square meter                              |
| `eccentricity_deg`       | Eccentricity of the stimulus in visual degrees                              |
| `color_direction`        | Categorical variable indicating chromatic direction (1: achromatic, 2: red-green, 3: yellow-violet)                       |
| `stimulus_size_deg`      | Size of the stimulus in degrees of visual angle                             |
| `stimulus_shape`         | Shape type of the stimulus encoded as categorical variable (1: square-wave grating modulated with a Gaussian envelope)                |
| `num_trials`             | Number of trials included in the threshold estimation                       |
| `threshold`              | Estimated resolution threshold in units of cube-root spatial frequency                                 |
| `threshold_ci_low`       | Lower bound of the confidence interval for the threshold estimate           |
| `threshold_ci_high`      | Upper bound of the confidence interval for the threshold estimate           |
| `threshold_ppd`          | Resolution threshold in pixels per degree (ppd)                               |
| `threshold_ppd_ci_low`   | Lower bound of ppd-based threshold confidence interval                      |
| `threshold_ppd_ci_high`  | Upper bound of ppd-based threshold confidence interval                      |
| `contrast_sensitivity`   | Inverse of the grating contrast used (cone contrast units)                  |

* `data/resolution_limit_median_data.csv` - Summarised median threshold results across participants, used to generate predicted heatmaps.

* `data/text_resolution_limit_data.csv` - Resolution limit data from the text stimuli experiment. The description of the variables is the same as in the main experiment data. 

## Code Files

MATLAB scripts to generate the figures in the paper are provided in the directory `matlab`.
* `plot_res_limit.m` - Generates **Figure 1a** (Figure 2a in Nature Communications paper): Violin plots showing resolution limit (in pixels per degree) across eccentricities and chromatic directions from the data in `data/resolution_limit_data.csv`.
* `plot_heatmap_prediction_interval.m` - Generates **Figure 1b** (Figure 2b in Nature Communications paper): A heatmap representation of the modeled spatial resolution limits across eccentricity. Confidence intervals and predicted contours are drawn using model fits based on Watson (2018).









