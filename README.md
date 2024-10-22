# Bathymetry Estimation Uncertainty

### Title
Quantifying uncertainty in flood predictions due to river bathymetry estimation

### Authors
Martin Nguyen, Prof. Matthew D. Wilson, Dr. Emily M. Lane, Prof. James Brasington, and Dr. Rose A. Pearson

## Overview
River bathymetry is important for accurate flood modelling but often unavailable due to the time-intensive and expensive nature of its acquisition. This leads to several proposed and implemented approaches for its estimation. However, the errors in measurements and estimations inherent in these methods, affecting the accuracy of the flood modelling outputs, are not extensively reserached. Hence, our research quantified this uncertainty using a Monte Carlo approach and flood model LISFLOOD-FP to generate multiple flood simulations for analysis. Validation of these simulations against observed flood levels reveals that inaccuracies in river bathymetry estimations can lead to significant deviations in flood predictions, underscoring the need for precise bathymetric data for reliable modeling.

<div align="center">
	<img width = "60%" src="https://github.com/Martin20494/Bathymetry_Uncertainty/blob/main/other_files/all_results/boxplots/S3_boxplot_RMSEs.jpg">
</div>

## Reproducibility

The project was written under modules format and has not been updated to class format yet. Therefore, the users might need to download all the modules to reproduce this project. After the instruction here, please click on each folder for further guidelines. If further information is needed please contact with the Github author.

### Environment installation

Download file ```bathymetry_uncertainty_packages.yaml``` and use the following command (reference from the answer of [@merv](https://stackoverflow.com/questions/76800978/conda-invalidversionspec-invalid-version-error-when-tryin-to-install-from-requi)) to recreate the anaconda environment (only support Windows at the moment).

```
conda env create -n bathymetryuncertainty -f bathymetryuncertainty.yaml
```

### Data

This project includes 2 main folders as described below.

- [simulation](https://github.com/Martin20494/Bathymetry_Uncertainty/tree/main/simulation): Including `Modelling_CMR` and `Modelling_UF` sub-folders used to generate simulations; `Analysis` sub-folder used to generate analysis results; `Executing` sub-folder includes all Jupyter Notebooks used to run the other subfolders to generate and analyse simulations.

- [other_files](https://github.com/Martin20494/Bathymetry_Uncertainty/tree/main/other_files): Including `all_results`, `flow_tide`, `map`, `simulation_process` sub-folders that store all figures, tables, and results in the paper and the research.

## Acknowledgement

This project is part of the NIWA-led national programme, ["Reducing flood inundation hazard and risk across Aotearoa - New Zealand"](https://niwa.co.nz/hazards/ma-te-haumaru-o-nga-puna-wai-o-rakaihautu-ka-ora-mo-ake-tonu-increasing-flood), funded by the Ministry for Business, Innovation and Employment (MBIE) Endeavour Programme. A PhD scholarship was provided from this programme to Martin Nguyen. The LISFLOOD-FP flood model was developed at the University of Bristol and is available [here](https://www.seamlesswave.com/LISFLOOD8.0).

## Contacts

Github author: Martin Nguyen martinnguyen20494@gmail.com
