Download `data` and `scripts` folders into your local computer. 

#### 1. File data

The `data` folder: Includes `CMR` and `UF` folders for the Conceptual Multivariate Regression and Uniform Flow formulas. Each formula folder includes `0_lidar_data` and `other_data` sub-folders. These files should be put into a simulation version folder (e.g: Put them into folder `cmr_slope_001` used to generate simulations of rotation and translation).

- `0_lidar_data`: Includes `flow_and_friction.csv.gz` and `rec2_3.geojson` files (river features) provided by NIWA to estimate the river bathymetry. These data can be found [here](https://data-niwa.opendata.arcgis.com/apps/NIWA::new-zealand-river-flood-statistics-app/explore).

- `other_data`: Includes shape file of observed flood debris (January-2005 flood event at Waikanae River) and linestring represents for tide boundary used for flood modelling.

- `cache`: Includes cache files from OpenStreetMap.

#### 2. File scripts

The `scripts` folder includes `Analysing`, `Modelling_CMR`, `Modelling_UF`, and `Executing` sub-folders. These sub-folders are used for generating simulations. 

- `Modelling_CMR`, `Modelling_UF`, and `Analysing`: Includes three separate groups of modules to generate and then analyse the simulations. `Modelling_CMR` and `Modelling_UF` will create multiple DEMs from LiDAR data downloaded from the [OpenTopography](https://portal.opentopography.org/datasets) with the river depths estimated by the Conceptual Multivariate Regression and Uniform Flow formulas. These DEMs will then be used in LISFLOOD-FP model to produce multiple water depths and water surfaces. These data will then be analysed in `Analysing`.

- `Executing`: Includes `Analysing`, `Comparing`, and `Modelling` files used to run modules in `Analysing` and `Modelling` folders to generate and compare simulations. These files were written in jupyter notebook format. Therefore, we suggest to also use jupyter notebook with the environment installed earlier. Please follow the instructions for generating each simulation version for different transformation types, resolutions, and flood return periods written in the script. Please change the `os.chdir` in each file into the folder containing `Analysing` or `Modelling` folders.

#### 3. Usage

After [installation](https://github.com/Martin20494/Bathy_Uncertainty?tab=readme-ov-file#environment-installation), download these subfolders to a local folder. In `Modelling_CMR`, `Modelling_UF`, `Analysing` (using the `forSimulation` to analyse the results, the `forRMSE` is for calculating RMSE distributions only), change the `MAIN_DIR` (path to your local folder) in module `folder.py`, change the name of simulation version in module `versionModule`. For example, I named `rupp_slope_001` for the river depths estimated by the Conceptual Multivariate Regression formula when the errors were added to the river slope parameter. You might also want to change [API key](https://www.linz.govt.nz/guidance/data-service/linz-data-service-guide/web-services/creating-api-key) to crawl data from LINZ in `dataPreparation.py` in `Modelling_CMR` and `Modelling_UF`.

There are 8 simulation versions, each of them represents fo a set of 50 simulations. These simulations are based on the river depths estimated by the Conceptual Multivariate Regression and Uniform Flow formulas. In these formulas, the errors were added into the river slope, bank-full flow, and width. To create 50 simulations for each version, open `Executing` folder and access to jupyter notebook `Modelling.ipynb`, change things as following bullet points before running all cells:

- Change the path in `0. Change directory` for cell directing to where the `Modelling.ipynb` is stored in your local computer. (Notice if you would like to use the Conceptual Multivariate Regression or Uniform Flow formulas)
- 
- Change the `2. Necessary variables` to the simulation version you want to run.

- Change the `2.2. Other basic variables` (probably the second cell) and `2.4. Preparing some flood model inputs` to the size and chunk to suit your machine.

- Change the `2.5. Random transformations` to the simulation version you want to run.

- Change the number of processors in all cells of `3. Execution` to suit your machine.

To generate analysis result, in `Executing` folder, access to jupyter notebook `Analysis.ipynb`, change things as following bullet points before running all cells:

- Change the path in `0. Change directory` for cell directing to where the `Analysis.ipynb` is stored in your local computer. (Notice if you would like to use the Conceptual Multivariate Regression or Uniform Flow formulas)

- Change the `2. Data preparation` (probably the resolution) to the simulation version you want to run.

- Change the number of processors in all cells of `3. Generate csv files` to suit your machine.

After generating all (8) subsets, to compare, in `Executing` folder, access to jupyter notebook `Comparison.ipynb`, change things as following bullet points before running all cells:

- Change the parth in the very first cell to let it direct to where the `Analysis.ipynb` is stored in your local computer. 

- Change the paths in `plot_cv`, `plot_areas`, and `plot_RMSE` into your local paths where you store all the subsets. (Notice that for the RMSE, you might want to see my versions as I have another folder in each version folder to include another simulation of when no errors were added into the parameters)

- Change the paths to store all the results, they are in cells with the last line of codes with comment line `# Save fig`.

**_Notice_**: This project was run on CPU machine with 8 cores and 64GB memory, you might want to choose different `size_of_processor` and `size_of_chunk` in `Modelling` files according to your machine features. It will take 1-2 days to generate a version but some might take more than a week (e.g. 2-meter resolution). Please contact with the Github author if further information is needed. 

#### 4. Result files explaining

The image below shows how files look like after running these modules.
  
<div align="center">
	<img width = "90%" src="https://github.com/Martin20494/Bathymetry_Uncertainty/blob/main/simulation/folder_example/folder_examples.jpg">
</div>

There are 8 folders representing for 8 simulation versions using the river depths estimated by the Conceptual Multivariate Regression and Uniform Flow formulas with errors in parameters - river slope, bank-full flow, and width - individually and collectively.

- `rupp_slope_001`: river depths estimated by the Conceptual Multivariate Regression formula with errors in the river slope.

- `rupp_flow_001`: river depths estimated by the Conceptual Multivariate Regression formula with errors in the river bank-full flow.

- `rupp_width_001`: river depths estimated by the Conceptual Multivariate Regression formula with errors in the river width.

- `rupp_combination_001`: river depths estimated by the Conceptual Multivariate Regression formula with errors in three parameters collectively.

- `neal_slope_001`: river depths estimated by the Uniform Flow formula with errors in the river slope.

- `neal_flow_001`: river depths estimated by the Uniform Flow formula with errors in the river bank-full flow.

- `neal_width_001`: river depths estimated by the Uniform Flow formula with errors in the river width.

- `neal_combination_001`: river depths estimated by the Uniform Flow formula with errors in three parameters collectively.

Each version includes 6 folders as explained below:

- `0_lidar_data`: Storing original LiDAR data and DEM

- `1_bathy`: Storing simulated river depths

- `2_raster`: Storing transformed DEMs and roughness length as well as Manning's n

- `3_LISFLOODD_FP`: Storing inputs and outputs of flood modelling using LISFLOOD-FP

- `4_untransformation`: Storing reversed outputs

- `5_analysis`: Storing analysis results

- `6_analysis_rmse`: Also storing analysis results but used for RMSE

- `cache`: Storing history of used data

- `other_data`: Other necessary data to generate and evaluate the simulations
