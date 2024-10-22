Download `data` and `scripts` folders into your local computer. 

#### 1. File data

The `data` folder: Includes `0_lidar_data` and `other_data` sub-folders. These files should be put into a simulation version folder (e.g: Put them into folder version_001 used to generate simulations of rotation and translation).

- `0_lidar_data`: Includes `flow_and_friction.csv.gz` and `rec2_3.geojson` files (river features) provided by NIWA to estimate the river bathymetry. These data can be found [here](https://data-niwa.opendata.arcgis.com/apps/NIWA::new-zealand-river-flood-statistics-app/explore).

- `other_data`: includes shape file of observed flood debris (January-2005 flood event at Waikanae River) and linestring represents for tide boundary used for flood modelling.

#### 2. File scripts

The `scripts` folder includes `Analysing`, `Modelling`, and `Executing` sub-folders. These sub-folders are used for generating simulations. 

- `Modelling` and `Analysing`: Includes two separate groups of modules to generate and then analyse the simulations. `Modelling` will create multiple DEMs from LiDAR data downloaded from the [OpenTopography](https://portal.opentopography.org/datasets). These DEMs will then be used in LISFLOOD-FP model to produce multiple water depths and water surfaces. These data will then be analysed in `Analysing`.

- `Executing`: Includes `Analysing`, `Comparing`, and `Modelling` files used to run modules in `Analysing` and `Modelling` folders to generate and compare simulations. These files were written in jupyter notebook format. Therefore, we suggest to also use jupyter notebook with the environment installed earlier. Please follow the instructions for generating each simulation version for different transformation types, resolutions, and flood return periods written in the script. Please change the `os.chdir` in each file into the folder containing `Analysing` or `Modelling` folders.

#### 3. Usage
