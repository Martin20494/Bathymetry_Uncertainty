# Prepare packages ---------------------------------------------------------------------------
import pathlib                          # For creating and manipulating the directory path
import os                               # For manipulating the directory and to execute commands
from versionModule import version       # For changing version
# --------------------------------------------------------------------------------------------

# Set up PATH --------------------------------------------------------------------------------
# Assign the path to the variable
MAIN_DIR = fr"S:\\Bathymetry\\versions\\{version}"

# Create header path
pathlib.Path(f"{MAIN_DIR}").mkdir(parents=True, exist_ok=True)

# Change the directory
os.chdir(f"{MAIN_DIR}")
# --------------------------------------------------------------------------------------------


# Create paths ###############################################################################
# --------------------------------- This is for TRANSFORMATION PART --------------------------

## General data: Folder stores all files that need to download
other_data = f"{MAIN_DIR}\\other_data"

## 0_lidar_data: Folder stores downloaded lidar point cloud data from Open Topography website
# Reference: https://portal.opentopography.org/datasets
original_lidar_path = f"{MAIN_DIR}\\0_lidar_data"

## 1_bathy: Folder stores paras and bathy
bathy_path = f"{MAIN_DIR}\\1_bathy"

## 2_raster: Folder stores DEM rasters
# DEM -------------
# Netcdf
dem_nc_path = f"{MAIN_DIR}\\2_raster\\dem\\netcdf"
# GeoTiff
dem_tiff_path = f"{MAIN_DIR}\\2_raster\\dem\\tiff"
# Ascii
dem_asc_path = f"{MAIN_DIR}\\2_raster\\dem\\ascii"

# MANNING'S N -------
# Roughness
roughness_nc_path = f"{MAIN_DIR}\\2_raster\\n\\roughness"
# Netcdf
n_nc_path = f"{MAIN_DIR}\\2_raster\\n\\netcdf"
# GeoTiff
n_tiff_path = f"{MAIN_DIR}\\2_raster\\n\\tiff"
# Ascii
n_asc_path = f"{MAIN_DIR}\\2_raster\\n\\ascii"

# STARTDEPTH ------
# Netcdf
startdepth_nc_path = f"{MAIN_DIR}\\2_raster\\startdepth\\netcdf"
# GeoTiff
startdepth_tiff_path = f"{MAIN_DIR}\\2_raster\\startdepth\\tiff"
# Ascii
startdepth_asc_path = f"{MAIN_DIR}\\2_raster\\startdepth\\ascii"

# SHAPEFILE of SEA
bathy_shapefile = f"{MAIN_DIR}\\2_raster\\shapefile"

## 3_LISFLOOD: Folder stores DEM rasters
# Parameters
FPpara_path = f"{MAIN_DIR}\\3_LISFLOOD_FP\\bathy_para"

# Outputs
FPoutput_path = f"{MAIN_DIR}\\3_LISFLOOD_FP\\bathy_output"

## 4_max_depth
extracted_flowdepth = f"{MAIN_DIR}\\4_extraction\\max_depth"
bathy_flowdepth = f"{MAIN_DIR}\\4_extraction\\poly_flowdepth"


## 5_results
other_bathy = f"{MAIN_DIR}\\5_analysis\\bathy_other"

# Variation
csv_bathy = f"{MAIN_DIR}\\5_analysis\\bathy_csv"
raster_bathy = f"{MAIN_DIR}\\5_analysis\\bathy_raster"
plot_bathy = f"{MAIN_DIR}\\5_analysis\\bathy_plot"

# Impact
onepolygon_bathy = f"{MAIN_DIR}\\5_analysis\\bathy_impact\\onepolygon_nobackground"
oneraster_bathy = f"{MAIN_DIR}\\5_analysis\\bathy_impact\\oneraster_nobackground"
# ---------------------------------------------------------------------------------------------------

# Ending path creation ###############################################################################


# Create a list of all sub-folders' paths ------------------------------------------------------------
necessary_files = [
    # Data needs downloading
    other_data,

    # 0_lidar_data
    original_lidar_path,

    # 1_bathy
    bathy_path,

    # 2_dem_raster
    # DEM
    dem_nc_path,
    dem_tiff_path,
    dem_asc_path,
    # MANNING'S N
    roughness_nc_path,
    n_nc_path,
    n_tiff_path,
    n_asc_path,
    # STARTDEPTH
    startdepth_nc_path,
    startdepth_tiff_path,
    startdepth_asc_path,

    # 3_LISFLOOD
    # Parameters
    FPpara_path,
    # Outputs
    FPoutput_path,

    # 4_max_depth
    extracted_flowdepth,
    bathy_flowdepth,

    # 5_results
    other_bathy,
    # Variation
    csv_bathy,
    raster_bathy,
    plot_bathy,
    # Impact
    onepolygon_bathy,
    oneraster_bathy
]
# --------------------------------------------------------------------

# Execute all sub-folders --------------------------------------------
for each_folder in necessary_files:
    pathlib.Path(f"{each_folder}").mkdir(parents=True, exist_ok=True)
# --------------------------------------------------------------------