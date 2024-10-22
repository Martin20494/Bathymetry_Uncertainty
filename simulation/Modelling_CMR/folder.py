# Prepare packages --------------------------------------------------------------------------------
import pathlib                              # For creating and manipulating the directory path
import os                                   # For manipulating the directory and to execute commands
from versionModule import version          # For changing version
# -------------------------------------------------------------------------------------------


# Set up PATH -------------------------------------------------------------------------------
# Assign the path to the variable
MAIN_DIR = f"S:\\Bathymetry\\version013\\{version}"

# Create header path
pathlib.Path(f"{MAIN_DIR}").mkdir(parents=True, exist_ok=True)

# Change the directory
os.chdir(f"{MAIN_DIR}")
# -------------------------------------------------------------------------------------------


# Create paths #############################################################################
# ------------------------------- This is for TRANSFORMATION PART --------------------------

## General data: Folder stores all files that need to download
other_data = f"{MAIN_DIR}\\other_data"

## 0_lidar_data: Folder stores downloaded lidar point cloud data from Open Topography website
# Reference: https://portal.opentopography.org/datasets
original_lidar_path = f"{MAIN_DIR}\\0_lidar_data"

## 1_transformation: Folder stores paras and bathy
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

# SHAPEFILE of PADDING
bathy_shapefile = f"{MAIN_DIR}\\2_raster\\shapefile"


## 3_LISFLOOD: Folder stores DEM rasters
# Parameters
FPpara_path = f"{MAIN_DIR}\\3_LISFLOOD_FP\\bathy_para"

# Outputs
FPoutput_path = f"{MAIN_DIR}\\3_LISFLOOD_FP\\bathy_output"


## 4_max_transformation
# Water depth (wd)
extracted_wd = f"{MAIN_DIR}\\4_extraction\\wd"
bathy_wd= f"{MAIN_DIR}\\4_extraction\\bathy_wd"
# -------
filter_extracted_wd_tiff = f"{MAIN_DIR}\\4_extraction\\filter_wd\\tiff"
filter_extracted_wd_shp = f"{MAIN_DIR}\\4_extraction\\filter_wd\\shp"
filter_extracted_wd_nc = f"{MAIN_DIR}\\4_extraction\\filter_wd\\nc"

# Water surface elevation (wse)
extracted_wse = f"{MAIN_DIR}\\4_extraction\\wse"
bathy_wse= f"{MAIN_DIR}\\4_extraction\\bathy_wse"
# -------
filter_extracted_wse_tiff = f"{MAIN_DIR}\\4_extraction\\filter_wse\\tiff"
filter_extracted_wse_shp = f"{MAIN_DIR}\\4_extraction\\filter_wse\\shp"
filter_extracted_wse_nc = f"{MAIN_DIR}\\4_extraction\\filter_wse\\nc"

# NEW water surface elevation (new_wse)
extracted_new_wse = f"{MAIN_DIR}\\4_extraction\\new_wse"
bathy_new_wse= f"{MAIN_DIR}\\4_extraction\\bathy_new_wse"
# -------
filter_extracted_new_wse_tiff = f"{MAIN_DIR}\\4_extraction\\filter_new_wse\\tiff"
filter_extracted_new_wse_shp = f"{MAIN_DIR}\\4_extraction\\filter_new_wse\\shp"
filter_extracted_new_wse_nc = f"{MAIN_DIR}\\4_extraction\\filter_new_wse\\nc"

# Elevation
bathy_elev = f"{MAIN_DIR}\\4_extraction\\bathy_elev"
# -------
filter_extracted_elev_tiff = f"{MAIN_DIR}\\4_extraction\\filter_elev\\tiff"
filter_extracted_elev_shp = f"{MAIN_DIR}\\4_extraction\\filter_elev\\shp"
filter_extracted_elev_nc = f"{MAIN_DIR}\\4_extraction\\filter_elev\\nc"

# Manning's n
bathy_n = f"{MAIN_DIR}\\4_extraction\\bathy_n"
# -------
filter_extracted_n_tiff = f"{MAIN_DIR}\\4_extraction\\filter_n\\tiff"
filter_extracted_n_shp = f"{MAIN_DIR}\\4_extraction\\filter_n\\shp"
filter_extracted_n_nc = f"{MAIN_DIR}\\4_extraction\\filter_n\\nc"

# ----------------------------------------- This is for ANALYSIS PART -----------------------
## 5_results
# Variation
# Water depth
wd_csv_bathy = f"{MAIN_DIR}\\5_analysis\\wd\\bathy_csv"
wd_raster_bathy = f"{MAIN_DIR}\\5_analysis\\wd\\bathy_raster"
wd_plot_bathy = f"{MAIN_DIR}\\5_analysis\\wd\\bathy_plot"
wd_onepolygon_bathy = f"{MAIN_DIR}\\5_analysis\\wd\\bathy_impact\\onepolygon_nobackground"
wd_oneraster_bathy = f"{MAIN_DIR}\\5_analysis\\wd\\bathy_impact\\oneraster_nobackground"

# Water surface elevation
wse_csv_bathy = f"{MAIN_DIR}\\5_analysis\\wse\\bathy_csv"
wse_raster_bathy = f"{MAIN_DIR}\\5_analysis\\wse\\bathy_raster"
wse_plot_bathy = f"{MAIN_DIR}\\5_analysis\\wse\\bathy_plot"
wse_onepolygon_bathy = f"{MAIN_DIR}\\5_analysis\\wse\\bathy_impact\\onepolygon_nobackground"
wse_oneraster_bathy = f"{MAIN_DIR}\\5_analysis\\wse\\bathy_impact\\oneraster_nobackground"

# Elevation
elev_csv_bathy = f"{MAIN_DIR}\\5_analysis\\elev\\bathy_csv"

# Manning's n
n_csv_bathy = f"{MAIN_DIR}\\5_analysis\\n\\bathy_csv"


# -------------------------------------------------------------------------------------------

# Ending path creation ######################################################################


# Create a list of all sub-folders' paths --------------------------------------------------
necessary_files = [
    # Data needs downloading
    other_data,

    # 0_lidar_data
    original_lidar_path,

    # 1_transformation
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
    # SHAPEFILE
    bathy_shapefile,

    # 3_LISFLOOD
    FPpara_path,
    FPoutput_path,

    # 4_untransformation
    # Water depth
    extracted_wd,
    extracted_wse,
    extracted_new_wse,
    # ---------
    filter_extracted_wd_tiff,
    filter_extracted_wd_shp,
    filter_extracted_wd_nc,
    # Water surface elevation
    bathy_wd,
    bathy_wse,
    bathy_new_wse,
    # ---------
    filter_extracted_wse_tiff,
    filter_extracted_wse_shp,
    filter_extracted_wse_nc,
    # New water surface elevation
    filter_extracted_new_wse_tiff,
    filter_extracted_new_wse_shp,
    filter_extracted_new_wse_nc,
    # Elevation
    bathy_elev,
    # ---------
    filter_extracted_elev_tiff,
    filter_extracted_elev_shp,
    filter_extracted_elev_nc,
    # Manning's n
    bathy_n,
    # ---------
    filter_extracted_n_tiff,
    filter_extracted_n_shp,
    filter_extracted_n_nc,

    # 5_results
    # Variation
    # Water depth
    wd_csv_bathy,
    wd_raster_bathy,
    wd_plot_bathy,
    # Water surface elevation
    wse_csv_bathy,
    wse_raster_bathy,
    wse_plot_bathy,
    # Elevation
    elev_csv_bathy,
    # Manning's n
    n_csv_bathy,

    # Impact
    # Water depth
    wd_onepolygon_bathy,
    wd_oneraster_bathy,
    # Water surface elevation
    wse_onepolygon_bathy,
    wse_oneraster_bathy
]
# -------------------------------------------------------------------------------------------


# Execute all sub-folders -------------------------------------------------------------------
for each_folder in necessary_files:
    pathlib.Path(f"{each_folder}").mkdir(parents=True, exist_ok=True)
# -------------------------------------------------------------------------------------------