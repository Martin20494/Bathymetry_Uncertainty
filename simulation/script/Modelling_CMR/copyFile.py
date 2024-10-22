# Prepare packages ----------------------------------------------------------------------------------------------------
from folder import *

import geopandas as gpd

import shutil
import pathlib
import glob

# ---------------------------------------------------------------------------------------------------------------------


# WRITING LAND AND BATHY FROM ZIP TO GEOJSON --------------------------------------------------------------------------
def land_and_bathy(
    path_in,
    path_out
):
    """
    @Definition:
                A function to write out land and bathymetry geojson file from zip file downloaded from LINZ
    @References:
                None.
    @Arguments:
                path_in (string):
                                Path of original land/bathy files
                path_out (string)
                                Path of writing out land/bathy geojson files
    @Returns:
                None.
    """
    # Read file in
    file = gpd.read_file(path_in)

    # Add crs
    geom_crs = file.to_crs(epsg=2193)

    # Write out geojson file
    geom_crs.to_file(path_out, crs=2193, driver="GeoJSON")

# END WRITING LAND AND BATHY FROM ZIP TO GEOJSON ----------------------------------------------------------------------

# COPYING ALL GEOMETRY FILES ------------------------------------------------------------------------------------------
def copy_geom(
    layers,
    number_simulation,
):
    """
    @Definition:
                A function to copy all geometries
    @References:
                None.
    @Arguments:
                layers (list):
                            List of id/number of land/bathy layers
                number_simulation (int):
                            Id of simulation
    @Returns:

    """
    # Bathy files
    # Create folder for writing out necessary files for DEM
    bathyfiles_folder_out = fr"{bathy_path}\\bathy_{number_simulation}"
    pathlib.Path(fr"{bathyfiles_folder_out}").mkdir(parents=True, exist_ok=True)
    # Create folder for writing out necessary files for roughness
    bathyfiles_roughness_folder_out = fr"{bathy_path}\\bathy_roughness_{number_simulation}"
    pathlib.Path(fr"{bathyfiles_roughness_folder_out}").mkdir(parents=True, exist_ok=True)


    # Land and bathy lines
    for layer_number in layers:
        # Get original file path
        zip_path = fr"{original_lidar_path}\\{layer_number}\\*.zip"

        # Get files from zip file
        landbathy_path_in = glob.glob(zip_path)[0]

        # Generate path out for land and bathy lines
        # Get names
        if layer_number == 51559:
            landbathy_path_out = fr"{bathyfiles_folder_out}\\land.geojson"
            landbathy_roughness_path_out = fr"{bathyfiles_roughness_folder_out}\\land.geojson"

        else:
            landbathy_path_out = fr"{bathyfiles_folder_out}\\bathymetry_contours.geojson"
            landbathy_roughness_path_out = fr"{bathyfiles_roughness_folder_out}\\bathymetry_contours.geojson"

        # Copy land and bathy lines
        land_and_bathy(landbathy_path_in, landbathy_path_out)
        land_and_bathy(landbathy_path_in, landbathy_roughness_path_out)

    # Get directory path of reference folder (folder stores original files)
    bathyfiles_ref_path = pathlib.Path(fr"{original_lidar_path}\\no_padding")

    # Create path for all files but exclude catchment_boundary
    filtered_bathyfiles_ref_path = [x for x in bathyfiles_ref_path.glob("**\*.geojson") if not x.name.startswith(("river_bathymetry", "catchment_boundary"))]

    # Copy all bathy files to a specific simulation
    for num in range(len(filtered_bathyfiles_ref_path)):
        bathyfiles_path_in = filtered_bathyfiles_ref_path[num]
        bathyfiles_path_out = fr"{bathyfiles_folder_out}\\{pathlib.Path(bathyfiles_path_in).stem}.geojson"
        shutil.copy2(bathyfiles_path_in, bathyfiles_path_out)
# END COPYING ALL GEOMETRY FILES ------------------------------------------------------------------------------------
