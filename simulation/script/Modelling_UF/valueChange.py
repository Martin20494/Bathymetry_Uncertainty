# Prepare packages -----------------------------------------------------------------------------------------------------
from folder import *

# Packages for spatial data manipulation
from shapely.geometry import Polygon                        # For creating polygons to change raster values
import geopandas as gpd                                     # For manipulating shape files
import rioxarray as rxr                                     # For manipulating spatial data under xarray array format

# Packages for data manipulation
import numpy as np
import pandas as pd
# ----------------------------------------------------------------------------------------------------------------------


def value_change(shapefile_func, file_need_changing_func, value_func, inside=True):
    """
    @Definition:
                A function to change pixel values inside or outside polygons
    @References:
                https://corteva.github.io/rioxarray/html/rioxarray.html
                https://corteva.github.io/rioxarray/stable/examples/convert_to_raster.html
                https://gis.stackexchange.com/questions/414194/changing-raster-pixel-values-outside-of-polygon-box-using-rioxarray
                https://automating-gis-processes.github.io/CSC/notebooks/L2/geopandas-basics.html
                https://corteva.github.io/rioxarray/stable/getting_started/nodata_management.html
                https://gdal.org/programs/gdal_translate.html
                https://gis.stackexchange.com/questions/390438/replacing-nodata-values-by-a-constant-using-gdal-cli-tools
    @Arguments:
                shapefile_func (polygon):
                                            Polygon boundaries
                file_need_changing_func (string):
                                            File contains values that need changing
                value_func (int or float):
                                            Values used to replace
                inside (boolean):
                                            If True, change values inside, else, change values outside
    @Returns:
                None.
    """
    # Set up value changing command
    if inside:
        # Change values inside polygons
        inside_command = fr"gdal_rasterize -burn {value_func} {shapefile_func} {file_need_changing_func}"
        os.system(inside_command)
    else:
        # Change values outside polygons
        outside_command = fr"gdal_rasterize -i -burn {value_func} {shapefile_func} {file_need_changing_func}"
        os.system(outside_command)


def polygon_generation(
        number_simulation
):

    """
    @Definition:
                A function to transform the polygon and write out shapefile
    @References:
                None.
    @Arguments:
                number_simulation (string):
                            A string to identify the order of simulation (should be angle, x, y)
                angle_func, x_translation_func, y_translation_func (float):
                            Values to rotate and translate
                center_x_func (float):
                            X coordinate of original/reference DEM
                center_y_func (float):
                            Y coordinate of original/reference DEM
    @Return:
                None.
    """

    # Get extent information of original no-padding raster
    raster_no_padding = rxr.open_rasterio(fr"{original_lidar_path}\\no_padding\\no_padding.nc")
    no_padding_xmin, no_padding_ymin, no_padding_xmax, no_padding_ymax = raster_no_padding.rio.bounds()

    # Create polygon to remove SEA values
    # Set up polygon coordinates into a list
    polygon_padding = [
        (no_padding_xmin, no_padding_ymax),
        (no_padding_xmin + 4115, no_padding_ymax), # 3100/4000
        (no_padding_xmin, no_padding_ymax - 3375) # 2170
    ]
    output_name = "sea"

    # Convert the list into array to do transformation
    boundary_coordinates = np.array(polygon_padding).astype('float64')

    # Convert into shapely geometry Polygon
    shapely_polygon = Polygon(boundary_coordinates)

    # Create geopandas dataframe
    polygon_dataframe = gpd.GeoDataFrame(geometry=[shapely_polygon], crs=2193)

    # Write into shape file
    # Notice that the folder <bathy_{number_simulation}> was created in module demRaster,
    # function <necessary_files>
    pathlib.Path(
        fr"{bathy_shapefile}\\{output_name}\\bathy_{number_simulation}"
    ).mkdir(parents=True, exist_ok=True)
    polygon_dataframe.to_file(
        fr"{bathy_shapefile}\\{output_name}\\bathy_{number_simulation}\\polygon_{number_simulation}.shp", crs=2193
    )


def nc_to_tiff(
    number_simulation
):
    """
    @Definition:
                A function to convert netCDF file to GeoTiff file
    @References:
                None.
    @Arguments:
                number_simulation (string):
                            A string to identify the order of simulation (should be angle, x, y)
    @Returns:
                None.
    """
    # DEM ----------------------------------------------
    # Convert raster from NetCDF file into GeoTiff file
    raster_nc = rxr.open_rasterio(
        fr"{dem_nc_path}\\generated_dem_bathy_{number_simulation}.nc"
    )

    raster_nc.z.rio.to_raster(
        fr"{dem_tiff_path}\\generated_dem_bathy_{number_simulation}.tif", cache=False
    )

    # MANNING'S N ----------------------------------------
    # Convert raster from NetCDF file into GeoTiff file
    raster_nc = rxr.open_rasterio(
        fr"{n_nc_path}\\generated_n_bathy_{number_simulation}.nc"
    )

    raster_nc.rio.to_raster(
        fr"{n_tiff_path}\\generated_n_bathy_{number_simulation}.tif", cache=False
    )


def value_sea_change(
        number_simulation,
        polygon=True
):
    """
    @Definition:
                A function to change the values of sea area
    @References:
                None.
    @Arguments:
                number_simulation (string):
                            A string to identify the order of simulation (should be angle, x, y)
                polygon (boolean):
                            True means using polygon to change the values of sea area
                            False means using tide level to change the values of sea area

    @Return:
                None.
    """

    if polygon:
        # Create polygon shapefile for sea
        polygon_generation(number_simulation)

        # Get directory of shapefile and file that has values need changing
        shapefile_sea_dir = fr"{bathy_shapefile}\\sea\\bathy_{number_simulation}\\polygon_{number_simulation}.shp"

        # DEM ---------------------------
        dem_need_changing_dir = fr"{dem_tiff_path}\\generated_dem_bathy_{number_simulation}.tif"
        value_change(shapefile_sea_dir, dem_need_changing_dir, -999, True)

        # MANNING'S N ---------------------
        n_need_changing_dir = fr"{n_tiff_path}\\generated_n_bathy_{number_simulation}.tif"
        value_change(shapefile_sea_dir, n_need_changing_dir, -999, True)

        # STARTDEPTH --------------------
        startdepth_need_changing_dir = fr"{startdepth_tiff_path}\\startdepth_{number_simulation}.tif"
        value_change(shapefile_sea_dir, startdepth_need_changing_dir, 0, True)

    else:
        # Get tide flow data
        tide_flow_df = pd.read_csv(fr"{other_data}\\tide_flow_data.csv", sep=',')
        tide_flow_df['DateTime'] = pd.to_datetime(tide_flow_df['DateTime'], format="%Y-%m-%d %H:%M:%S", utc=False)

        # Get dem file from GeoTiff file and make a copy of it
        dem = rxr.open_rasterio(fr"{dem_tiff_path}\\generated_dem_bathy_{number_simulation}.tif")
        dem_copy = dem.copy(deep=True)

        # Get Manning's n file from GeoTiff file and make a copy of it
        n = rxr.open_rasterio(fr"{n_tiff_path}\\generated_n_bathy_{number_simulation}.tif")
        n_copy = n.copy(deep=True)

        # Get startdepth from GeoTiff file and make a copy of it
        startdepth = rxr.open_rasterio(fr"{startdepth_tiff_path}\\startdepth_{number_simulation}.tif")
        startdepth_copy = startdepth.copy(deep=True)

        # Get sea area
        seaarea = dem_copy * np.nan
        seaarea = seaarea.where(dem_copy >= (np.min(tide_flow_df.Level)+0.2), other=1) # value 0.2 here can be
        # adjusted manually and accordingly

        # Remove sea area
        dem_xr_domain_copy = dem_copy.where(seaarea.values != 1, other=np.nan) # DEM
        startdepth_copy = startdepth_copy.where(seaarea.values != 1, other=np.nan) # STARTDEPTH
        n_copy = n_copy.where(seaarea.values != 1, other=np.nan) # MANNING'S N

        # Write out file
        dem_xr_domain_copy.rio.to_raster(fr"{dem_tiff_path}\\generated_dem_bathy_{number_simulation}.tif", cache=False)
        startdepth_copy.rio.to_raster(fr"{startdepth_tiff_path}\\startdepth_{number_simulation}.tif", cache=False)
        n_copy.rio.to_raster(fr"{n_tiff_path}\\generated_n_bathy_{number_simulation}.tif", cache=False)


def convert_to_asc(number_simulation):

    """
    @Definition:
                A function to change the values of padding
    @References:
                None.
    @Arguments:
                number_simulation (string):
                            A string to identify the order of simulation (should be angle, x, y)
    @Return:
                None.
    """

    # DEM ---------------------------------
    # Convert GeoTiff file into ASCII file
    tiff_dem_func = rxr.open_rasterio(
        fr"{dem_tiff_path}\\generated_dem_bathy_{number_simulation}.tif"
    )
    nodata_tiff_dem_func = tiff_dem_func.rio.write_nodata(-999, inplace=True)
    crs_tiff_dem_func = nodata_tiff_dem_func.rio.write_crs(2193)
    crs_tiff_dem_func.rio.to_raster(fr"{dem_asc_path}\\generated_dem_bathy_{number_simulation}.asc")

    # MANNING'S N ---------------------------
    # Convert GeoTiff file into ASCII file
    tiff_n_func = rxr.open_rasterio(
        fr"{n_tiff_path}\\generated_n_bathy_{number_simulation}.tif"
    )
    tiff_n_func_upperlimit = tiff_n_func.where(tiff_n_func.values < 1000, 0)
    tiff_n_func_lowerlimit = tiff_n_func_upperlimit.where(tiff_n_func.values > -1000, -999)
    nodata_tiff_n_func = tiff_n_func_lowerlimit.rio.write_nodata(-999)
    nodata_tiff_n_func.rio.to_raster(fr"{n_asc_path}\\generated_n_bathy_{number_simulation}.asc")

    # STARTDEPTH --------------------------
    # Convert GeoTiff file into ASCII file
    tiff_startdepth_func = rxr.open_rasterio(
        fr"{startdepth_tiff_path}\\startdepth_{number_simulation}.tif"
    )
    nodata_tiff_startdepth_func = tiff_startdepth_func.rio.write_nodata(-999, inplace=True)
    crs_tiff_startdepth_func = nodata_tiff_startdepth_func.rio.write_crs(2193)
    crs_tiff_startdepth_func.rio.to_raster(fr"{startdepth_asc_path}\\startdepth_{number_simulation}.asc")