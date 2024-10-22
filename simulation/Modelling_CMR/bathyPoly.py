# Prepare packages -----------------------------------------------------------------------------------------------------
# Packages for paths manipulation
from folder import *                              # For paths of sub-folders
import os                                         # For manipulating the directory path and to execute commands
import pathlib                                    # For manipulating the directory path

# Packages for untransformation
import numpy as np                                # For all calculation and data array/matrices manipulation

# Packages for unrotating and untranslating
import xarray as xr
import rioxarray as rxr                           # For reading flowdepth files
import rasterio                                   # For reading and manipulating spatial data
import rasterio.features                          # For vectorising features in array
from shapely.geometry import shape                # For manipulating spatial information (geometry) under GeoJSON format
import geopandas as gpd                           # For manipulating shape files


# Packages for multiprocessing
from functools import partial                       # For containing many variables
import multiprocessing                              # For parallelising

from valueChange import value_change

# ----------------------------------------------------------------------------------------------------------------------


def water_extraction(
    number_simulation,
    extract_name
):
    """
    @Definition:
                A function to extract and add crs into a specific flowdepth file
    @References:
                None.
    @Arguments:
                number_simulation (string):
                                            A string to identify the order of simulation (should be angle, x, y)
                extract_name (string):
                                            Name of a specific output among all flood model outputs
    @Returns:
                None.
    """
    # Extract required water file
    water_asc = rxr.open_rasterio(
        fr"{FPoutput_path}\\bathy_{number_simulation}\\{extract_name}"
    )

    # Add crs
    new_water = water_asc.rio.write_crs(2193)

    # Write out new flowdepth file
    if extract_name == 'out.max':
        new_water.rio.to_raster(fr"{extracted_wd}\\bathy_{extract_name}_{number_simulation}.nc")
    else:
        new_water.rio.to_raster(fr"{extracted_wse}\\bathy_{extract_name}_{number_simulation}.nc")


def wse_elev_combination(
    extract_name,
    number_simulation
):
    """
    @Definition:
                A function to create new elevation by combining water depth and elevation
    @References:
                None
    @Arguments:
                number_simulation (string):
                                            A string to identify the order of simulation (should be angle, x, y)
                extract_name (string):
                                            Name of a specific output among all flood model outputs
    @Returns:
                None.
    """

    # Get elev
    raster_elev = rxr.open_rasterio(fr"{dem_nc_path}\\generated_dem_bathy_{number_simulation}.nc").z

    # Get water depth
    raster_wd = rxr.open_rasterio(fr"{extracted_wd}\\bathy_out.max_{number_simulation}.nc")

    # Get water surface elevation
    wse_values = raster_elev.values + raster_wd.values

    # Create raster for new water surface elevation
    new_wse = xr.DataArray(
        data=wse_values[0],
        dims=['y', 'x'],
        coords={
            'x': (['x'], raster_wd.x.values),
            'y': (['y'], raster_wd.y.values)
        },
        attrs=raster_wd.attrs
    )

    # Set up crs and nodata
    new_wse.rio.write_crs("epsg:2193", inplace=True)
    new_wse.rio.write_nodata(-9999, inplace=True)

    # Write out
    new_wse.rio.to_raster(fr"{extracted_new_wse}\\bathy_new_wse_{number_simulation}.nc")



def river_filter(
    number_simulation,
    extract_name
):
    """
    @Definition:
                A function to filter river by polygon in flood model output
    @References:
                None.
    @Arguments:
                number_simulation (string):
                                            A string to identify the order of simulation (should be angle, x, y)
                extract_name (string):
                                            Name of a specific output among all flood model outputs
    @Returns:
                None.
    """
    if extract_name == "out.max":
        water_origin = rxr.open_rasterio(fr"{extracted_wd}\\bathy_{extract_name}_{number_simulation}.nc")
        filter_tiff = filter_extracted_wd_tiff
        filter_shp = filter_extracted_wd_shp
        filter_nc = filter_extracted_wd_nc
        filter_value = 0
    elif extract_name == 'out.mxe':
        water_origin = rxr.open_rasterio(fr"{extracted_wse}\\bathy_{extract_name}_{number_simulation}.nc")
        filter_tiff = filter_extracted_wse_tiff
        filter_shp = filter_extracted_wse_shp
        filter_nc = filter_extracted_wse_nc
        filter_value = -9999
    elif extract_name == 'elev':
        dem_elev = rxr.open_rasterio(fr"{dem_nc_path}\\generated_dem_bathy_{number_simulation}.nc")
        water_origin = dem_elev.z
        filter_tiff = filter_extracted_elev_tiff
        filter_shp = filter_extracted_elev_shp
        filter_nc = filter_extracted_elev_nc
        filter_value = -9999
    elif extract_name == 'new_wse':
        # Create new wse by combining water depth and elevation
        wse_elev_combination(
            extract_name,
            number_simulation
        )
        water_origin = rxr.open_rasterio(fr"{extracted_new_wse}\\bathy_{extract_name}_{number_simulation}.nc")
        filter_tiff = filter_extracted_new_wse_tiff
        filter_shp = filter_extracted_new_wse_shp
        filter_nc = filter_extracted_new_wse_nc
        filter_value = -9999
    else:
        water_origin = rxr.open_rasterio(fr"{n_nc_path}\\generated_n_bathy_{number_simulation}.nc")
        filter_tiff = filter_extracted_n_tiff
        filter_shp = filter_extracted_n_shp
        filter_nc = filter_extracted_n_nc
        filter_value = -9999

    # Convert netCDF file to Tiff
    water_origin.rio.to_raster(fr"{filter_tiff}\\bathy_{extract_name}_{number_simulation}.tif")

    # Convert geojson into shapefile
    river_geo = gpd.read_file(fr"{bathy_path}\\bathy_{number_simulation}\\river_polygon.geojson")
    river_geo.to_file(fr"{filter_shp}\\river_polygon_{number_simulation}.shp", crs=2193)

    # Remove river
    watertiff_path = fr"{filter_tiff}\\bathy_{extract_name}_{number_simulation}.tif"
    rivershp_path = fr"{filter_shp}\\river_polygon_{number_simulation}.shp"
    value_change(rivershp_path, watertiff_path, filter_value, True)

    # Convert back to netcdf
    watertiff_filter = rxr.open_rasterio(fr"{filter_tiff}\\bathy_{extract_name}_{number_simulation}.tif")
    if extract_name == 'elev':
        watertiff_filter = watertiff_filter.rio.write_nodata(-9999)
        watertiff_filter.rio.to_raster(fr"{filter_nc}\\bathy_{extract_name}_{number_simulation}.nc")
    elif extract_name == 'n':
        watertiff_filter = watertiff_filter.rio.write_nodata(-9999)
        watertiff_filter.rio.to_raster(fr"{filter_nc}\\bathy_{extract_name}_{number_simulation}.nc")
    else:
        watertiff_filter.rio.to_raster(fr"{filter_nc}\\bathy_{extract_name}_{number_simulation}.nc")


def polygon_bathy(
    number_simulation,
    extract_name
):
    """
    @Definition:
                A function to convert raster into polygons and un-transform them
    @References:
                https://docs.dea.ga.gov.au/notebooks/Frequently_used_code/Polygonise_pixel_edges.html
                https://gis.stackexchange.com/questions/187877/how-to-polygonize-raster-to-shapely-polygons/

                https://github.com/sgillies/affine/blob/master/affine/__init__.py#L178
                https://gis.stackexchange.com/questions/408386/rotating-raster-using-python
                https://pythonhosted.org/PyAgg/affine.m.html

                https://gis.stackexchange.com/questions/408386/rotating-raster-using-python
                https://gis.stackexchange.com/questions/350526/gdal-setgeotransform-issue
                https://corteva.github.io/rioxarray/stable/examples/convert_to_raster.html
                https://gdal.org/tutorials/geotransforms_tut.html#geotransforms-tut
                https://gdal.org/tutorials/raster_api_tut.html
    @Arguments:
                number_simulation (string):
                                            A string to identify the order of simulation (should be angle, x, y)
                extract_name (string):
                                            Name of a specific output among all flood model outputs
    @Returns:
                None.
    """
    # Extract array from flowdepth raster
    if extract_name == 'elev':
        bathy_water = bathy_elev # Choose the output
        raster_poly = rasterio.open(fr"{filter_extracted_elev_nc}\\bathy_{extract_name}_{number_simulation}.nc")

    elif extract_name == 'n':
        bathy_water = bathy_n # Choose the output
        raster_poly = rasterio.open(fr"{filter_extracted_n_nc}\\bathy_{extract_name}_{number_simulation}.nc")

    elif extract_name == 'out.max':
        bathy_water = bathy_wd # Choose the output
        raster_poly = rasterio.open(fr"{filter_extracted_wd_nc}\\bathy_{extract_name}_{number_simulation}.nc")

    elif extract_name == 'new_wse':
        bathy_water = bathy_new_wse  # Choose the output
        raster_poly = rasterio.open(fr"{filter_extracted_new_wse_nc}\\bathy_{extract_name}_{number_simulation}.nc")

    else:
        bathy_water = bathy_wse # Choose the output
        raster_poly = rasterio.open(fr"{filter_extracted_wse_nc}\\bathy_{extract_name}_{number_simulation}.nc")

    raster_array = raster_poly.read(1)
    raster_transform = raster_poly.transform
    raster_crs = raster_poly.crs

    # Extract parameters: id and depth
    id_pixels = np.arange(raster_array.size).reshape(raster_array.shape)
    depth_list = list(raster_array.flatten())

    # Find out raster centers
    bathy_vectors = rasterio.features.shapes(source=id_pixels.astype(np.int16),
                                             transform=raster_transform)
    # Vectorise features
    bathy_vectors_list = list(bathy_vectors)

    # Get geometry
    polygons_geometry_values = [
        shape(polygon_geometry) for polygon_geometry, value_geometry in bathy_vectors_list
    ]

    # Get id
    id_bathy_pixels_values = [
        id_value for id_polygon, id_value in bathy_vectors_list
    ]

    # Create BATHY database under geopandas dataframe (gdf) format
    bathy_data = {
        "id": id_bathy_pixels_values,
        "depth": depth_list
    }
    bathy_raster_poly_gdf = gpd.GeoDataFrame(
        data=bathy_data,
        geometry=polygons_geometry_values,
        crs=raster_crs
    )

    # Write out csv files
    bathy_raster_poly_gdf.to_csv(
        fr"{bathy_water}\\bathy_{extract_name}_{number_simulation}.csv", index=False
    )


def bathy_simulation(
    extract_name,
    ran_trans_i
):
    """
    @Definition:
                A function to extract, convert raster into polygons, and un-transform them
    @References:
                None.
    @Arguments:
                extract_name (string):
                                        Name of a specific output among all flood model outputs
                ran_trans_i (array):
                                        A 3D array contains values of angle, x, and y coordinates of points in tiles
                                        (iterating variable generated from multiprocessing)
    @Returns:
                None.
    """

    if extract_name == 'elev':
        # Remove river
        river_filter(
            ran_trans_i,
            extract_name
        )

        # Convert to and untransform polygons
        polygon_bathy(
            ran_trans_i,
            extract_name
        )

    elif extract_name == 'n':
        # Convert to and untransform polygons
        # Remove river
        river_filter(
            ran_trans_i,
            extract_name
        )

        # Convert to and untransform polygons
        polygon_bathy(
            ran_trans_i,
            extract_name
        )

    elif extract_name == 'new_wse':
        # Remove river
        river_filter(
            ran_trans_i,
            extract_name
        )

        # Convert to and untransform polygons
        polygon_bathy(
            ran_trans_i,
            extract_name
        )

    else:
        # # Get specific flowdepth
        # water_extraction(
        #     ran_trans_i,
        #     extract_name
        # )

        # Remove river
        river_filter(
            ran_trans_i,
            extract_name
        )

        # Convert to and untransform polygons
        polygon_bathy(
            ran_trans_i,
            extract_name
        )

def bathy_parallelism(
    extract_name,
    ran_trans,
    num_processes
):
    """
    @Definition:
                A function to extract, convert raster into polygons, and un-transform them by applying nested multiprocessing
    @References:
                None.
    @Arguments:
                extract_name (string):
                                        Name of a specific output among all flood model outputs
                ran_trans (array):
                                        A big array of small 3D arrays of transformation values (angle, x, y)
                num_processes (int):
                                        A number of process for the parallelism
    @Returns:
                None.
    """
    # List parameters
    extract_name = extract_name

    # Design a func to be used in multiprocessing
    func = partial(
        bathy_simulation,
        extract_name
    )

    # Design the pool and execute the multiprocessing
    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.map(func, [ran for ran in ran_trans])
    pool.close()
    pool.join()


