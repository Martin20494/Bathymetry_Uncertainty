# Prepare packages -----------------------------------------------------------------------------------------------------
# Packages for paths manipulation
from folder import *                              # For paths of sub-folders
import os                                         # For manipulating the directory path and to execute commands
import pathlib                                    # For manipulating the directory path

# Packages for untransformation
import numpy as np                                # For all calculation and data array/matrices manipulation

# Packages for unrotating and untranslating
import rioxarray as rxr                           # For reading flowdepth files
import rasterio                                   # For reading and manipulating spatial data
import rasterio.features                          # For vectorising features in array
from shapely.geometry import shape                # For manipulating spatial information (geometry) under GeoJSON format
from shapely.ops import unary_union               # For combining all polygons into one
import geopandas as gpd                           # For manipulating shape files


# Packages for multiprocessing
from functools import partial                       # For containing many variables
import multiprocessing                              # For parallelising

# ----------------------------------------------------------------------------------------------------------------------




def flowdepth_extraction(
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
    # Extract required flowdepth file
    flowdepth_asc = rxr.open_rasterio(
        fr"{FPoutput_path}\\bathy_{number_simulation}\\{extract_name}"
    )

    # Add crs
    new_flowdepth = flowdepth_asc.rio.write_crs(2193)

    # Write out new flowdepth file
    new_flowdepth.rio.to_raster(fr"{extracted_flowdepth}\\bathy_{extract_name}_{number_simulation}.nc")


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
                angle_func, x_translation_func, y_translation_func (float):
                                            Values to rotate and translate
                number_simulation (string):
                                            A string to identify the order of simulation (should be angle, x, y)
                extract_name (string):
                                            Name of a specific output among all flood model outputs
    @Returns:
                None.
    """
    # Extract array from flowdepth raster
    raster_poly = rasterio.open(fr"{extracted_flowdepth}\\bathy_{extract_name}_{number_simulation}.nc")
    raster_array = raster_poly.read(1)
    raster_transform = raster_poly.transform
    raster_crs = raster_poly.crs

    # Extract parameters: id and depth
    id_pixels = np.arange(raster_array.size).reshape(raster_array.shape)
    depth_list = list(raster_array.flatten())

    # Vectorise features
    bathy_vectors = rasterio.features.shapes(source=id_pixels.astype(np.int16),
                                             transform=raster_transform)

    # List the BATHY polygons to extract necessary parameters
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
        fr"{bathy_flowdepth}\\bathy_{extract_name}_{number_simulation}.csv", index=False
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
    # Get specific flowdepth
    flowdepth_extraction(
        ran_trans_i,
        extract_name
    )

    # Convert to polygons
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


