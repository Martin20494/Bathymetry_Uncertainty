# Prepare packages -----------------------------------------------------------------------------------------------------
# For controlling the paths
from folder import *                                # For paths of sub-folders

# For data handling generally
import numpy as np                                  # For all calculation and data array/matrices manipulation
import pandas as pd                                 # For manipulating dataframe

# For raster and vector handling
import rioxarray as rxr                             # For manipulating pixel values, spatial attributes, and raster files
import xarray as xr                                 # For writing arrays into raster files
import rasterio                                     # For manipulating raster files (mainly crs and transformation)
import rasterio.features                            # For manipulating raster features


# For polygon handling
import geopandas as gpd                             # For manipulating shape files
from shapely.geometry import shape                  # For manipulating spatial information (geometry) under GeoJSON format
from shapely.ops import unary_union                 # For combining all polygons into one

# Packages for multiprocessing
from functools import partial                       # For containing many variables
import multiprocessing                              # For parallelising
# ----------------------------------------------------------------------------------------------------------------------


# RASTER WRITING #######################################################################################################

def raster_conversion(
    x_func, y_func, z_func
):
    """
    @Definition:
                A function to convert !D dimension of array into raster array.
                Particularly, convert 1D x, 1D y, and 1D z into multi-dimensional arrays of x, y, z
    @References:
                https://stackoverflow.com/questions/41897544/make-a-contour-plot-by-using-three-1d-arrays-in-python
                https://stackoverflow.com/questions/41815079/pandas-merge-join-two-data-frames-on-multiple-columns
    @Arguments:
                x_func, y_func, z_func (arrays):
                                            1D x, 1D y, 1D z datasets
    @Returns:
                x_func, y_func, z_func (arrays):
                                            Multi-dimensional arrays of x, y, z
    """
    # Gather x, y, z datasets into a pandas dataframe (DATAFRAME of ACTUAL DATA)
    pd_dataframe_actual = pd.DataFrame(dict(
        x=x_func,
        y=y_func,
        z=z_func
    ))

    # Assign dataframe column names into variables
    xcol, ycol, zcol = 'x', 'y', 'z'

    # Sort actual dataframe according to x then y values
    pd_dataframe_sorted = pd_dataframe_actual.sort_values(by=[xcol, ycol])

    # Getting values of x, y, z under raster array format
    # unique() function is used to remove duplicates
    x_values_func = pd_dataframe_sorted[xcol].unique()
    y_values_func = pd_dataframe_sorted[ycol].unique()
    z_values_func = pd_dataframe_sorted[zcol].values.reshape(len(x_values_func), len(y_values_func)).T

    return x_values_func, y_values_func, z_values_func


def raster_generation(
    x_func, y_func, z_func,
    filename,
    savepath=None
):
    """
    @Definition:
                A function to convert 1D x, 1D y, 1D z arrays into multi-dimensional arrays of x, y, z.
                Then write out netCDF files
    @References:
                https://stackoverflow.com/questions/41897544/make-a-contour-plot-by-using-three-1d-arrays-in-python
                https://stackoverflow.com/questions/41815079/pandas-merge-join-two-data-frames-on-multiple-columns
    @Arguments:
                x_func, y_func, z_func (arrays):
                                            1D x, 1D y, 1D z datasets
                filename (string):
                                            Name of raster file
                savepath (string):
                                            None when saving statistical files such as mean, sd, cv, etc.
                                            If not None, a specific path for saving
    @Returns:
                None.
    """
    # Read original DEM raster without padding
    raster_origin = rxr.open_rasterio(
        fr"{original_lidar_path}\\no_padding\\no_padding.nc"
    )

    # Round x and y to get duplicates
    x_round = np.round(x_func, 3)
    y_round = np.round(y_func, 3)

    # Get x, y, z values under raster array (gridded array)
    x_values_func, y_values_func, z_values_func = raster_conversion(x_round, y_round, z_func)

    # Write x, y, z values into raster array
    raster_array = xr.DataArray(
        data=z_values_func,
        dims=['y', 'x'],
        coords={
            'x': (['x'], x_values_func),
            'y': (['y'], y_values_func)
        },
        attrs=raster_origin.attrs
    )

    # Set up crs and nodata
    raster_array.rio.write_crs("epsg:2193", inplace=True)
    raster_array.rio.write_nodata(-999, inplace=True)

    # Write into raster file (nc)
    if savepath is None:
        # Saving statistical files such as mean, sd, cv, etc.
        raster_array.rio.to_raster(
            fr"{raster_bathy}\\{filename}_flowdepth.nc"
        )
    else:
        # Saving files (raster) from merging flowdepth polygons and building polygons
        raster_array.rio.to_raster(fr"{savepath}\\{filename}.nc")

# END RASTER WRITING ###################################################################################################


# MERGE FLOWDEPTH POLYGONS AND BUILDING POLYGONS #######################################################################
def flowdepth_raster_nobackground(
    dataset_origin,
    flood_rate
):
    """
    @Definition:
                A function to generate raster of flowdepth map (excluding background). This flowdepth map will be
                used to convert into polygon (to get only one big polygon to merge with buildings rather than many
                pixel polygons)
    @References:
                https://stackoverflow.com/questions/41897544/make-a-contour-plot-by-using-three-1d-arrays-in-python
                https://stackoverflow.com/questions/41815079/pandas-merge-join-two-data-frames-on-multiple-columns
    @Arguments:
                dataset_origin (pandas DataFrame):
                                    A pandas dataframe contains x, y coordinates and depth values of all simulation
                                    after un-transformation
                flood_rate (float):
                                    A rate where the flood was defined
    @Returns:
                None.
    """
    # Copy dataset to avoid the original dataset being changed
    dataset = dataset_origin.copy(deep=True)

    for num_file in range(len(dataset.columns) - 2): # Minus 2 due to two columns of coordinates x and y
        # Convert background to -999
        col = dataset.columns[num_file + 2] # Plus 2 due to two columns of coordinates x and y
        dataset.loc[dataset[col] < flood_rate, [col]] = -999

        # Generate flowdepth raster excluding background
        raster_generation(
            dataset['x_coord'],
            dataset['y_coord'],
            dataset[f'{col}'],
            f"nobackground_{col}",
            oneraster_bathy
        )

def flowdepthmap_to_onepolygon(
    dataset_origin,
    flood_rate,
    column_num
):
    """
    @Definition:
                A function to convert flowdepth map (after excluding background) into one big polygon rather than
                many pixel polygons
    @References:
                https://numpy.org/doc/stable/reference/generated/numpy.copy.html
                https://stackoverflow.com/questions/19666626/replace-all-elements-of-python-numpy-array-that-are-greater-than-some-value
    @Arguments:
                dataset_origin (pandas DataFrame):
                                    A pandas dataframe contains x, y coordinates and depth values of all simulation
                                    after un-transformation
                flood_rate (float):
                                    A rate where the flood was defined
                column_num (int):
                                    Ordinal number of column of simulation defined in the dataset (pandas dataframe)
    @Returns:
                None.
    """
    # Copy dataset to avoid the original dataset being changed
    dataset = dataset_origin.copy(deep=True)

    # Get column value
    column = dataset.columns[column_num + 2]

    # Read flowdepth map (raster) after excluding background already
    raster_poly = rasterio.open(
        fr"{oneraster_bathy}\\nobackground_{column}.nc"
    )

    # Get information from original DEM raster
    raster_array = raster_poly.read(1)
    raster_transform = raster_poly.transform
    raster_crs = raster_poly.crs

    # Extract parameters: id
    id_pixels = np.arange(raster_array.size).reshape(raster_array.shape)

    # Extract parameters: depth
    value_list = dataset[f"{column}"].tolist()

    # Vectorise features
    vectors = rasterio.features.shapes(
        source=id_pixels.astype(np.int16),
        transform=raster_transform
    )
    vectors_list = list(vectors)

    # Get geometry
    poly_geometry = [shape(polygon) for polygon, value in vectors_list]

    # Get id
    id_poly = [id_val for id_poly, id_val in vectors_list]

    # Create database
    poly_database = {
        'id': id_poly,
        "depth": value_list
    }
    poly_dataframe = gpd.GeoDataFrame(
        data=poly_database,
        geometry=poly_geometry,
        crs=raster_crs
    )

    # Select only FLOODED area by using the filter rate
    flood_polygons_list = poly_dataframe[poly_dataframe['depth'] >= flood_rate]['geometry'].tolist()

    # Combine all polygons into one using function unary_union of shapely package
    one_flood_polygon = gpd.GeoSeries(unary_union(flood_polygons_list))

    # Write that one polygon again into geopandas dataframe
    one_poly_gdf = gpd.GeoDataFrame(
        data={"id": [1]},
        geometry=one_flood_polygon.tolist(),
        crs=2193
    )

    # Write out polygon
    one_poly_gdf.to_csv(
        fr"{onepolygon_bathy}\\flowdepth_polygon_{column}.csv", index=False
    )

def flowdepthmap_onepolygon_parallelism(
    dataset,
    flood_rate,
    num_processes
):
    """
    @Definition:
                A function to convert all flowdepth maps (after excluding background) into one-big polygons rather than
                many pixel polygons by multiprocessing
    @References:
                https://numpy.org/doc/stable/reference/generated/numpy.copy.html
                https://stackoverflow.com/questions/19666626/replace-all-elements-of-python-numpy-array-that-are-greater-than-some-value
    @Arguments:
                dataset (pandas DataFrame):
                                    A pandas dataframe contains x, y coordinates and depth values of all simulation
                                    after un-transformation
                flood_rate (float):
                                    A rate where the flood was defined
                num_processes (int):
                                    A number of process for the parallelism
    @Returns:
                None.
    """
    # Get list of all flowdepth files
    all_flowdepth_files = list(range(len(dataset.columns) - 2)) # Minus 2 for two columns of x and y coordinates

    # Design func parameters
    func = partial(
        flowdepthmap_to_onepolygon,
        dataset,
        flood_rate
    )

    # Design the pool and execute the multiprocessing
    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.map(func, all_flowdepth_files)
    pool.close()
    pool.join()

# END MERGE FLOWDEPTH POLYGONS AND BUILDING POLYGONS ###################################################################


