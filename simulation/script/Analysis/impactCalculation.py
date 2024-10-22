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
    extract_name,
    savepath=None,
    nodata=-9999
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
                extract_name (string):
                                            Name of a specific output among all flood model outputs
                savepath (string):
                                            None when saving statistical files such as mean, sd, cv, etc.
                                            If not None, a specific path for saving
                nodata (float):
                                            No data value
    @Returns:
                None.
    """
    # Choose outputs
    if extract_name == 'out.max':
        raster_bathy = wd_raster_bathy
    elif extract_name == 'out.mxe':
        raster_bathy = wse_raster_bathy
    elif extract_name == 'elev':
        raster_bathy = elev_raster_bathy
    else:
        raster_bathy = n_raster_bathy

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
    raster_array.rio.write_nodata(nodata, inplace=True)

    # Write into raster file (nc)
    if savepath is None:
        # Saving statistical files such as mean, sd, cv, etc.
        raster_array.rio.to_raster(
            fr"{raster_bathy}\\{filename}_water.nc"
        )
    else:
        # Saving files (raster) from merging water polygons and building polygons
        raster_array.rio.to_raster(fr"{savepath}\\{filename}.nc")

# END RASTER WRITING ###################################################################################################


# MERGE WATER POLYGONS AND BUILDING POLYGONS #######################################################################
def apply_zero_values(dataset_func, extract_name):
    """
    @Definition:
                A function to get zero values
    @References:

    @Arguments:
                dataset_func (pandas dataframe):
                                         A full dataset of water depth, water surface
                                         elevation, elevation including x and y coordinates
                extract_name (string):
                                        Name of a specific output among all flood model outputs

    @Returns:


    """
    # Copy dataset
    copy_dataset = dataset_func.copy(deep=True)

    # Replace -9999 with 0
    if extract_name == 'out.mxe':
        replaced_dataset = copy_dataset.replace([-9999], 0)
    else:
        replaced_dataset = copy_dataset

    # Remove geometry data (x, y coordinates)
    nogeo_dataset = replaced_dataset.drop([
        'x_coord',
        'y_coord'
    ], axis=1)

    # Calculate mean data to filter rows which are not flooded at all
    # and rows which have some simulations have flooding
    # and some do not have

    # Calculate mean
    mean_dataset = nogeo_dataset.copy(deep=True)
    mean_dataset['mean'] = nogeo_dataset.mean(axis=1)

    # Replace rows which are full of 0 values or mean = 0 with -9999
    mean_dataset.loc[mean_dataset['mean']==0, :] = -9999

    # Drop mean column after using it to filter
    new_dataset = mean_dataset.drop(mean_dataset.columns[-1], axis=1)

    # Return
    return new_dataset


def get_floodrate_index(flood_rate, csv_path=None):
    """
    @Definition:
                A function to collect floodrate index
    @References:

    @Arguments:
                flood_rate (float):
                                        A rate where the flood was defined
    @Returns:
                new_dataframe (pandas dataframe):
                                        A dataframe includes x, y coordinates and filtered water depth values
    """
    # Get water depth data
    if csv_path is None:
        df_wd = pd.read_csv(fr"{wd_csv_bathy}\\all_simulations.csv")
    else:
        df_wd = pd.read_csv(fr"{csv_path[1]}\\all_simulations.csv")

    # Get zero values. These dataframes have no geometry/coordinates
    df_wd_zero = apply_zero_values(df_wd, 'out.max')

    # Calculate the mean of water depth
    # to filter the flood rate for both
    # water surface elevation and water depth

    # Mean of water depth for filtering by water depth flood rate
    floodrate_df = df_wd_zero.copy(deep=True)
    floodrate_df['mean'] = df_wd_zero.mean(axis=1)

    # Filter flood rate and get its index (list format).
    # Indices in this list will be used to filter out
    # and replace values having the same indices in other dataframes
    # with -9999
    floodrate_index = floodrate_df['mean'].index[
        (floodrate_df['mean'] > 0)  # Larger than 0 to filter out -9999
        & (floodrate_df['mean'] < flood_rate)].tolist()  # Smaller than flood rate

    return floodrate_index


def water_raster_nobackground(
    dataset_origin,
    flood_rate,
    extract_name
):
    """
    @Definition:
                A function to generate raster of water map (excluding background). This water map will be
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
                extract_name (string):
                                    Name of a specific output among all flood model outputs
    @Returns:
                None.
    """
    # Choose flood outputs
    if extract_name == 'out.max':
        oneraster_bathy = wd_oneraster_bathy
    else:
        oneraster_bathy = wse_oneraster_bathy

    # Copy dataset to avoid the original dataset being changed and replace 0 with -9999
    # to create raster with nodata as -9999
    dataset = dataset_origin.replace([0], -9999)

    # Get flood rate index
    floodrate_index = get_floodrate_index(flood_rate)

    for num_file in range(len(dataset.columns) - 2): # Minus 2 due to two columns of coordinates x and y

        # Convert background to -9999
        col = dataset.columns[num_file + 2] # Plus 2 due to two columns of coordinates x and y
        dataset.loc[floodrate_index, [col]] = -9999

        # Generate water raster excluding background
        # Water surface elevation
        raster_generation(
            dataset['x_coord'],
            dataset['y_coord'],
            dataset[f'{col}'],
            f"nobackground_{col}",
            extract_name,
            oneraster_bathy
        )

def watermap_to_onepolygon(
    dataset_origin,
    flood_rate,
    extract_name,
    column_num
):
    """
    @Definition:
                A function to convert water map (after excluding background) into one big polygon rather than
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
                extract_name (string):
                                    Name of a specific output among all flood model outputs
                column_num (int):
                                    Ordinal number of column of simulation defined in the dataset (pandas dataframe)
    @Returns:
                None.
    """
    # Choose flood outputs
    if extract_name == 'out.max':
        onepolygon_bathy = wd_onepolygon_bathy
        oneraster_bathy = wd_oneraster_bathy
    else:
        onepolygon_bathy = wse_onepolygon_bathy
        oneraster_bathy = wse_oneraster_bathy


    # Copy dataset to avoid the original dataset being changed and replace 0 with -9999
    # to create raster with nodata as -9999
    dataset = dataset_origin.replace([0], -9999)

    # Get column value
    column = dataset.columns[column_num + 2]

    # Read water map (raster) after excluding background already
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
    floodrate_index = get_floodrate_index(flood_rate)
    poly_dataframe.loc[floodrate_index, 'depth'] = -9999
    flood_polygons_list = poly_dataframe.loc[poly_dataframe['depth']!=-9999, :]['geometry'].tolist()

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
        fr"{onepolygon_bathy}\\water_polygon_{column}.csv", index=False
    )

def watermap_onepolygon_parallelism(
    dataset,
    flood_rate,
    num_processes,
    extract_name
):
    """
    @Definition:
                A function to convert all water maps (after excluding background) into one-big polygons rather than
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
                extract_name (string):
                                    Name of a specific output among all flood model outputs
    @Returns:
                None.
    """
    # Get list of all water files
    all_water_files = list(range(len(dataset.columns) - 2)) # Minus 2 for two columns of x and y coordinates

    # Design func parameters
    func = partial(
        watermap_to_onepolygon,
        dataset,
        flood_rate,
        extract_name
    )

    # Design the pool and execute the multiprocessing
    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.map(func, all_water_files)
    pool.close()
    pool.join()

# END MERGE WATER POLYGONS AND BUILDING POLYGONS ###################################################################


