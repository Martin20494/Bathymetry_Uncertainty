# Prepare packages -----------------------------------------------------------------------------------------------------
from folder import *                                # For paths of sub-folders

# Packages for transformation
from numba import guvectorize, float64              # For speeding up the time of the code running

# Packages for general data manipulation
import numpy as np                                  # For all calculation and data array/matrices manipulation
import pandas as pd                                 # For dealing with csv file within pandas dataframe

# For handling raster
import rioxarray as rxr                             # For manipulating pixel values, spatial attributes, and raster files
import xarray as xr                                 # For reading raster file

# Packages for geometry manipulation
import geopandas as gpd                             # For creating polygon
from shapely.geometry import Point                  # For creating point from coordinates of pixel middle

# Packages for searching depth values of point sample's coordinates
from shapely.strtree import STRtree                 # For building up a tree for searching points within this tree
from shapely import wkt                             # For reading geometry in csv files

# Packages for multiprocessing
from functools import partial                       # For containing many variables
import multiprocessing                              # For parallelising

# Packages for directory manipulation
import glob                                         # For collect a list of all files stored in a folder
from pathlib import Path                            # For manipulating paths' names
# ----------------------------------------------------------------------------------------------------------------------


# Remove all warnings  -------------------------------------------------------------------------------------------------
# Reference: https://numba.pydata.org/numba-doc/dev/reference/deprecation.html#suppressing-deprecation-warnings
from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning, NumbaWarning
import warnings

warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaWarning)
# ----------------------------------------------------------------------------------------------------------------------


# CONVERT RECTANGLE ARRAY INTO 3D XYZ ARRAY ############################################################################
def array_creation(data_array, value, switch=False):
    """
    @Definition:
                A function to create an array of coordinate
    @References:
                None
    @Arguments:
                data_array (array):
                                        Original array with full data
                value (string):
                                        Name of value - 'x' or 'y'
                switch (boolean):
                                        Switch the shape of array (shape x,y or shape y,x)
    @Returns:
                new_array (array):
                                        An array of coordinate
    """
    # Read x, y values into arrays
    arr_x = data_array.x.values
    arr_y = data_array.y.values

    if switch:
        # Create zero array
        new_array = np.zeros((arr_y.shape[0], arr_x.shape[0]))
    else:
        # Create zero array
        new_array = np.zeros((arr_x.shape[0], arr_y.shape[0]))

    # Get full number of x or y values
    if value == "x":
        for i in range(arr_x.shape[0]):
            for j in range(arr_y.shape[0]):
                new_array[j, i] = arr_x[i]
        return new_array

    else:
        for i in range(arr_x.shape[0]):
            for j in range(arr_y.shape[0]):
                new_array[j, i] = arr_y[j]
        return new_array

def xyz_array(
    dataset_z_rxr,
    switch=False
):
    """
    @Definition:
                A function to create an array of coordinate
    @References:
                None
    @Arguments:
                dataset_z_rxr (array):
                                    Rioxarray array that contains elevation or flowdepth values including padding
    @Returns:
                new_array (array):
                                        An array of coordinate
    """

    # Create full number of values of x, y, z coordinates
    array_x = array_creation(dataset_z_rxr, 'x', switch)
    array_y = array_creation(dataset_z_rxr, 'y', switch)
    array_z = dataset_z_rxr.isel(band=0).values

    # Flatten x, y, z arrays
    flatten_x = array_x.flatten()
    flatten_y = array_y.flatten()
    flatten_z = array_z.flatten()

    # Put all x, y, z into one array
    full_dataset = np.vstack((flatten_x, flatten_y, flatten_z)).transpose()

    return full_dataset

# END CONVERT RECTANGLE ARRAY INTO 3D XYZ ARRAY ########################################################################

# POINT SAMPLE #########################################################################################################
def clip_padding_xyzdataset(
    dataset
):
    """
    @Definition:
                A function to clip 3D xyz array
    @References:
                None.
    @Arguments:
                dataset (array):
                                    An array with x, y, z values or xyz dataset
    @Returns:
                adjusted_dataset_clip2 (3D array):
                                    A dataset without padding.
                                    The coordinates are adjusted 0.000001 degree and meter for value extraction
    """
    # Get boundary coordinates
    raster_origin = rxr.open_rasterio(fr"{original_lidar_path}\\no_padding\\no_padding.nc")
    xmin, ymin, xmax, ymax = raster_origin.z.rio.bounds()

    # if clipping using rectangle shape
    dataset_clip = dataset[(xmin <= dataset[:, 0])
                           & (xmax >= dataset[:, 0])
                           & (ymin <= dataset[:, 1])
                           & (ymax >= dataset[:, 1])]

    return dataset_clip

def point_sample_generation(
    switch=False
):
    """
    @Definition:
                A function to create a point sample file. The coordinates will be extracted from this file
                to collect the elevation data with the same coordinates in all simulations
    @References:
                https://geopandas.org/docs/reference/api/geopandas.GeoDataFrame.to_crs.html
                https://shapely.readthedocs.io/en/stable/manual.html
                https://geopandas.org/docs/user_guide/io.html
                https://sgillies.net/2014/01/18/getting-shapes-of-raster-features-with-rasterio.html
    @Arguments:
                extract_name (string):
                                        Name of a specific output among all flood model outputs
                number_simulation (string):
                                        A string to identify the order of simulation (should be angle, x, y)
                switch (boolean):
                                        Switch the shape of array (shape x,y or shape y,x)
    @Returns:
                point_df (geopandas DataFrame):
                                        A geospatial dataframe containing coordinates (geometry) of point sample
    """
    # Get dataset x rioxarray
    dataset_full_rxr = rxr.open_rasterio(
        fr"{original_lidar_path}\\no_padding\\no_padding.nc"
    )
    dataset_z_rxr = dataset_full_rxr.z

    # Get xyz dataset and clip the padding
    raster_xyz = xyz_array(dataset_z_rxr, switch=switch)
    clipped_padding_xyz = clip_padding_xyzdataset(raster_xyz)

    # Convert x, y coordinates array into shapely geometry
    # General idea here is to convert pixel middle point (containing x, y coordinates and elevation data)
    # to shapely geometry Point (only contains x, y coordinates)
    point_geo_values = [
        Point(clipped_padding_xyz[i, 0], clipped_padding_xyz[i, 1]) for i in range(clipped_padding_xyz.shape[0])
    ]

    # Develop geopandas dataframe
    data_coords = {
        "x_coord": clipped_padding_xyz[:, 0],
        "y_coord": clipped_padding_xyz[:, 1]
    }
    point_df = gpd.GeoDataFrame(
        data=data_coords,
        geometry=point_geo_values,
        crs=2193
    )

    return point_df

# END POINT SAMPLE #####################################################################################################



# DEPTH VALUES #########################################################################################################
def get_flowdepth_values(point_sample_csv, flowdepth_file):
    """
    @Definition:
                A function to extract flowdepth at given coordinates from any geometry files
    @References:
                https://gis.stackexchange.com/questions/102933/more-efficient-spatial-join-in-python-without-qgis-arcgis-postgis-etc/103066#103066
                https://gis.stackexchange.com/questions/121469/get-shapefile-polygon-attribute-value-at-a-specific-point-using-python-e-g-via
                https://stackoverflow.com/questions/59030022/checking-whether-point-is-within-polygon-returns-wrong-results-in-shapely
                https://gis.stackexchange.com/questions/119919/maximizing-code-performance-for-shapely
                https://gis.stackexchange.com/questions/42931/rtree-python-polygon-index
                https://gis.stackexchange.com/questions/227474/rtree-spatial-index-does-not-result-in-faster-intersection-computation
                https://rtree.readthedocs.io/en/latest/tutorial.html
                https://sgillies.net/2014/01/18/getting-shapes-of-raster-features-with-rasterio.html
                https://stackoverflow.com/questions/20474549/extract-points-coordinates-from-a-polygon-in-shapely

                As rtree cannot be used with multiprocessing, STRtree should be used instead

                https://gis.stackexchange.com/questions/353619/shapely-with-rtree-versus-strtree
                https://shapely.readthedocs.io/en/stable/manual.html
                https://gis.stackexchange.com/questions/396615/shapely-geometry-how-to-preserve-refrence-to-original-feature
    @Arguments:
                point_sample_csv (geopandas DataFrame):
                                        A geospatial dataframe contains point sample's coordinates
                flowdepth_file (string):
                                        A directory of a simulation csv file of flowdepth
    @Returns:
                new_array (array):
                                        An array of coordinate
    """
    # Get the files
    infile_csv = pd.read_csv(flowdepth_file, engine='pyarrow')
    infile_csv['geometry'] = infile_csv['geometry'].apply(wkt.loads)
    infile_csv = gpd.GeoDataFrame(infile_csv, crs='epsg:2193')

    # Convert dataframe into lists
    list_file_polygon = infile_csv['geometry'].tolist()
    list_file_point = point_sample_csv['geometry'].tolist()

    # Build up the tree
    tree = STRtree(list_file_polygon)

    # Get index (id) of polygon list
    index_by_id = dict((id(pt), i) for i, pt in enumerate(list_file_polygon))

    # Empty list
    flowdepth_list = []

    # Get intersected tuple
    for num_point in range(len(list_file_point)):
        intersect_tuple = [(index_by_id[id(pts)], pts.wkt) for pts in tree.query(list_file_point[num_point]) if
                           list_file_point[num_point].within(pts)]
        flowdepth_list.append(infile_csv.iloc[intersect_tuple[0][0]]['depth'])

    return flowdepth_list


def get_flowdepth_parallelism(
    num_processes
):
    """
    @Definition:
                A function to extract flowdepth at given coordinates from any geometry files by multiprocessing
    @References:
                https://stackoverflow.com/questions/43175382/python-create-a-pandas-data-frame-from-a-list
    @Arguments:
                num_processes (int):
                                        A number of process for the parallelism
    @Returns:
                new_array (array):
                                        An array of coordinate
    """
    # Get necessary parameters -----------------------

    # Get point sample dataframe
    point_sample_df = point_sample_generation(True)

    # Get all csv folders
    all_csv_folders = glob.glob(fr"{bathy_flowdepth}\\bathy_*")

    # Get all column names
    column_names = [Path(all_csv_folders[i]).stem for i in range(len(all_csv_folders))]

    # Get flowdepth list -----------------------------

    # List all parameters
    point_sample_csv = point_sample_df

    # Design func parameters
    func = partial(
        get_flowdepth_values,
        point_sample_csv,
    )

    # Design the pool and execute the multiprocessing
    with multiprocessing.Pool(processes=num_processes) as pool:
        flowdepth_df = pd.DataFrame(
            pool.map(func, all_csv_folders),
            index=column_names
        ).T
    pool.close()
    pool.join()

    # Add x coordinate column
    flowdepth_df.insert(0, 'y_coord', value=point_sample_csv['y_coord'])
    flowdepth_df.insert(0, 'x_coord', value=point_sample_csv['x_coord'])

    # Write out file
    flowdepth_df.to_csv(fr"{csv_bathy}\\all_simulations.csv", index=False)

    return flowdepth_df


# END DEPTH VALUES #####################################################################################################