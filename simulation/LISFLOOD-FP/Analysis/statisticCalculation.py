# Prepare packages -----------------------------------------------------------------------------------------------------
from folder import *                                # For paths of sub-folders
import geopandas as gpd                             # For manipulating geopandas dataframe
import pandas as pd                                 # For manipulating dataframe
from shapely import wkt                             # For reading geometry in csv files
from impactCalculation import raster_generation     # For raster generation
# ----------------------------------------------------------------------------------------------------------------------

def statistic_calculation(
    dataset_func,
    calculation_option,
    flood_rate
):
    """
    @Definition:
                A function to calculate mean, standard deviation, coefficient of variation, and proportion of
                simulations of each cell being inundated (cell)
    @References:
                https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.std.html
                https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.mean.html
    @Arguments:
                dataset_func (pandas dataframe):
                                         A full dataset including unnecessary information
                calculation_option (string):
                                        "mean" means to calculate mean
                                        "sd" means to calculate standard deviation
                                        "cv" means to calculate coefficient of variation
                                        "cell" means to calculate probability of each pixel being inundated
                flood_rate (float):
                                        A rate where the flood was defined
    @Returns:
                new_dataframe (pandas dataframe):
                                        A dataframe includes x, y coordinates and filtered water depth values
    """
    # Remove geometry data (x, y coordinates)
    nogeo_data_func = dataset_func.drop([
        'x_coord',
        'y_coord'
    ], axis=1)

    # Calculate mean data
    mean_data_func = nogeo_data_func.copy()
    mean_data_func['mean'] = mean_data_func.mean(axis=1)

    # Remove pixels having values lower than <flood_rate> values
    copy_mean_data_func = mean_data_func.copy()
    copy_mean_data_func.loc[
        copy_mean_data_func['mean'] < flood_rate, ['mean']
        ] = -999

    # Manipulate and calculate necessary information
    if calculation_option == 'mean':
        filter_data_func = copy_mean_data_func['mean']

    elif calculation_option == 'sd':
        # Calculate standard deviation
        sd_data_func = nogeo_data_func.copy()
        sd_data_func['sd'] = sd_data_func.std(axis=1)

        # Copy standard deviation and filter data
        copy_sd_data_func = sd_data_func.copy()
        copy_sd_data_func.loc[
            copy_mean_data_func['mean'] == -999, ['sd']
        ] = -999
        filter_data_func = copy_sd_data_func['sd']

    elif calculation_option == 'cv':
        # Calculate coefficient of variation
        cv_data_func = nogeo_data_func.copy()
        cv_data_func['mean'] = cv_data_func.mean(axis=1)
        cv_data_func['sd'] = cv_data_func.std(axis=1)
        cv_data_func['cv'] = cv_data_func['sd'] / cv_data_func['mean'] * 100

        # Copy coefficient of variation and filter data
        copy_cv_data_func = cv_data_func.copy()
        copy_cv_data_func.loc[copy_mean_data_func['mean'] == -999, ['cv']] = -999
        filter_data_func = copy_cv_data_func['cv']

    else:
        # Calculate probability of each location getting inundated
        cell_data_func = nogeo_data_func.copy()
        cell_data_func['cell'] = (
                nogeo_data_func.shape[1] - (nogeo_data_func <= flood_rate).sum(axis=1)
        ) / nogeo_data_func.shape[1] * 100

        # Filter data
        cell_data_func.loc[copy_mean_data_func['mean'] == -999, ['cell']] = -999
        filter_data_func = cell_data_func['cell']

    # Create new dataframe with calculated information
    new_database = {'x': dataset_func['x_coord'],
                    'y': dataset_func['y_coord'],
                    f"{calculation_option}": filter_data_func}
    new_dataframe = pd.DataFrame(new_database)

    # Return new dataframe
    return new_dataframe


def area_calculation(dataset_func, resolution_func, flood_rate):
    """
    @Definition:
                A function to calculate area being flooded
    @References:
                None.
    @Arguments:
                dataset_func (pandas dataframe):
                                        A full dataset including unnecessary information
                resolution_func (int):
                                        Resolution value in meter
                flood_rate (float):
                                        A rate where the flood was defined
    @Returns:
                area_dataframe (pandas dataframe):
                                        Dataframe of areas of simulations
    """
    # Remove geometry data (x, y coordinates)
    nogeo_data_func = dataset_func.drop([
        'x_coord', 'y_coord'
    ], axis=1)

    # Get area dictionary
    area_dict = {}

    # Get area of each simulation
    for each_col in range(len(nogeo_data_func.columns)):
        # Get column num of each simulation
        column = nogeo_data_func[nogeo_data_func.columns[each_col]]

        # Get quantity of flooded cells
        num_flooded_cell = len(column[column > flood_rate])

        # Get area of each cell
        area_each_cell = resolution_func ** 2

        # Calculate area of each simulation
        area_func = num_flooded_cell * area_each_cell

        # Add each area to area dictionary
        area_dict[f"{nogeo_data_func.columns[each_col]}"] = area_func

    # Write dictionary into dataframe
    area_dataframe = pd.DataFrame(data=area_dict, index=[0])

    return area_dataframe


def building_calculation(dataset_func, building_path,
                         onepolygon_bathy_path=None):
    """
    @Definition:
                A function to calculate number of buildings being inundated
    @References:
                https://shapely.readthedocs.io/en/stable/manual.html
                https://stackoverflow.com/questions/40385782/make-a-union-of-polygons-in-geopandas-or-shapely-into-a-single-geometry
                https://shapely.readthedocs.io/en/stable/manual.html#shapely.ops.cascaded_union
                https://gis.stackexchange.com/questions/224496/creating-spatial-join-between-points-and-polygons-in-geopandas
                https://stackoverflow.com/questions/15943769/how-do-i-get-the-row-count-of-a-pandas-dataframe
                https://gis.stackexchange.com/questions/52705/how-to-write-shapely-geometries-to-shapefiles
                https://data.linz.govt.nz/layer/53353-nz-street-address/
    @Arguments:
                dataset_func (pandas dataframe):
                                        A full dataset including unnecessary information
                building_path (string):
                                        Path of file containing building polygons
                onepolygon_bathy_path (string):
                                        If specify, a string of MAINDIR will be provided to change the version path.
                                        Otherwise, one version will be kept over iteration and causes error file not
                                        existed.
                                        If None, a path will be automatically chosen.
    @Returns:
                building_dataframe (pandas dataframe):
                                        Dataframe of simulations' flooded buildings
    """

    # Get dictionary of number of buildings being inundated
    building_dict = {}

    # Get building polygons dataframe
    building_gpf = gpd.GeoDataFrame.from_file(building_path)

    # Remove geometry data (x, y coordinates)
    nogeo_data_func = dataset_func.drop(['x_coord', 'y_coord'], axis=1)

    # Get area of each simulation
    for each_col in range(len(nogeo_data_func.columns)):
        # Get column num of each simulation
        column = nogeo_data_func.columns[each_col]

        # Get flowdepth map which were converted to one big polygon
        if onepolygon_bathy_path is None:
            flowdepthmap_gpf = pd.read_csv(
                fr"{onepolygon_bathy}\\flowdepth_polygon_{column}.csv", engine='pyarrow'
            )
        else:
            flowdepthmap_gpf = pd.read_csv(
                fr"{onepolygon_bathy_path}\\5_analysis\\untransformed_impact\\onepolygon_nobackground\\flowdepth_polygon_{column}.csv", engine='pyarrow'
            )
        flowdepthmap_gpf['geometry'] = flowdepthmap_gpf['geometry'].apply(wkt.loads)
        flowdepthmap_gpf = gpd.GeoDataFrame(flowdepthmap_gpf, crs='epsg:2193')

        # Get number of buildings being inundated by merging one big polygon with building polygons
        building_in_polygon = gpd.sjoin(building_gpf, flowdepthmap_gpf, predicate='within')

        # Add that number to dictionary of buildings being inundated
        building_dict[f"{column}"] = len(building_in_polygon.index)

    # Write dictionary into dataframe
    building_dataframe = pd.DataFrame(data=building_dict, index=[0])

    return building_dataframe

def calculation_dict(dataset_func, resolution,
                     building_path, flood_rate,
                     raster_generation_command=True,
                     onepolygon_bathy_path=None,
                     ):
    """
    @Definition:
                A function to calculate number of buildings being inundated
    @References:
                https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.std.html
                https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.mean.html
    @Arguments:
                dataset_func (pandas dataframe):
                                        A full dataset including unnecessary information
                resolution (int):
                                        Resolution value in meter
                building_path (string):
                                        Path of file containing building polygons
                flood_rate (float):
                                        A rate where the flood was defined
                raster_generation_command (boolean):
                                        If True, raster_generation function will be processed
                                        If False, raster_generation function will not be processed
                onepolygon_bathy_path (string):
                                        If specify, a string of MAINDIR will be provided to change the version path.
                                        Otherwise, one version will be kept over iteration and causes error file not
                                        existed.
                                        If None, a path will be automatically chosen.
    @Returns:
                calculation_dict_set (dictionary):
                                        A dictionary contains all statistical results
    """
    statistic_method = ["mean", "sd", "cv", "cell"]

    calculation_dict_set = {}

    for each_method in statistic_method:
        # Statistical calculation
        statistic_result = statistic_calculation(dataset_func, each_method, flood_rate)

        # Store in dictionary
        calculation_dict_set[each_method] = statistic_result

        # If raster_generation_command = True, we will generate rasters of mean, sd, ...
        # If False, we will not generate rasters of mean, sd, ...
        # In comparison modules, we might not need to generate the raster again
        if raster_generation_command:
            raster_generation(
                statistic_result['x'],
                statistic_result['y'],
                statistic_result[each_method],
                filename=each_method
            )
        else:
            pass

    # Add area result
    calculation_dict_set['area'] = area_calculation(dataset_func, resolution, flood_rate)

    # Add building result
    calculation_dict_set['building'] = building_calculation(dataset_func, building_path,
                                                            onepolygon_bathy_path)

    return calculation_dict_set