# Prepare packages -----------------------------------------------------------------------------------------------------
from folder import *                                # For paths of sub-folders
import geopandas as gpd                             # For manipulating geopandas dataframe
import pandas as pd                                 # For manipulating dataframe
import numpy as np                                  # For manipulating data
from sklearn.metrics import mean_squared_error      # For calculating RMSE
from shapely import wkt                             # For reading geometry in csv files
from impactCalculation import raster_generation, \
                              apply_zero_values, \
                              get_floodrate_index   # For raster generation and statistical measurements
# ----------------------------------------------------------------------------------------------------------------------

def statistic_calculation(
    extract_name,
    calculation_option,
    flood_rate,
    replace_value,
    csv_path=None
):
    """
    @Definition:
                A function to calculate mean, standard deviation, coefficient of variation, and proportion of
                simulations of each cell being inundated (cell)
    @References:
                https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.std.html
                https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.mean.html
    @Arguments:
                extract_name (string):
                                        Name of a specific output among all flood model outputs
                calculation_option (string):
                                        "mean" means to calculate mean
                                        "sd" means to calculate standard deviation
                                        "cv" means to calculate coefficient of variation
                                        "cell" means to calculate probability of each pixel being inundated
                flood_rate (float):
                                        A rate where the flood was defined
                replace_value (float):
                                        A value used to replace values lower than flood_rate
                 csv_path (string):
                                        Path link to csv files
    @Returns:
                new_dataframe (pandas dataframe):
                                        A dataframe includes x, y coordinates and filtered water depth values
    """
    # PREPARATION --------------------------------------------------------------

    # Get wse, wd, elev dataframes
    if csv_path is None:
        df_wse = pd.read_csv(fr"{wse_csv_bathy}\\all_simulations.csv")
        df_wd = pd.read_csv(fr"{wd_csv_bathy}\\all_simulations.csv")
        df_elev = pd.read_csv(fr"{elev_csv_bathy}\\all_simulations.csv")
        # df_n = pd.read_csv(fr"{n_csv_bathy}\\all_simulations.csv")
    else:
        df_wse = pd.read_csv(fr"{csv_path[0]}\\all_simulations.csv")
        df_wd = pd.read_csv(fr"{csv_path[1]}\\all_simulations.csv")
        df_elev = pd.read_csv(fr"{csv_path[2]}\\all_simulations.csv")

    # Get zero values. These dataframes have no geometry/coordinates
    df_wse_zero = apply_zero_values(df_wse, 'out.mxe')
    df_wd_zero = apply_zero_values(df_wd, 'out.max')

    # Get flood rate index (list format).
    # Indices in this list will be used to filter out
    # and replace values having the same indices in other dataframes
    # with -9999
    floodrate_index = get_floodrate_index(flood_rate, csv_path)

    # CALCULATION ----------------------------------------------------------------
    # Calculate mean, sd, cv, and proportion for water surface elevation
    # ------------------------------------------------------------------
    if extract_name == 'out.mxe':
        # Fill 0 values in water surface elevation
        # with elevation ---------------------------------

        # Make a copy of filled-zero dataframe
        df_wse_zero_copy = df_wse_zero.copy(deep=True)

        # Drop coordinate columns in elevation dataframe
        df_elev_nogeo = df_elev.drop([
            'x_coord',
            'y_coord'
        ], axis=1)

        # Fill 0 values
        for i in range(df_wse_zero_copy.shape[1]):
            # Get index of 0 values
            index_zero_elev = df_wse_zero_copy[df_wse_zero_copy.columns[i]].index[df_wse_zero_copy[df_wse_zero_copy.columns[i]] == 0].tolist()
            df_wse_zero_copy.loc[index_zero_elev, df_wse_zero_copy.columns[i]] = df_elev_nogeo.loc[index_zero_elev,
                                                                                         df_elev_nogeo.columns[i]]

        # Convert to new dataframe
        df_wse_elev = df_wse_zero_copy.copy(deep=True)

        # Write out csv file
        # Copy df_wse_elev
        df_wse_elev_copy = df_wse_elev.copy(deep=True)
        # Insert y
        df_wse_elev_copy.insert(loc=0, column='y', value=df_elev.y_coord.to_numpy())
        # Insert x
        df_wse_elev_copy.insert(loc=0, column='x', value=df_elev.x_coord.to_numpy())
        # Write out
        df_wse_elev_copy.to_csv(fr"{wse_csv_bathy}\\all_simulations_including_elev.csv", index=False)

        # Caclulate -------------------------------------

        # Calculate mean
        if calculation_option == 'mean':
            mean_wse_df = df_wse_elev.copy(deep=True)
            mean_wse_df['mean'] = df_wse_elev.mean(axis=1)

            # Filter flood rate
            mean_wse_df.loc[floodrate_index, 'mean'] = replace_value

            # Get values
            calculated_dataset = mean_wse_df['mean']

        # Calculate standard deviation
        elif calculation_option == 'sd':
            sd_wse_df = df_wse_elev.copy(deep=True)
            sd_wse_df['mean'] = df_wse_elev.mean(axis=1)
            sd_wse_df['sd'] = df_wse_elev.std(axis=1)

            # Filter unchanged values
            sd_wse_df.loc[sd_wse_df['mean'] == replace_value, 'sd'] = replace_value

            # Filter flood rate
            sd_wse_df.loc[floodrate_index, 'sd'] = replace_value

            # Get values
            calculated_dataset = sd_wse_df['sd']

        # Calculate coefficient of variation
        elif calculation_option == 'cv':
            cv_wse_df = df_wse_elev.copy(deep=True)
            cv_wse_df['mean'] = df_wse_elev.mean(axis=1)
            cv_wse_df['sd'] = df_wse_elev.std(axis=1)
            cv_wse_df['cv'] = cv_wse_df['sd'] / cv_wse_df['mean'] * 100

            # Filter unchanged values
            cv_wse_df.loc[cv_wse_df['mean'] == replace_value, 'cv'] = replace_value

            # Filter flood rate
            cv_wse_df.loc[floodrate_index, 'cv'] = replace_value

            # Get values
            calculated_dataset = cv_wse_df['cv']

        # Calculate proportion of simulations of each cell being inundated
        else:
            cell_wse_df = df_wse_elev.copy(deep=True)
            cell_wse_df['mean'] = df_wse_elev.mean(axis=1)
            cell_wse_df['cell'] = (df_wse_zero.shape[1] - (df_wse_zero.eq(0)).sum(axis=1)) / df_wse_zero.shape[1] * 100

            # Filter unchanged values
            cell_wse_df.loc[cell_wse_df['mean'] == replace_value, 'cell'] = replace_value

            # Filter flood rate
            cell_wse_df.loc[floodrate_index, 'cell'] = replace_value

            # Get values
            calculated_dataset = cell_wse_df['cell']

    # Calculate mean, sd, cv, and proportion for water depth
    elif extract_name == 'out.max':
        # Calculate mean
        if calculation_option == 'mean':
            mean_wd_df = df_wd_zero.copy(deep=True)
            mean_wd_df['mean'] = df_wd_zero.mean(axis=1)

            # Filter flood rate
            mean_wd_df.loc[floodrate_index, 'mean'] = replace_value

            # Get values
            calculated_dataset = mean_wd_df['mean']

        # Calculate standard deviation
        elif calculation_option == 'sd':
            sd_wd_df = df_wd_zero.copy(deep=True)
            sd_wd_df['mean'] = df_wd_zero.mean(axis=1)
            sd_wd_df['sd'] = df_wd_zero.std(axis=1)

            # Filter unchanged values
            sd_wd_df.loc[sd_wd_df['mean'] == replace_value, 'sd'] = replace_value

            # Filter flood rate
            sd_wd_df.loc[floodrate_index, 'sd'] = replace_value

            # Get values
            calculated_dataset = sd_wd_df['sd']

        # Calculate coefficient of variation
        elif calculation_option == 'cv':
            cv_wd_df = df_wd_zero.copy(deep=True)
            cv_wd_df['mean'] = df_wd_zero.mean(axis=1)
            cv_wd_df['sd'] = df_wd_zero.std(axis=1)
            cv_wd_df['cv'] = cv_wd_df['sd'] / cv_wd_df['mean'] * 100

            # Filter unchanged values
            cv_wd_df.loc[cv_wd_df['mean'] == replace_value, 'cv'] = replace_value

            # Filter flood rate
            cv_wd_df.loc[floodrate_index, 'cv'] = replace_value

            # Get values
            calculated_dataset = cv_wd_df['cv']

        # Calculate proportion of simulations of each cell being inundated
        else:
            cell_wd_df = df_wd_zero.copy(deep=True)
            cell_wd_df['mean'] = df_wd_zero.mean(axis=1)
            cell_wd_df['cell'] = (df_wd_zero.shape[1] - (df_wd_zero.eq(0)).sum(axis=1)) / df_wd_zero.shape[1] * 100

            # Filter unchanged values
            cell_wd_df.loc[cell_wd_df['mean'] == replace_value, 'cell'] = replace_value

            # Filter flood rate
            cell_wd_df.loc[floodrate_index, 'cell'] = replace_value

            # Get values
            calculated_dataset = cell_wd_df['cell']

    elif extract_name == 'elev':
        # Calculate mean of water depth to filter
        mean_wd_df = df_wd_zero.copy(deep=True)
        mean_wd_df['mean'] = df_wd_zero.mean(axis=1)

        # Filter flood rate
        mean_wd_df.loc[floodrate_index, 'mean'] = replace_value

        # Drop x, y coordinates
        df_elev_nocoords = df_elev.drop(columns=['x_coord', 'y_coord'])

        if calculation_option == 'mean':
            mean_elev_df = df_elev_nocoords.copy(deep=True)
            mean_elev_df['mean'] = df_elev_nocoords.mean(axis=1)

            # Filter by mean of water depth
            mean_elev_df.loc[mean_wd_df['mean'] == replace_value, 'mean'] = replace_value

            # Get values
            calculated_dataset = mean_elev_df['mean']

        elif calculation_option == 'sd':
            sd_elev_df = df_elev_nocoords.copy(deep=True)
            sd_elev_df['sd'] = df_elev_nocoords.std(axis=1)

            # Filter by mean of water depth
            sd_elev_df.loc[mean_wd_df['mean'] == replace_value, 'sd'] = replace_value

            # Get values
            calculated_dataset = sd_elev_df['sd']

        else:
            cv_elev_df = df_elev_nocoords.copy(deep=True)
            cv_elev_df['mean'] = df_elev_nocoords.mean(axis=1)
            cv_elev_df['sd'] = df_elev_nocoords.std(axis=1)
            cv_elev_df['cv'] = cv_elev_df['sd'] / cv_elev_df['mean'] * 100

            # Filter unchanged values (background, otherwise the background will be colored)
            cv_elev_df.loc[mean_wd_df['mean'] == replace_value, 'cv'] = replace_value

            # Filter flood rate
            cv_elev_df.loc[floodrate_index, 'cv'] = replace_value

            # Get values
            calculated_dataset = cv_elev_df['cv']

    else:
        # Calculate mean of water depth to filter
        mean_wd_df = df_wd_zero.copy(deep=True)
        mean_wd_df['mean'] = df_wd_zero.mean(axis=1)

        # Filter flood rate
        mean_wd_df.loc[floodrate_index, 'mean'] = replace_value

        # Drop x, y coordinates
        df_n_nocoords = df_n.drop(columns=['x_coord', 'y_coord'])

        # Calculate mean to filter
        mean_n_df = df_n_nocoords.copy(deep=True)
        mean_n_df['mean'] = df_n_nocoords.mean(axis=1)
        # Filter by mean of water depth
        mean_n_df.loc[mean_wd_df['mean'] == replace_value, 'mean'] = replace_value

        if calculation_option == 'mean':
            # # Filter one more time for any values larger than 1
            # mean_n_df['mean'][mean_n_df['mean'] >= 1] = replace_value

            # Get values
            calculated_dataset = mean_n_df['mean']

        elif calculation_option == 'sd':
            sd_n_df = df_n_nocoords.copy(deep=True)
            sd_n_df['sd'] = df_n_nocoords.std(axis=1)

            # Filter by mean of water depth
            sd_n_df.loc[mean_wd_df['mean'] == replace_value, 'sd'] = replace_value

            # # Filter one more time for any values larger than 1
            # sd_n_df['sd'][mean_n_df['mean'] >= 1] = replace_value

            # Get values
            calculated_dataset = sd_n_df['sd']

        else:
            cv_n_df = df_n_nocoords.copy(deep=True)
            cv_n_df['mean'] = df_n_nocoords.mean(axis=1)
            cv_n_df['sd'] = df_n_nocoords.std(axis=1)
            cv_n_df['cv'] = cv_n_df['sd'] / cv_n_df['mean'] * 100

            # Filter unchanged values (background, otherwise the background will be colored)
            cv_n_df.loc[mean_wd_df['mean'] == replace_value, 'cv'] = replace_value

            # Filter flood rate
            cv_n_df.loc[floodrate_index, 'cv'] = replace_value

            # Filter one more time for any values larger than 1
            cv_n_df['cv'][mean_n_df['mean'] >= 1] = replace_value

            # Get values
            calculated_dataset = cv_n_df['cv']

    # Create new dataset with calculated data and coordinates
    new_database = {'x': df_wse['x_coord'],
                    'y': df_wse['y_coord'],
                    f"{calculation_option}": calculated_dataset}
    new_dataframe = pd.DataFrame(new_database)

    # Return
    return new_dataframe


def area_calculation(dataset_func, resolution_func, flood_rate, csv_path):
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
    # Replace data
    newdataset_func = dataset_func.replace([0], -9999)

    # Remove geometry data (x, y coordinates)
    nogeo_data_func = newdataset_func.drop([
        'x_coord', 'y_coord'
    ], axis=1)

    # Get flood rate index
    floodrate_index = get_floodrate_index(flood_rate, csv_path)

    # Filter flood rate
    nogeo_data_func.loc[floodrate_index, :] = -9999

    # Get area dictionary
    area_dict = {}

    # Get area of each simulation
    for each_col in range(len(nogeo_data_func.columns)):
        # Get column num of each simulation
        column = nogeo_data_func[nogeo_data_func.columns[each_col]]

        # Get quantity of flooded cells
        num_flooded_cell = len(column[column != -9999])

        # Get area of each cell
        area_each_cell = resolution_func ** 2

        # Calculate area of each simulation
        area_func = num_flooded_cell * area_each_cell

        # Add each area to area dictionary
        area_dict[f"{nogeo_data_func.columns[each_col]}"] = area_func

    # Write dictionary into dataframe
    area_dataframe = pd.DataFrame(data=area_dict, index=[0])

    return area_dataframe

def area_calculation2(dataset_func, resolution_func, flood_rate, csv_path):
    """
    @Definition:
                A function to calculate area being flooded
    @References:
                None.
    @Arguments:
                dataset_func (pandas dataframe):
                                        A full dataset including unnecessary information. Dataset here is only water
                                        depth.
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
        num_flooded_cell = len(column[column != 0])

        # Get area of each cell
        area_each_cell = resolution_func ** 2

        # Calculate area of each simulation
        area_func = num_flooded_cell * area_each_cell

        # Add each area to area dictionary
        area_dict[f"{nogeo_data_func.columns[each_col]}"] = area_func

    # Write dictionary into dataframe
    area_dataframe = pd.DataFrame(data=area_dict, index=[0])

    return area_dataframe


def building_calculation(dataset_func, building_path, extract_name, onepolygon_bathy_path=None):
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
                extract_name (string):
                                        Name of a specific output among all flood model outputs
                onepolygon_bathy_path (string):
                                        If specify, a string of MAINDIR will be provided to change the version path.
                                        Otherwise, one version will be kept over iteration and causes error file not
                                        existed.
                                        If None, a path will be automatically chosen.
    @Returns:
                building_dataframe (pandas dataframe):
                                        Dataframe of simulations' flooded buildings
    """
    # Choose flood outputs
    if extract_name == 'out.max':
        onepolygon_bathy = wd_onepolygon_bathy
    else:
        onepolygon_bathy = wse_onepolygon_bathy

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

        # Get water map which were converted to one big polygon
        if onepolygon_bathy_path is None:
            watermap_gpf = pd.read_csv(
                fr"{onepolygon_bathy}\\water_polygon_{column}.csv", engine='pyarrow'
            )
        else:
            watermap_gpf = pd.read_csv(
                fr"{onepolygon_bathy_path}\\5_analysis\\wd\\bathy_impact\\onepolygon_nobackground\\water_polygon_{column}.csv", engine='pyarrow'
            )
        watermap_gpf['geometry'] = watermap_gpf['geometry'].apply(wkt.loads)
        watermap_gpf = gpd.GeoDataFrame(watermap_gpf, crs='epsg:2193')

        # Get number of buildings being inundated by merging one big polygon with building polygons
        building_in_polygon = gpd.sjoin(building_gpf, watermap_gpf, predicate='within')

        # Add that number to dictionary of buildings being inundated
        building_dict[f"{column}"] = len(building_in_polygon.index)

    # Write dictionary into dataframe
    building_dataframe = pd.DataFrame(data=building_dict, index=[0])

    return building_dataframe

def error_calculation(pts_extraction_df, validation_df, rmse_calculation=True):
    """
    @Definition:
                A function to calculate RMSEs
    @References:
                None.
    @Arguments:
                pts_extraction_df (pandas dataframe):
                                        A dataframe contains values of observed points
                validation_df (pandas dataframe):
                                        A dataframe contains observed data
    @Returns:
                rmse_dataframe (pandas dataframe):
                                        A dataframe contains RMSEs of each simulation
    """

    # Remove nan
    pts_extraction_df_nonan = pts_extraction_df.dropna()

    # Get indices from dataframe in which nan is removed
    ind = pts_extraction_df_nonan.index.tolist()


    if rmse_calculation:
        # Get empty RMSE dict
        rmse_dict = {}

        # Get RMSE for each simulation
        for i in range(pts_extraction_df_nonan.shape[1]):
            rmse = mean_squared_error(validation_df.iloc[ind, 3], pts_extraction_df_nonan.iloc[:, i], squared=False)
            rmse_dict[f'{pts_extraction_df_nonan.columns[i]}'] = rmse

        # Convert to dataframe
        rmse_dataframe = pd.DataFrame(rmse_dict, index=[0])

        return rmse_dataframe

    else:
        # Get empty MAE dict
        mae_dict = {}

        # Get MAE for each simulation
        for i in range(pts_extraction_df_nonan.shape[1]):
            mae = (validation_df.iloc[ind, 3] - pts_extraction_df_nonan.iloc[:, i]).abs().mean()
            mae_dict[f"{pts_extraction_df_nonan.columns[i]}"] = mae

        # Convert to dataframe
        mae_dataframe = pd.DataFrame(mae_dict, index=[0])

        return mae_dataframe


def calculation_dict(dataset_func, resolution,
                     building_path, flood_rate,
                     extract_name, replace_value,
                     raster_generation_command=True,
                     onepolygon_bathy_path=None,
                     csv_path=None
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
                extract_name (string):
                                        Name of a specific output among all flood model outputs
                replace_value (float):
                                        A value used to replace values lower than flood_rate
                onepolygon_bathy_path (string):
                                        If specify, a string of MAINDIR will be provided to change the version path.
                                        Otherwise, one version will be kept over iteration and causes error file not
                                        existed.
                                        If None, a path will be automatically chosen.
                csv_path (string):
                                        Path link to csv files

    @Returns:
                calculation_dict_set (dictionary):
                                        A dictionary contains all statistical results
    """

    if extract_name == 'elev':

        statistic_method = ['mean', 'sd', 'cv']
        calculation_dict_set = {}

        for each_method in statistic_method:
            # Statistical calculation
            statistic_result = statistic_calculation(extract_name, each_method, flood_rate, replace_value, csv_path)

            # Store in dictionary
            calculation_dict_set[each_method] = statistic_result

            if raster_generation_command:
                raster_generation(
                    statistic_result['x'],
                    statistic_result['y'],
                    statistic_result[each_method],
                    filename=each_method,
                    extract_name='elev'
                )
            else:
                pass

    elif extract_name == 'n':

        statistic_method = ['mean', 'sd', 'cv']
        calculation_dict_set = {}

        for each_method in statistic_method:
            # Statistical calculation
            statistic_result = statistic_calculation(extract_name, each_method, flood_rate, replace_value, csv_path)

            # Store in dictionary
            calculation_dict_set[each_method] = statistic_result

            if raster_generation_command:
                raster_generation(
                    statistic_result['x'],
                    statistic_result['y'],
                    statistic_result[each_method],
                    filename=each_method,
                    extract_name='n'
                )
            else:
                pass


    else:

        statistic_method = ["mean", "sd", "cv", "cell"]

        calculation_dict_set = {}

        for each_method in statistic_method:
            # Statistical calculation
            statistic_result = statistic_calculation(extract_name, each_method, flood_rate, replace_value, csv_path)

            # Store in dictionary
            calculation_dict_set[each_method] = statistic_result

            # If raster_generation_command = True, we will generate rasters of mean, sd, ...
            # If False, we will not generate rasters of mean, sd, ...
            # In comparison modules, we might not need to generate the raster again
            if raster_generation_command:
                if extract_name == 'out.max':
                    # Water depth
                    raster_generation(
                        statistic_result['x'],
                        statistic_result['y'],
                        statistic_result[each_method],
                        filename=each_method,
                        extract_name='out.max',
                    )
                else:
                    # Water surface elevation
                    raster_generation(
                        statistic_result['x'],
                        statistic_result['y'],
                        statistic_result[each_method],
                        filename=each_method,
                        extract_name='out.mxe'
                    )

            else:
                pass

        # Add area result
        calculation_dict_set['area'] = area_calculation2(dataset_func, resolution, flood_rate, csv_path)

        # Add building result
        calculation_dict_set['building'] = building_calculation(dataset_func, building_path, extract_name,
                                                                onepolygon_bathy_path)

    # Add RMSEs
    if onepolygon_bathy_path is None:
        pts_extraction_df = pd.read_csv(fr"{new_wse_csv_bathy}\\rmse_all_simulations.csv")
        validation_df = pd.read_csv(fr"{new_wse_csv_bathy}\\rmse_validation.csv")
    else:
        pts_extraction_df = pd.read_csv(fr"{onepolygon_bathy_path}\\5_analysis\\new_wse\\bathy_csv\\rmse_all_simulations.csv")
        validation_df = pd.read_csv(fr"{onepolygon_bathy_path}\\5_analysis\\new_wse\\bathy_csv\\rmse_validation.csv")
    calculation_dict_set['rmse'] = error_calculation(pts_extraction_df, validation_df)
    calculation_dict_set['mae'] = error_calculation(pts_extraction_df, validation_df, False)


    return calculation_dict_set