# Prepare packages ----------------------------------------------------------------------------------------------------
from randomNumber import uniform_dist_number, normal_dist_number
from depth import *
from copyFile import copy_geom
# ---------------------------------------------------------------------------------------------------------------------

# GET ORIGINAL WATER SURFACE ELEVATION --------------------------------------------------------------------------------
def valueLIDAR_list_generation(kind):
    """
    @Definition:
                A function to generate estimated water surface elevation
    @References:
                None.
    @Arguments:
                kind (string):
                                            Kind of depth calculation ('neal' or 'rupp')
    @Returns:
                wse (numpy array):
                                            Estimated water surface elevation from Rose's codes
    """

    # Get empty list of depth
    depth_list = []

    # Get each row information from river characteristics and calculate the depth
    for i in range(river_characteristics.shape[0]):
        row = river_characteristics.iloc[i]
        if kind == 'neal':
            depth = calculate_neal_et_al(row[slope_name], row[width_name], row[flow_name], row[friction_name])

        else:
            depth = calculate_rupp_et_al(row[slope_name], row[width_name], row[flow_name])
        depth_list.append(depth)

    wse = river_bathymetry.bed_elevation_Neal_et_al.to_numpy() + np.array(depth_list).astype('float64')

    return wse

# Get wse from Neal
wse_neal = valueLIDAR_list_generation("neal")
# END GET ORIGINAL WATER SURFACE ELEVATION ----------------------------------------------------------------------------


# UNIFORM DISTRIBUTION ##################################################################################################
# PARAMETERS ----------------------------------------------------------------------------------------------------------
def generate_para_uniform(
    dataset_id_number,
    chosen_para, random_lower_upper
):
    """
    @Definition:
                A function to generate random value of chosen parameter
    @References:
                None.
    @Arguments:
                dataset_id_number (int):
                                ID of random uniform distribution
                chosen_para (string):
                                Chosen para to be varied (width, slope, friction, and flow)
                random_lower_upper (list)
                                Upper and lower limits of random uniform distribution
    @Returns:
                para_arr (numpy array):
                                Array of random value of chosen parameters being varied - cross-sections of a river
                depth_neal_arr (numpy array):
                                Array of neal depth from chosen parameter being varied - cross-sections of a river
                bed_neal_arr (numpy array):
                                Array of neal bed elevation from chosen parameter being varied - cross-sections of a river
    """
    # Get empty list
    para_list = []
    depth_neal_list = []

    # Extract a number for each cross-section from normal distributions
    for i in range(river_characteristics.shape[0]):
        # Get each row in river data
        row = river_characteristics.iloc[i]
        # Get chosen addition value
        selected_addition = uniform_dist_number(random_lower_upper[0], random_lower_upper[1], dataset_id_number)

        # Vary para
        if chosen_para == 'width':
            # Width
            selected_para = row[width_name] + selected_addition
            para_list.append(selected_para)

            # Get depth
            depth_neal = calculate_neal_et_al(row[slope_name], selected_para, row[flow_name], row[friction_name])
            depth_neal_list.append(depth_neal)

        elif chosen_para == 'slope':
            # Slope
            selected_para = row[slope_name] + selected_addition
            para_list.append(selected_para)

            # Get depth
            depth_neal = calculate_neal_et_al(selected_para, row[width_name], row[flow_name], row[friction_name])
            depth_neal_list.append(depth_neal)

        elif chosen_para == 'friction':
            # Friction
            selected_para = row[friction_name] + selected_addition
            para_list.append(selected_para)

            # Get depth
            depth_neal = calculate_neal_et_al(row[slope_name], row[width_name], row[flow_name], selected_para)
            depth_neal_list.append(depth_neal)

        elif chosen_para == 'flow':
            # Flow
            selected_para = row[flow_name] + selected_addition
            para_list.append(selected_para)

            # Get depth
            depth_neal = calculate_neal_et_al(row[slope_name], row[width_name], selected_para, row[friction_name])
            depth_neal_list.append(depth_neal)

    # Convert list to array
    para_arr = np.array(para_list).astype('float64')
    depth_neal_arr = np.array(depth_neal_list).astype('float64')

    # Calculate bed elev
    bed_neal_arr = wse_neal - depth_neal_arr

    return para_arr, depth_neal_arr, bed_neal_arr

def generate_para_uniform_df(
    dataset_id_size_list,
    chosen_para, random_lower_upper
):
    """
    @Definition:
                A function to generate dataframe of number of simulation of chosen parameters
    @References:
                None.
    @Arguments:
                sample_size (int):
                                Size of random uniform distribution
                dataset_id_number (int):
                                ID of random uniform distribution
                chosen_para (string):
                                Chosen para to be varied (width, slope, friction, and flow)
                random_lower_upper (list)
                                Upper and lower limits of random uniform distribution
    @Returns:
                para_df (dataframe):
                                Dataframe of number of varied parameters
                depth_df (dataframe):
                                Dataframe of number of depths when varying parameters
                bed_df (dataframe):
                                Dataframe of number of bed elevations when varying parameters

    """
    # Create geopandas dataframe containing geomtry
    para_df = gpd.GeoDataFrame(data=river_characteristics.geometry)
    depth_df = gpd.GeoDataFrame(data=river_characteristics.geometry)
    bed_df = gpd.GeoDataFrame(data=river_characteristics.geometry)

    # Loop for getting values
    for ids in dataset_id_size_list:

        # Get para and depth values for dataframes
        para_arr, depth_neal_arr, bed_neal_arr = generate_para_uniform(ids, chosen_para, random_lower_upper)

        # Append para values
        para_df[f'id_{ids}'] = para_arr

        # Append depth values
        depth_df[f'id_{ids}'] = depth_neal_arr

        # Append bed values
        bed_df[f'id_{ids}'] = bed_neal_arr

        # Copy other files
        copy_geom([51559, 50554], ids)

        # Add simulated bed values into river bathymetry
        copy_river_bathymetry = river_bathymetry.copy(deep=True)
        copy_river_bathymetry['bed_elevation_Neal_et_al'] = bed_neal_arr
        copy_river_bathymetry.to_file(fr"{bathy_path}\\bathy_{ids}\\river_bathymetry.geojson")

    # Write out dataframes
    para_df.to_file(fr"{bathy_path}\\para_{chosen_para}_df.geojson")
    depth_df.to_file(fr"{bathy_path}\\depth_{chosen_para}_df.geojson")
    bed_df.to_file(fr"{bathy_path}\\bed_{chosen_para}_df.geojson")

# END PARAMETERS -----------------------------------------------------------------------------------------------------


# EXPONENTS ----------------------------------------------------------------------------------------------------------
def generate_expo_uniform(
    expo, expo_id,
    lower, upper
):
    """
    @Definition:
                A function to generate random value of chosen exponent
    @References:
                None.
    @Arguments:
                expo (string):
                            'alpha' or 'beta'
                expo_id (int):
                            ID of chosen exponent
                lower and upper (float):
                            lower and upper limits of random uniform distribution
    @Returns:
                expo_number (float):
                            Random exponent from uniform distribution
                bed_neal_arr (numpy array):
                            Array of Neal bed elevation when varying chosen exponent
    """
    # Get empty list
    depth_neal_list = []

    # Get exponent
    expo_number = uniform_dist_number(
        lower, upper, expo_id
    )

    # Extract a number for each cross-section from normal distrbution
    for i in range(river_characteristics.shape[0]):
        # Get each row in river data
        row = river_characteristics.iloc[i]

        # Calculate depth
        depth = calculate_neal_et_al_expo(
            expo, expo_number,
            row[slope_name], row[width_name], row[flow_name], row[friction_name]
        )
        depth_neal_list.append(depth)

    # Convert list to array
    depth_neal_arr = np.array(depth_neal_list).astype('float64')

    # Calculate bed elev
    bed_neal_arr = wse_neal - depth_neal_arr

    return expo_number, depth_neal_arr, bed_neal_arr


def generate_expo_uniform_df(
    expo, expo_id_set,
    lower, upper
):
    """
    @Definition:
                A function to generate dataframes of expo, depth, and bed
    @References:
                None.
    @Arguments:
                expo (string):
                            'alpha' or 'beta'
                expo_id_set (int):
                            sets of ID of chosen exponent
                lower and upper (float):
                            lower and upper limits of random uniform distribution
    @Returns:
                para_df (dataframe):
                                Dataframe of number of varied exponent
                depth_df (dataframe):
                                Dataframe of number of depths when varying exponent
                bed_df (dataframe):
                                Dataframe of number of bed elevations when varying exponent
    """
    # Generate exponent, depth, and bed dataframe
    expo_dict = {}
    depth_df = gpd.GeoDataFrame(data=river_characteristics.geometry)
    bed_df = gpd.GeoDataFrame(data=river_characteristics.geometry)

    # Loop to create dataframe
    for expo_id in range(expo_id_set):

        # Get exponent and depth values for dataframes
        expo_number, depth_arr, bed_arr = generate_expo_uniform(
            expo, expo_id, lower, upper
        )

        # Get exponent list
        expo_dict[f'id_{expo_id}'] = expo_number

        # Append depth values
        depth_df[f'id_{expo_id}'] = depth_arr

        # Append bed values
        bed_df[f'id_{expo_id}'] = bed_arr

        # Copy other files
        copy_geom([51559, 50554], expo_id)

        # Add simulated bed values into river bathymetry
        copy_river_bathymetry = river_bathymetry.copy(deep=True)
        copy_river_bathymetry['bed_elevation_Neal_et_al'] = bed_arr
        copy_river_bathymetry.to_file(fr"{bathy_path}\\bathy_{expo_id}\\river_bathymetry.geojson")


    # Write out dataframes
    depth_df.to_file(fr"{bathy_path}\\depth_{expo}_df.geojson")
    bed_df.to_file(fr"{bathy_path}\\bed_{expo}_df.geojson")

    # Get exponent dataframe
    expo_df = pd.DataFrame(data=expo_dict, index=[0])
    expo_df.to_file(fr"{bathy_path}\\expo_{expo}_df.geojson")

# END EXPONENTS ------------------------------------------------------------------------------------------------------
# END UNIFORM DISTRIBUTION ###########################################################################################



# NORMAL DISTRIBUTION ################################################################################################
def generate_para_normal(
    dataset_id_number,
    chosen_para,
    random_sd, sample_size
):
    """
    @Definition:
                A function to generate random value of chosen parameter
    @References:
                None.
    @Arguments:
                dataset_id_number (int):
                                ID of random uniform distribution
                chosen_para (string):
                                Chosen para to be varied (width, slope, friction, and flow)
                random_lower_upper (list)
                                Upper and lower limits of random uniform distribution
    @Returns:
                para_arr (numpy array):
                                Array of random value of chosen parameters being varied - cross-sections of a river
                depth_neal_arr (numpy array):
                                Array of neal depth from chosen parameter being varied - cross-sections of a river
                bed_neal_arr (numpy array):
                                Array of neal bed elevation from chosen parameter being varied - cross-sections of a river
    """
    # Get empty list
    para_list = []
    depth_neal_list = []

    # Extract a number for each cross-section from normal distributions
    for i in range(river_characteristics.shape[0]):
        # Get each row in river data
        row = river_characteristics.iloc[0]

        # Vary para
        if chosen_para == 'width':
            # Width
            selected_para = normal_dist_number(row[width_name], random_sd, sample_size, dataset_id_number, i)
            para_list.append(selected_para)

            # Get depth
            depth_neal = calculate_neal_et_al(row[slope_name], selected_para, row[flow_name], row[friction_name])
            depth_neal_list.append(depth_neal)

        elif chosen_para == 'slope':
            # Slope
            selected_para = normal_dist_number(row[slope_name], random_sd, sample_size, dataset_id_number, i)
            para_list.append(selected_para)

            # Get depth
            depth_neal = calculate_neal_et_al(selected_para, row[width_name], row[flow_name], row[friction_name])
            depth_neal_list.append(depth_neal)

        elif chosen_para == 'friction':
            # Friction
            selected_para = normal_dist_number(row[friction_name], random_sd, sample_size, dataset_id_number, i)
            para_list.append(selected_para)

            # Get depth
            depth_neal = calculate_neal_et_al(row[slope_name], row[width_name], row[flow_name], selected_para)
            depth_neal_list.append(depth_neal)

        elif chosen_para == 'flow':
            # Flow
            selected_para = normal_dist_number(row[flow_name], random_sd, sample_size, dataset_id_number,
                                               int(row[flow_name]))
            para_list.append(selected_para)

            # Get depth
            depth_neal = calculate_neal_et_al(row[slope_name], row[width_name], selected_para, row[friction_name])
            depth_neal_list.append(depth_neal)

    # Convert list to array
    para_arr = np.array(para_list).astype('float64')
    depth_neal_arr = np.array(depth_neal_list).astype('float64')

    # Calculate bed elev
    bed_neal_arr = wse_neal - depth_neal_arr

    return para_arr, depth_neal_arr, bed_neal_arr

def generate_para_normal_df(
    chosen_para,
    random_sd, sample_size,
    dataset_id_size_list
):
    """

    """
    # Create geopandas dataframe containing geomtry
    para_df = gpd.GeoDataFrame(data=river_characteristics.geometry)
    depth_df = gpd.GeoDataFrame(data=river_characteristics.geometry)
    bed_df = gpd.GeoDataFrame(data=river_characteristics.geometry)

    # Loop for getting values
    for ids in dataset_id_size_list:
        # Get para and depth values for dataframes
        para_arr, depth_neal_arr, bed_neal_arr = generate_para_normal(
            ids, chosen_para, random_sd, sample_size)

        # Append para values
        para_df[f'id_{ids}'] = para_arr

        # Append depth values
        depth_df[f'id_{ids}'] = depth_neal_arr

        # Append bed values
        bed_df[f'id_{ids}'] = bed_neal_arr

        # Copy other files
        copy_geom([51559, 50554], ids)

        # Add simulated bed values into river bathymetry
        copy_river_bathymetry = river_bathymetry.copy(deep=True)
        copy_river_bathymetry['bed_elevation_Neal_et_al'] = bed_neal_arr
        copy_river_bathymetry.to_file(fr"{bathy_path}\\bathy_{ids}\\river_bathymetry.geojson")

    # Write out dataframes
    para_df.to_file(fr"{bathy_path}\\para_{chosen_para}_df.geojson")
    depth_df.to_file(fr"{bathy_path}\\depth_{chosen_para}_df.geojson")
    bed_df.to_file(fr"{bathy_path}\\bed_{chosen_para}_df.geojson")

# END NORMAL DISTRIBUTION ############################################################################################













