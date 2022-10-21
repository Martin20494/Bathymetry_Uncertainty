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

    if kind == 'neal':
        wse = river_bathymetry.bed_elevation_Neal_et_al.to_numpy() + np.array(depth_list).astype('float64')
    else:
        wse = river_bathymetry.bed_elevation_Rupp_and_Smart.to_numpy() + np.array(depth_list).astype('float64')

    return wse

# Get wse from Neal
wse_rupp = valueLIDAR_list_generation("rupp")
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
    depth_rupp_list = []

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
            depth_rupp = calculate_rupp_et_al(row[slope_name], selected_para, row[flow_name])
            depth_rupp_list.append(depth_rupp)

        elif chosen_para == 'slope':
            # Slope
            selected_para = row[slope_name] + selected_addition
            para_list.append(selected_para)

            # Get depth
            depth_rupp = calculate_rupp_et_al(selected_para, row[width_name], row[flow_name])
            depth_rupp_list.append(depth_rupp)

        elif chosen_para == 'flow':
            # Flow
            selected_para = row[flow_name] + selected_addition
            para_list.append(selected_para)

            # Get depth
            depth_rupp = calculate_rupp_et_al(row[slope_name], row[width_name], selected_para)
            depth_rupp_list.append(depth_rupp)

    # Convert list to array
    para_arr = np.array(para_list).astype('float64')
    depth_rupp_arr = np.array(depth_rupp_list).astype('float64')

    # Calculate bed elev
    bed_rupp_arr = wse_rupp - depth_rupp_arr

    return para_arr, depth_rupp_arr, bed_rupp_arr

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
        para_arr, depth_rupp_arr, bed_rupp_arr = generate_para_uniform(ids, chosen_para, random_lower_upper)

        # Append para values
        para_df[f'id_{ids}'] = para_arr

        # Append depth values
        depth_df[f'id_{ids}'] = depth_rupp_arr

        # Append bed values
        bed_df[f'id_{ids}'] = bed_rupp_arr

        # Copy other files
        copy_geom([51559, 50554], ids)

        # Add simulated bed values into river bathymetry
        copy_river_bathymetry = river_bathymetry.copy(deep=True)
        copy_river_bathymetry['bed_elevation_Rupp_and_Smart'] = bed_rupp_arr
        copy_river_bathymetry.to_file(fr"{bathy_path}\\bathy_{ids}\\river_bathymetry.geojson")

    # Write out dataframes
    para_df.to_file(fr"{bathy_path}\\para_{chosen_para}_df.geojson")
    depth_df.to_file(fr"{bathy_path}\\depth_{chosen_para}_df.geojson")
    bed_df.to_file(fr"{bathy_path}\\bed_{chosen_para}_df.geojson")

# END PARAMETERS -----------------------------------------------------------------------------------------------------


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
                depth_rupp_arr (numpy array):
                                Array of neal depth from chosen parameter being varied - cross-sections of a river
                bed_rupp_arr (numpy array):
                                Array of neal bed elevation from chosen parameter being varied - cross-sections of a river
    """
    # Get empty list
    para_list = []
    depth_rupp_list = []

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
            depth_rupp = calculate_rupp_et_al(row[slope_name], selected_para, row[flow_name])
            depth_rupp_list.append(depth_rupp)

        elif chosen_para == 'slope':
            # Slope
            selected_para = normal_dist_number(row[slope_name], random_sd, sample_size, dataset_id_number, i)
            para_list.append(selected_para)

            # Get depth
            depth_rupp = calculate_rupp_et_al(selected_para, row[width_name], row[flow_name])
            depth_rupp_list.append(depth_rupp)

        elif chosen_para == 'flow':
            # Flow
            selected_para = normal_dist_number(row[flow_name], random_sd, sample_size, dataset_id_number, i)
            para_list.append(selected_para)

            # Get depth
            depth_rupp = calculate_rupp_et_al(row[slope_name], row[width_name], selected_para)
            depth_rupp_list.append(depth_rupp)

    # Convert list to array
    para_arr = np.array(para_list).astype('float64')
    depth_rupp_arr = np.array(depth_rupp_list).astype('float64')

    # Calculate bed elev
    bed_rupp_arr = wse_rupp - depth_rupp_arr

    return para_arr, depth_rupp_arr, bed_rupp_arr

def generate_para_normal_df(
    chosen_para,
    random_sd, sample_size,
    dataset_id_size_list
):
    """
    @Definition:
                A function to generate dataframe of random value of chosen parameter
    @References:
                None.
    @Arguments:

    @Returns:

    """
    # Create geopandas dataframe containing geomtry
    para_df = gpd.GeoDataFrame(data=river_characteristics.geometry)
    depth_df = gpd.GeoDataFrame(data=river_characteristics.geometry)
    bed_df = gpd.GeoDataFrame(data=river_characteristics.geometry)

    # Loop for getting values
    for ids in dataset_id_size_list:
        # Get para and depth values for dataframes
        para_arr, depth_rupp_arr, bed_rupp_arr = generate_para_normal(
            ids, chosen_para, random_sd, sample_size)

        # Append para values
        para_df[f'id_{ids}'] = para_arr

        # Append depth values
        depth_df[f'id_{ids}'] = depth_rupp_arr

        # Append bed values
        bed_df[f'id_{ids}'] = bed_rupp_arr

        # Copy other files
        copy_geom([51559, 50554], ids)

        # Add simulated bed values into river bathymetry
        copy_river_bathymetry = river_bathymetry.copy(deep=True)
        copy_river_bathymetry['bed_elevation_Rupp_and_Smart'] = bed_rupp_arr
        copy_river_bathymetry.to_file(fr"{bathy_path}\\bathy_{ids}\\river_bathymetry.geojson")

    # Write out dataframes
    para_df.to_file(fr"{bathy_path}\\para_{chosen_para}_df.geojson")
    depth_df.to_file(fr"{bathy_path}\\depth_{chosen_para}_df.geojson")
    bed_df.to_file(fr"{bathy_path}\\bed_{chosen_para}_df.geojson")

# END NORMAL DISTRIBUTION ############################################################################################


# FULL NORMAL DISTRIBUTION ###########################################################################################
def generate_para_normal_full(
    dataset_id_number,
    sd_slope, sd_width, sd_flow,
    sample_size
):
    """
    @Definition:
                A function to generate random value of chosen parameter
    @References:
                None.
    @Arguments:
                dataset_id_number (int):
                                ID of random uniform distribution
                sd_slope, sd_width, sd_flow (float):
                                Random standard deviation for slope, width, and flow
                sample_size (int):
                                Sample size for normal distribution
    @Returns:
                slope_arr, width_arr, flow_arr (numpy array):
                                Array of random value of chosen parameters being varied - cross-sections of a river
                depth_rupp_arr (numpy array):
                                Array of neal depth from chosen parameter being varied - cross-sections of a river
                bed_rupp_arr (numpy array):
                                Array of neal bed elevation from chosen parameter being varied - cross-sections of a river
    """
    # Get empty list
    width_list = []
    slope_list = []
    flow_list = []

    depth_rupp_list = []

    # Extract a number for each cross-section from normal distributions
    for i in range(river_characteristics.shape[0]):
        # Get each row in river data
        row = river_characteristics.iloc[0]

        # Para ---------------------------------------------------------
        # Slope
        slope_para = normal_dist_number(row[slope_name], sd_slope, sample_size, dataset_id_number, i)
        slope_list.append(slope_para)

        # Width
        width_para = normal_dist_number(row[width_name], sd_width, sample_size, dataset_id_number, i)
        width_list.append(width_para)

        # Flow
        flow_para = normal_dist_number(row[flow_name], sd_flow, sample_size, dataset_id_number, i)
        flow_list.append(flow_para)


        # Depth --------------------------------------------------------
        depth_rupp = calculate_rupp_et_al(slope_para, width_para, flow_para)
        depth_rupp_list.append(depth_rupp)

    # Convert list to array
    # Para ----------------
    slope_arr = np.array(slope_list).astype('float64')
    width_arr = np.array(width_list).astype('float64')
    flow_arr = np.array(flow_list).astype('float64')
    # Depth --------------
    depth_rupp_arr = np.array(depth_rupp_list).astype('float64')

    # Calculate bed elev
    bed_rupp_arr = wse_rupp - depth_rupp_arr

    return slope_arr, width_arr, flow_arr, depth_rupp_arr, bed_rupp_arr

def generate_para_normal_full_df(
    sd_slope, sd_width, sd_flow,
    sample_size,
    dataset_id_size_list
):
    """
    @Definition:
                A function to generate dataframe of random value of chosen parameter
    @References:
                None.
    @Arguments:

    @Returns:

    """
    # Create geopandas dataframe containing geometry
    # Para ----------------------------------------
    slope_df = gpd.GeoDataFrame(data=river_characteristics.geometry)
    width_df = gpd.GeoDataFrame(data=river_characteristics.geometry)
    flow_df = gpd.GeoDataFrame(data=river_characteristics.geometry)
    # Depth and bed ------------------------------
    depth_df = gpd.GeoDataFrame(data=river_characteristics.geometry)
    bed_df = gpd.GeoDataFrame(data=river_characteristics.geometry)

    # Loop for getting values
    for ids in dataset_id_size_list:
        # Get para and depth values for dataframes
        slope_arr, width_arr, flow_arr, depth_rupp_arr, bed_rupp_arr = generate_para_normal_full(
            ids,
            sd_slope, sd_width, sd_flow,
            sample_size)

        # Append para values
        slope_df[f'id_{ids}'] = slope_arr
        width_df[f'id_{ids}'] = width_arr
        flow_df[f'id_{ids}'] = flow_arr

        # Append depth values
        depth_df[f'id_{ids}'] = depth_rupp_arr

        # Append bed values
        bed_df[f'id_{ids}'] = bed_rupp_arr

        # Copy other files
        copy_geom([51559, 50554], ids)

        # Add simulated bed values into river bathymetry
        copy_river_bathymetry = river_bathymetry.copy(deep=True)
        copy_river_bathymetry['bed_elevation_Rupp_and_Smart'] = bed_rupp_arr
        copy_river_bathymetry.to_file(fr"{bathy_path}\\bathy_{ids}\\river_bathymetry.geojson")

    # Write out dataframes
    slope_df.to_file(fr"{bathy_path}\\para_slope_df.geojson")
    width_df.to_file(fr"{bathy_path}\\para_width_df.geojson")
    flow_df.to_file(fr"{bathy_path}\\para_width_df.geojson")
    depth_df.to_file(fr"{bathy_path}\\depth_df.geojson")
    bed_df.to_file(fr"{bathy_path}\\bed_df.geojson")


# END FULL NORMAL DISTRIBUTION #######################################################################################











