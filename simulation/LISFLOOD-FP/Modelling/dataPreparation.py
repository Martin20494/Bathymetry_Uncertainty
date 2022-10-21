# Prepare packages ------------------------------------------------
# Folder packages
from folder import *

# Packages for path creation
import pathlib
import shutil

# Packages for downloading lidar
import geoapis
import geoapis.lidar

# Packages for converting lidar point cloud las/laz file into DEM raster netCDF file
import numpy as np
import shapely.geometry
import geopandas
import json
from geofabrics import processor

# Packages for manipulate raster
import rioxarray as rxr
from osgeo import gdal

import multiprocessing
# -----------------------------------------------------------------------

# NESTED PARALLELISM --------------------------------------------------------------------------------------------------
"""
@Definition:
            A class to perform nested parallelism
@References:
            https://stackoverflow.com/questions/50937362/multiprocessing-on-python-3-jupyter
            https://stackoverflow.com/questions/66420735/simple-multiprocessing-pool-hangs-in-jupyter-notebook
            https://stackoverflow.com/questions/23641475/multiprocessing-working-in-python-but-not-in-ipython/23641560#23641560
            https://stackoverflow.com/questions/5442910/how-to-use-multiprocessing-pool-map-with-multiple-arguments
            https://stackoverflow.com/questions/20886565/using-multiprocessing-process-with-a-maximum-number-of-simultaneous-processes

            https://stackoverflow.com/questions/6974695/python-process-pool-non-daemonic
            https://www.geeksforgeeks.org/partial-functions-python/
"""


class NoDaemonProcess(multiprocessing.Process):
    @property
    def daemon(self):
        return False

    @daemon.setter
    def daemon(self, value):
        pass


class NoDaemonContext(type(multiprocessing.get_context())):
    Process = NoDaemonProcess


# We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
# because the latter is only a wrapper function, not a proper class.
class MyPool(multiprocessing.pool.Pool):
    def __init__(self, *args, **kwargs):
        kwargs['context'] = NoDaemonContext()
        super(MyPool, self).__init__(*args, **kwargs)


# END NESTED PARALLELISM ----------------------------------------------------------------------------------------------


def download_lidar(boundary_coordinates_func, name_dataset, simulation=False):
    """
    @Definition:
                A function to download las/laz files from Open Topography
    @References:
                https://github.com/niwa/geoapis/wiki/Basic-Usage
    @Arguments:
                boundary_coordinates_func (list):
                                                    A list of xmin, ymin, xmax, ymax of rectangular boundary
                name_dataset (string):
                                                    Name of the dataset
                simulation (boolean):
                                                    False if downloading data for preparation
                                                    True if downloading data for simulation
    @Returns:
                None.
    """
    # Get crs
    h_crs = 2193

    # Get xmin, ymin, xmax, ymax of rectangular boundary
    x0_func = boundary_coordinates_func[0]
    y0_func = boundary_coordinates_func[1]
    x1_func = boundary_coordinates_func[2]
    y1_func = boundary_coordinates_func[3]

    # Assign those coordinates into shapely geometry polygon and store under geopandas format
    lidar_bound_coordinates = shapely.geometry.Polygon(
        [(x0_func, y0_func), (x1_func, y0_func), (x1_func, y1_func), (x0_func, y1_func)])
    lidar_bound_coordinates = geopandas.GeoSeries([lidar_bound_coordinates])
    lidar_bound_coordinates = lidar_bound_coordinates.set_crs(h_crs)

    # Create test_catchment.zip
    if simulation:
        test_path = pathlib.Path(f"{bathy_path}\\test_catchment")
    else:
        test_path = pathlib.Path(f"{original_lidar_path}\\test_catchment")
    lidar_bound_coordinates.to_file(test_path)
    shutil.make_archive(base_name=test_path, format='zip', root_dir=test_path)
    shutil.rmtree(test_path)
    lidar_bound_path = pathlib.Path(str(test_path) + ".zip")

    # Create polygon
    lidar_polygon = geopandas.read_file(lidar_bound_path)
    lidar_polygon.to_crs(h_crs)

    # Download lidar
    if simulation:
        lidar_fetcher = geoapis.lidar.OpenTopography(cache_path=fr"{bathy_path}",
                                                     search_polygon=lidar_polygon,
                                                     verbose=True)
    else:
        lidar_fetcher = geoapis.lidar.OpenTopography(cache_path=fr"{original_lidar_path}",
                                                     search_polygon=lidar_polygon,
                                                     verbose=True)
    lidar_fetcher.run(name_dataset)


def dem_raster_reference(resolution_func,
                         chunk_size_func, processor_func,
                         padding_func, lidar_dataset_name, element_name):
    """
    @Definition:
                A function to create a reference raster file from a las/laz file used for center calculation and
                padding reference
    @References:
                https://github.com/rosepearson/GeoFabrics
    @Arguments:
                resolution_func (int or float):
                                            Resolution value in meter
                chunk_size_func (int):
                                            Size value of chunk
                processor_func (int):
                                            Number of processor
                padding_func (list):
                                            A list of x min, x max, y min and y max
                lidar_dataset_name (string):
                                            LiDAR name
                element_name (string):
                                            Name of element folder
    @Returns:
                None.
    """
    # Set crs
    h_crs = 2193
    v_crs = 7839

    # Paths
    basepath = pathlib.Path(original_lidar_path)
    result_folder = element_name
    data_dir = basepath / result_folder
    data_dir.mkdir(parents=True, exist_ok=True)

    # Bounding box/ Catchment boundary
    x0 = padding_func[0]
    y0 = padding_func[1]
    x1 = padding_func[2]
    y1 = padding_func[3]
    catchment = shapely.geometry.Polygon([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])
    catchment = geopandas.GeoSeries([catchment])
    catchment = catchment.set_crs(h_crs)
    catchment.to_file(basepath / result_folder / "catchment_boundary.geojson", crs="EPSG:2193", driver="GeoJSON")

    # Design JSON instructions
    instruction_json = {}

    # DRAINS
    instruction_json['drains'] = {
        # drain - output
        "output": {
            "crs": {
                "horizontal": h_crs,
                "vertical": v_crs
            },
            "grid_params": {
                "resolution": 1
            }
        },

        # drain - processing
        "processing": {
            "chunk_size": 1000,
            "number_of_cores": 4
        },

        # drain - data_paths
        "data_paths": {
            "local_cache": str(basepath),
            "subfolder": result_folder,
            "catchment_boundary": "catchment_boundary.geojson"
        },

        # drain - apis
        "apis": {
            "open_topography": {
                "Wellington_2013": {
                    "crs": {
                        "horizontal": h_crs,
                        "vertical": v_crs
                    }
                }
            },
            "linz": {
                "key": "857413f41ce446ed8961e2f1e960a24b",
                "land": {
                    "layers": [51559]
                }
            }
        },

        # drain - general
        "general": {
            "lidar_classifications_to_keep": [2, 9]
        },

        # drain - drains
        "drains": {
            "width": 5
        }
    }

    # RIVERS
    instruction_json['rivers'] = {
        # river - output
        "output": {
            "crs": {
                "horizontal": h_crs,
                "vertical": v_crs
            },
            "grid_params": {
                "resolution": 1
            }
        },

        # river - processing
        "processing": {
            "chunk_size": 1000,
            "number_of_cores": 4
        },

        # river - data_paths
        "data_paths": {
            "local_cache": str(basepath),
            "subfolder": result_folder,
        },

        # river - apis
        "apis": {
            "open_topography": {
                "Wellington_2013": {
                    "crs": {
                        "horizontal": h_crs,
                        "vertical": v_crs
                    }
                }
            },
            "linz": {
                "key": "857413f41ce446ed8961e2f1e960a24b",
                "land": {
                    "layers": [51559]
                },
                "bathymetry_contours": {
                    "layers": [50554]
                }
            }
        },

        # river - general
        "general": {
            "set_dem_shoreline": True,
            "bathymetry_contours_z_label": "valdco",
            "drop_offshore_lidar": False,
            "lidar_classifications_to_keep": [2, 9],
            "interpolate_missing_values": True
        },

        # river - rivers
        "rivers": {
            "osm_id": 132793862,
            "veg_lidar_classifications_to_keep": [2, 3, 4, 5, 9],
            "max_channel_width": 120,
            "min_channel_width": 10,
            "max_bank_height": 2,
            "rec_alignment_tolerance": 65,
            "width_centre_smoothing": 10,
            "channel_area_threshold": 125000000, #125000000 or 141000000
            "channel_rec_id": 9253579,
            "cross_section_spacing": 10,
            "min_bank_height": 0.75,
            "rec_file": str(basepath / "rec2_3.geojson"),
            "flow_file": str(basepath / "flow_and_friction.csv.gz")
        }
    }

    # DEM
    instruction_json['dem'] = {
        # outputs
        "output": {
            "crs": {
                "horizontal": h_crs,
                "vertical": v_crs
            },
            "grid_params": {
                "resolution": resolution_func
            }
        },

        # processing
        "processing": {
            "chunk_size": chunk_size_func,
            "number_of_cores": processor_func
        },

        # apis
        "apis": {
            "open_topography": {
                f"{lidar_dataset_name}": {
                    "crs": {
                        "horizontal": h_crs,
                        "vertical": v_crs
                    }
                }
            },
            "linz": {
                "key": "857413f41ce446ed8961e2f1e960a24b",
                "land": {
                    "layers": [51559]
                },
                "bathymetry_contours": {
                    "layers": [50554]
                }
            }
        },

        # data paths
        "data_paths": {
            "local_cache": str(basepath),
            "subfolder": result_folder,
            "catchment_boundary": "catchment_boundary.geojson",
            "raw_dem": "raw_dem.nc",

            # For bathymetry
            "river_bathymetry": ["river_bathymetry.geojson", "fan_bathymetry.geojson",
                                 "open_drain_elevation_5m_width.geojson", "closed_drain_elevation_5m_width.geojson"],
            "river_polygons": ["river_polygon.geojson", "fan_polygon.geojson",
                               "open_drain_polygon_5m_width.geojson", "closed_drain_polygon_5m_width.geojson"],
            # Result
            "result_dem": f"{element_name}.nc"
        },

        # general
        "general": {
            # For lidar
            "set_dem_shoreline": True,
            "drop_offshore_lidar": False,
            "lidar_classification_to_keep": [2, 9],
            "interpolation_method": "nearest",
            "lidar_interpolation_method": "idw",

            # For bathymetry
            "bathymetry_contours_z_label": "valdco",
            "bathymetry_points_type": ['rivers', 'rivers', 'drains', 'drains'],
            "bathymetry_points_z_label": ['bed_elevation_Rupp_and_Smart', "depths", "elevation", "elevation"]
        }
    }

    # Save the instructions
    with open(basepath / result_folder / "instructions.json", "w") as instruction:
        json.dump(instruction_json, instruction)

    # Create DEM raster
    runner = processor.RiverBathymetryGenerator(instruction_json['rivers'])
    runner.run()
    runner = processor.DrainBathymetryGenerator(instruction_json['drains'])
    runner.run()
    runner = processor.RawLidarDemGenerator(instruction_json['dem'])
    runner.run()
    runner = processor.HydrologicDemGenerator(instruction_json['dem'])
    runner.run()

def terrain_shading(
    altitude,
    azimuth
):
    """
    @Definition:
                A function to create terrain shading
    @References:
                https://www.youtube.com/watch?v=5dDZeEXws9Q
                https://blog.datawrapper.de/shaded-relief-with-gdal-python/
                https://www.geophysique.be/2014/02/25/shaded-relief-map-in-python/
                https://www.l3harrisgeospatial.com/docs/topographicshading.html
    @Arguments:
                azimuth and altitude (float):
                            Values to customise the terrain shading. Ex: altitude = 45, azimuth = 355
    @Returns:
                None.
    """
    # Convert nc to tiff
    terrain = rxr.open_rasterio(fr"{original_lidar_path}\\shading\\shading.nc")
    terrain.z.rio.to_raster(fr"{original_lidar_path}\\shading\\shading.tiff")

    # Create shading
    terrain_tiff = gdal.Open(fr"{original_lidar_path}\\shading\\shading.tiff")
    terrain_shade = gdal.DEMProcessing(
        fr"{original_lidar_path}\\shading\\terrain_shading.tiff",
        terrain_tiff,
        "hillshade",
        computeEdges=True,
        altitude=altitude, azimuth=azimuth
    )

    # Close the file
    terrain_shade = None
    terrain_tiff = None
