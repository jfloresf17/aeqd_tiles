import pathlib
from typing import Dict, List, Optional, Union

from pydantic import BaseModel

class SensorProduct(BaseModel):
    snippet: Union[Dict, str]
    bands: Optional[List[str]]
    resolution: float
    start_date: Optional[str]
    end_date: Optional[str]

class Sensors(BaseModel):
    sensors: Dict[str, SensorProduct]


def get_sensor_products(
    product: Optional[List[str]] = [
        "MODIS_BRDF",
        "WorldCLIM",
        "MODIS_LAND_SURFACE",
        "CCNL",
        "fabdem",
        "geomorpho90m",
        "s1gbm",
        "isric",
        "soil_salinity",
        "ai0",
        "gwa",
        "WorldCover",
        "OpenLandMap",
        "malaria",
        "habitat",
        "GLO_DEM",
    ]
) -> Sensors:
    # Create a dictionary of sensors
    sensors = dict()

    # Add the sensors to the dictionary
    if "MODIS_BRDF" in product:
        sensors["MODIS_BRDF"] = {
            "snippet": "MODIS/061/MCD43A4",
            "bands": [
                "Nadir_Reflectance_Band1",
                "Nadir_Reflectance_Band2",
                "Nadir_Reflectance_Band3",
                "Nadir_Reflectance_Band4",
                "Nadir_Reflectance_Band5",
                "Nadir_Reflectance_Band6",
                "Nadir_Reflectance_Band7",
                "BRDF_Albedo_Band_Mandatory_Quality_Band1",
                "BRDF_Albedo_Band_Mandatory_Quality_Band2",
                "BRDF_Albedo_Band_Mandatory_Quality_Band3",
                "BRDF_Albedo_Band_Mandatory_Quality_Band4",
                "BRDF_Albedo_Band_Mandatory_Quality_Band5",
                "BRDF_Albedo_Band_Mandatory_Quality_Band6",
                "BRDF_Albedo_Band_Mandatory_Quality_Band7",
            ],
            "resolution": 512,
            "start_date": "2000-02-24",
            "end_date": None,
        }

    if "WorldCLIM" in product:
        sensors["WorldCLIM"] = {
            "snippet": "WORLDCLIM/V1/MONTHLY",
            "bands": ["tavg", "tmin", "tmax", "prec"],
            "resolution": 1024,
            "start_date": "1960-01-01",
            "end_date": "1990-01-01",
        }

    if "MODIS_LAND_SURFACE" in product:
        sensors["MODIS_LAND_SURFACE"] = {
            "snippet": "MODIS/061/MOD11A1",
            "bands": ["LST_Day_1km", "LST_Night_1km", "QC_Day", "QC_Night"],
            "resolution": 1024,
            "start_date": "2000-02-24",
            "end_date": None,
        }

    if "CCNL" in product:
        sensors["CCNL"] = {
            "snippet": "BNU/FGS/CCNL/v1",
            "bands": ["b1"],
            "resolution": 1024,
            "start_date": "1992-01-01",
            "end_date": "2014-01-01",
        }

    if "fabdem" in product:
        sensors["fabdem"] = {
            "snippet": "projects/sat-io/open-datasets/FABDEM",
            "bands": [],
            "resolution": 128,
            "start_date": None,
            "end_date": None,
        }

    if "geomorpho90m" in product:
        sensors["geomorpho90m"] = {
            "snippet": "projects/sat-io/open-datasets/Geomorpho90m/geom",
            "bands": [],
            "resolution": 256,
            "start_date": None,
            "end_date": None,
        }

    if "s1gbm" in product:
        sensors["s1gbm"] = {
            "snippet": {
                "VH": "projects/sat-io/open-datasets/S1GBM/normalized_s1_backscatter_VH",
                "VV": "projects/sat-io/open-datasets/S1GBM/normalized_s1_backscatter_VV",
            },
            "bands": [],
            "resolution": 128,
            "start_date": None,
            "end_date": None,
        }

    if "isric" in product:
        sensors["isric"] = {
            "snippet": {
                key: f"projects/soilgrids-isric/{key}"
                for key in [
                    "bdod_mean",
                    "cec_mean",
                    "cfvo_mean",
                    "clay_mean",
                    "sand_mean",
                    "silt_mean",
                    "nitrogen_mean",
                    "phh2o_mean",
                    "soc_mean",
                    "ocd_mean",
                    "ocs_mean",
                ]
            },
            "bands": [],
            "resolution": 256,
            "start_date": None,
            "end_date": None,
        }

    if "soil_salinity" in product:
        sensors["soil_salinity"] = {
            "snippet": "projects/sat-io/open-datasets/global_soil_salinity",
            "bands": [],
            "resolution": 256,
            "start_date": None,
            "end_date": None,
        }

    if "ai0" in product:
        sensors["ai0"] = {
            "snippet": {
                key: f"projects/sat-io/open-datasets/global_ai/{key}"
                for key in ["global_ai_yearly", "global_ai_monthly"]
            },
            "bands": [],
            "resolution": 1024,
            "start_date": None,
            "end_date": None,
        }

    if "gwa" in product:
        sensors["gwa"] = {
            "snippet": {
                key: f"projects/earthengine-legacy/assets/projects/sat-io/open-datasets/global_wind_atlas/{key}"
                for key in [
                    "air-density",
                    "capacity-factor",
                    "power-density",
                    "ruggedness-index",
                    "wind-speed",
                ]
            },
            "bands": [],
            "resolution": 1024,
            "start_date": None,
            "end_date": None,
        }

    if "WorldCover" in product:
        sensors["WorldCover"] = {
            "snippet": "ESA/WorldCover/v200",
            "bands": ["Map"],
            "resolution": 128,
            "start_date": "2021-01-01",
            "end_date": "2022-01-01",
        }

    if "OpenLandMap" in product:
        sensors["OpenLandMap"] = {
            "snippet": "OpenLandMap/PNV/PNV_BIOME-TYPE_BIOME00K_C/v01",
            "bands": ["biome_type"],
            "resolution": 1024,
            "start_date": "2001-01-01",
            "end_date": "2002-01-01",
        }

    if "malaria" in product:
        sensors["malaria"] = {
            "snippet": {
                "level1": "projects/sat-io/open-datasets/IUCN_HABITAT/iucn_habitatclassification_composite_lvl1_ver004",
                "level2": "projects/sat-io/open-datasets/IUCN_HABITAT/iucn_habitatclassification_composite_lvl2_ver004",
            },
            "bands": [],
            "resolution": 5120,
            "start_date": "2001-01-01",
            "end_date": "2002-01-01",
        }

    if "habitat" in product:
        sensors["habitat"] = {
            "snippet": "Oxford/MAP/IGBP_Fractional_Landcover_5km_Annual",
            "bands": [
                "Overall_Class",
                "Water",
                "Evergreen_Needleleaf_Forest",
                "Evergreen_Broadleaf_Forest",
                "Deciduous_Needleleaf_Forest",
                "Deciduous_Broadleaf_Forest",
                "Mixed_Forest",
                "Closed_Shrublands",
                "Open_Shrublands",
                "Woody_Savannas",
                "Savannas",
                "Grasslands",
                "Permanent_Wetlands",
                "Croplands",
                "Urban_And_Built_Up",
                "Cropland_Natural_Vegetation_Mosaic",
                "Snow_And_Ice",
                "Barren_Or_Sparsely_Populated",
                "Unclassified",
                "No_Data",
            ],
            "resolution": 1024,
            "start_date": None,
            "end_date": None,
        }

    if "GLO_DEM" in product:
        sensors["GLO_DEM"] = {
            "snippet": "COPERNICUS/DEM/GLO30",
            "bands": ["DEM", "EDM", "FLM", "HEM", "WBM"],
            "resolution": 128,
            "start_date": None,
            "end_date": None,
        }

    return Sensors(sensors=sensors)
