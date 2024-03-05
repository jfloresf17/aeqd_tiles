from utils import download

## ALL PRODUCTS
# products = ["MODIS_BRDF", #date and bands
#             "WorldCLIM", #date and bands
#             "MODIS_LAND_SURFACE", #date and bands
#             "CCNL", #date and bands
#             "fabdem", #no date and no bands
#             "geomorpho90m", #no date and no bands
#             "s1gbm", #no date and no bands
#             "isric", #no date and no bands
#             "soil_salinity", #no date and no bands
#             "ai0", #no date and no bands
#             "gwa", #no date and no bands
#             "WorldCover", #date and bands
#             "OpenLandMap", #date and bands
#             "malaria", #no bands
#             "habitat", #no date
#             "GLO_DEM" #no date
#           ]

sample_products = [
    "MODIS_BRDF",  # date and bands
    "MODIS_LAND_SURFACE",  # date and bands
    "GLO_DEM",  # no date
]

download(sample_products, zone="SA", T1_tile="E0589N0557T1")
