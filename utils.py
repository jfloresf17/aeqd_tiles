import pathlib
from calendar import monthrange
from datetime import datetime
from typing import Dict, List, Tuple

import ee
import geopandas as gpd
import requests
import tqdm

from datamodels import get_sensor_products

# Initialize Earth Engine with high volume
ee.Initialize(opt_url="https://earthengine-highvolume.googleapis.com")


def filter_by_month(
    start_date: str = "2019-01-01", end_date: str = "2022-12-31"
) -> List[List[Tuple[datetime, datetime]]]:
    """
    Filter the date range by month.
    Args:
    start_date (str): The start date.
    end_date (str): The end date.

    Returns:
    List[List[Tuple[datetime, datetime]]]: The list of list of tuples with the date range.
    """

    # Divide the range of dates by scale
    start_date = datetime.strptime(start_date, "%Y-%m-%d")
    end_date = datetime.strptime(end_date, "%Y-%m-%d")

    # List to store the dates
    filters = []

    current_date = datetime(start_date.year, start_date.month, start_date.day)
    while current_date <= end_date:
        last_day = monthrange(current_date.year, current_date.month)[1]
        interval_end = datetime(current_date.year, current_date.month, last_day)

        # To avoid end date greater than global end date
        if interval_end > end_date:
            interval_end = end_date

        filters.append((current_date, interval_end))

        # Get the next month
        current_date = datetime(
            current_date.year + current_date.month // 12,
            ((current_date.month % 12) + 1),
            1,
        )
    # Group by month
    months_filter = [
        [intervalo for intervalo in filters if intervalo[0].month == month]
        for month in range(1, 13)
    ]

    return months_filter


def apply_quality_filter(
    image: ee.Image, data_bands: List[str], quality_bands: List[str]
) -> ee.Image:
    """
    Apply quality filter to the image.
    Args:
    image (ee.Image): The image to apply the filter.
    data_bands (List[str]): The data bands to filter.
    quality_bands (List[str]): The quality bands to filter.

    Returns:
    ee.Image: The filtered image.
    """

    # Apply filter to each quality bands
    for dband, qband in zip(data_bands, quality_bands):
        # Get the quality band
        quality = image.select(qband)
        # Mask the data band
        data = image.select(dband).updateMask(quality.eq(0))
        # Replace the data band
        image = image.addBands(data, overwrite=True)

    return image


def save_composite(
    composite: ee.Image,
    aoi: ee.Geometry,
    crs: str,
    resolution: int,
    out_path: pathlib.Path,
    data_bands: List[str],
):
    """
    Save the composite to disk.
    Args:
    composite (ee.Image): The composite to save.
    aoi (ee.Geometry): The area of interest.
    crs (str): The coordinate reference system.
    resolution (int): The resolution of the composite.
    out_path (pathlib.Path): The output path to save the composite.
    data_bands (List[str]): The data bands to save.
    """

    # Get the composite as a download URL
    url = composite.getDownloadURL(
        {
            "scale": resolution,
            "region": aoi,
            "filePerBand": False,
            "format": "GeoTIFF",
            "crs": crs,
            "bands": data_bands,
        }
    )
    response = requests.get(url)

    # Save the composite
    if response.status_code == 200:
        with open(out_path, "wb") as file:
            file.write(response.content)
    else:
        raise ValueError(
            f"Error to download the composite: {response.status_code}"
        )


def download_modis(
    product: str,
    snippet: str,
    bands: str,
    start_date: str,
    end_date: str,
    aoi: ee.Geometry.Polygon,
    crs: str,
    resolution: int,
    folder: pathlib.Path,
):
    """
    Retrieves a MODIS image collection based on specified parameters and filters.

    Parameters:
    - product (str): The MODIS product type, e.g., "MODIS_BRDF" or "MODIS_LAND_SURFACE".
    - snippet (str): The snippet string for accessing the MODIS image collection.
    - bands (str): The bands to select from the MODIS images.
    - start_date (str): The start date for filtering the MODIS image collection (format: "YYYY-MM-DD").
    - end_date (str): The end date for filtering the MODIS image collection (format: "YYYY-MM-DD").
    - aoi (ee.Geometry.Polygon): The area of interest (AOI) as an Earth Engine geometry.
    - crs (str): The coordinate reference system (CRS) for the MODIS images.
    - resolution (int): The resolution for the MODIS images.
    - folder (pathlib.Path): The folder to save the composite images.
    """

    # Filtering by diferent MODIS products
    if product == "MODIS_BRDF":
        # Filter by quality bands
        quality_bands = bands[7:]
        data_bands = bands[:7]

    elif product == "MODIS_LAND_SURFACE":
        # Filter by quality bands
        quality_bands = bands[2:]
        data_bands = bands[:2]

    # Get end_date
    if end_date is None:
        end_date = "2023-12-31"

    # Get the collection
    collection = (
        ee.ImageCollection(snippet)
        .filterBounds(aoi)
        .filterDate(start_date, end_date)
        .select(bands)
    )

    # Apply quality filter
    quality_collection = collection.map(
        lambda image: apply_quality_filter(image, data_bands, quality_bands)
    )

    # Generate a compsity by months
    filter_months = filter_by_month(start_date, end_date)

    for i, months in enumerate(filter_months):
        monthly_composites = []
        for month in months:
            month_collection = quality_collection.filterDate(
                month[0].strftime("%Y-%m-%d"), month[1].strftime("%Y-%m-%d")
            )
            # Composite the collection
            monthly_composite = month_collection.median()
            monthly_composites.append(monthly_composite)

        # Combine the monthly composites to generate a new collection
        composite = ee.ImageCollection(monthly_composites).median()

        # Generate the out path for each month
        month_name = datetime(2000, i + 1, 1).strftime("%B")
        out_path = folder / f"{product}_{month_name}.tif"

        # Save the composite
        save_composite(composite, aoi, crs, resolution, out_path, data_bands)


def download_product(
    product: str,
    snippet: str,
    bands: List[str],
    start_date: str,
    end_date: str,
    aoi: ee.Geometry.Polygon,
    crs: str,
    resolution: int,
    folder: pathlib.Path,
):
    """
    Gets products from an Earth Engine image collection based on specified parameters.

    Parameters:
    - product (str): The name of the product.
    - snippet (str): The snippet for accessing the image collection.
    - bands (List[str]): The bands to select from the images.
    - start_date (str): The start date for filtering the image collection (format: "YYYY-MM-DD").
    - end_date (str): The end date for filtering the image collection (format: "YYYY-MM-DD").
    - aoi (ee.Geometry.Polygon): The area of interest (AOI) as an Earth Engine geometry.
    - crs (str): The coordinate reference system (CRS) for the images.
    - resolution (int): The resolution for the images.
    - folder (pathlib.Path): The folder to save the image.
    """

    # Check if there are no bands or any date is missing
    if not bands or start_date is None or end_date is None:
        # Filter the image by the region of interest and get the first image
        image = ee.ImageCollection(snippet).filterBounds(aoi).median()
    else:
        # If the product is "OpenLandMap", clip the image by the region of interest
        # Otherwise, filter the collection by the region of interest, date, and bands
        if product == "OpenLandMap":
            image = ee.Image(snippet).clip(aoi)
        else:
            image = (
                ee.ImageCollection(snippet)
                .filterBounds(aoi)
                .filterDate(start_date, end_date)
                .select(bands)
                .median()
            )

    out_path = folder / f"{product}.tif"
    save_composite(image, aoi, crs, resolution, out_path, bands)


def download_subproduct(
    product: str,
    snippet: Dict[str, str],
    bands: List[str],
    aoi: ee.Geometry.Polygon,
    crs: str,
    resolution: int,
    folder: pathlib.Path,
):
    """
    Gets sub-products from an Earth Engine image collection based on specified parameters.

    Parameters:
    - product (str): The name of the product.
    - snippet (Dict[str, str]): A dictionary containing the snippets for accessing the image collections.
    - bands (str): The bands to select from the images.
    - aoi (ee.Geometry.Polygon): The area of interest (AOI) as an Earth Engine geometry.
    - crs (str): The coordinate reference system (CRS) for the images.
    - resolution (int): The resolution for the images.
    - folder (pathlib.Path): The folder to save the image.
    """

    for subproduct, snippet2 in snippet.items():
        if product in ["malaria", "isric"] or subproduct in [
            "ruggedness-index",
            "global_ai_yearly",
        ]:
            # Clip the image by the region of interest
            composite = ee.Image(snippet2).clip(aoi)
        else:
            # Filter the image collection by the region of interest and compute the median
            composite = ee.ImageCollection(snippet2).filterBounds(aoi).median()

        subfolder = folder / product
        subfolder.mkdir(parents=True, exist_ok=True)

        out_path = subfolder / f"{subproduct}.tif"
        save_composite(composite, aoi, crs, resolution, out_path, bands)


def download(
    products: List[str], zone: str = "AF", T1_tile: str = "E0589N0557T1"
):
    """
    Downloads satellite imagery products for a specified T1 tile.

    Parameters:
    - products (List[str]): A list of product names to download.
    - zone (str, optional): The zone code. Default is "AF".
    - T1_tile (str, optional): The name of the T1 tile. Default is "E0589N0557T1".
    """

    # Open gpkg
    gdf_t1 = gpd.read_file(
        f"https://huggingface.co/datasets/jfloresf/landclip/resolve/main/grids/{zone}/{zone}.gpkg",
        layer="PROJ_T1",
    )
    shape = gdf_t1[gdf_t1["TILE"] == T1_tile]

    # Get tile codes
    search_path = f"./landclip/{zone}"
    matching_files = list(pathlib.Path(search_path).rglob("*T1*"))

    # Filtrar para obtener el archivo exacto que coincide con T1_tile
    code = next((file for file in matching_files if file.name == T1_tile), None)

    # Get PROJ crs
    crs = gdf_t1.crs.to_string()

    # Convert to lat/lon
    geog_shape = shape.to_crs("EPSG:4326")
    shape = geog_shape.geometry.values[0]

    # Get BBox geometry
    buffer_shape = list(shape.exterior.coords.xy)
    buffer_coords = [[x, y] for x, y in zip(buffer_shape[0], buffer_shape[1])]
    aoi = ee.Geometry.Polygon(buffer_coords)

    sensor_metadata = get_sensor_products()

    # Download the products
    for product in tqdm.tqdm(products):
        snippet = sensor_metadata.sensors[product].snippet
        resolution = sensor_metadata.sensors[product].resolution
        bands = sensor_metadata.sensors[product].bands
        start_date = sensor_metadata.sensors[product].start_date
        end_date = sensor_metadata.sensors[product].end_date

        # Append images collections
        if "MODIS" in product:
            download_modis(
                product,
                snippet,
                bands,
                start_date,
                end_date,
                aoi,
                crs,
                resolution,
                code,
            )

        elif product in ["s1gbm", "isric", "ai0", "gwa", "malaria"]:
            download_subproduct(product, snippet, aoi, crs, resolution, code)

        else:
            download_product(
                product,
                snippet,
                bands,
                start_date,
                end_date,
                aoi,
                crs,
                resolution,
                code,
            )

        print(f"{product} downloaded successfully!")
