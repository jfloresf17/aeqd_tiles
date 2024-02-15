import pathlib
from calendar import monthrange
from datetime import datetime
from typing import List, Optional, Tuple

import ee
import geopandas as gpd
import numpy as np
import pyproj
import requests
import shapely
import tqdm
import xarray as xr

from datamodels import get_equi7_tiles, get_sensor_products

# Initialize Earth Engine with high volume
ee.Initialize(opt_url="https://earthengine-highvolume.googleapis.com")


def generate_grids(
    tile_name: str = "AF",
    tile_extent: float = 500 * 128,
    dst_path: pathlib.Path = pathlib.Path("./output"),
) -> List[List[float]]:
    """
    Generate a grid of points for a given Equi7Grid tile name.

    Args:
    tile_name (str): The Equi7Grid tile name.
    tile_extent (float): The extent of the grid cell in meters.
    dst_path (pathlib.Path): The destination path to save the grid cells.

    Returns:
    List[List[float]]: A list of lists with the coordinates of the
    grid cells.
    """

    # Get the tile grid
    equi7 = get_equi7_tiles()
    tiles = equi7.tiles[tile_name]
    gdf = gpd.read_file(tiles.grid_path)

    # Get the PROJ4-string
    with open(tiles.prj_path, "r") as file:
        crs = file.read()

    # Make a grid from the layer
    xmin, ymin, xmax, ymax = gdf.total_bounds

    # create the cells in a loop
    grid_cells = []
    for x0 in np.arange(xmin, xmax + tile_extent, tile_extent):
        for y0 in np.arange(ymin, ymax + tile_extent, tile_extent):
            # bounds
            x1 = x0 - tile_extent
            y1 = y0 + tile_extent

            # create the shapely geometry
            cell_box = shapely.geometry.box(x0, y0, x1, y1)

            # append the centroid to the list
            grid_cells.append(cell_box)

    # Create a GeoDataFrame
    cell = gpd.GeoDataFrame(geometry=grid_cells)

    # Intersect with the land tile layer
    gdf_land = gdf[gdf["COVERSLAND"] == 1]
    cell_land = cell[cell.intersects(gdf_land.unary_union)]

    # Intersect with the countries layer
    countries = gpd.read_file(tiles.land_path)
    cell_countries = cell_land[cell_land.intersects(countries.unary_union)]

    # Get centroids
    cell_points = cell_countries.copy()
    cell_points["geometry"] = cell_points["geometry"].centroid

    # Set projection
    cell_points.set_crs(crs, inplace=True)

    # Generate a buffer around each point
    buffer = cell_points.buffer(tile_extent / 2).envelope
    buffer_list = []
    for geom in tqdm.tqdm(buffer):
        minx, miny, maxx, maxy = geom.bounds
        polygon = [[minx, miny], [minx, maxy], [maxx, maxy], [maxx, miny]]

        # Convert to latlon
        transformer = pyproj.Transformer.from_crs(
            crs, "EPSG:4326", always_xy=True
        )
        latlon_coords = [
            list(transformer.transform(x[0], x[1])) for x in polygon
        ]

        # append to list
        buffer_list.append(latlon_coords)

    # Save buffer list as GeoJSON
    gdf_buffer = gpd.GeoDataFrame(
        data={"id": range(len(buffer_list))},
        geometry=gpd.GeoSeries(
            [shapely.Polygon(coords) for coords in buffer_list]
        ),
        crs="EPSG:4326",
    )

    folder_path = dst_path / tile_name
    folder_path.mkdir(parents=True, exist_ok=True)

    gdf_buffer.to_file(
        f"{folder_path}/{tile_name}_buffers.geojson", driver="GeoJSON"
    )


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


def save_composite(
    composite: ee.Image,
    aoi: ee.Geometry,
    crs: str,
    resolution: int,
    out_path: str,
    data_bands: List[str],
):
    """
    Save the composite to disk.
    Args:
    composite (ee.Image): The composite to save.
    aoi (ee.Geometry): The area of interest.
    crs (str): The coordinate reference system.
    resolution (int): The resolution of the composite.
    out_path (str): The output path to save the composite.
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


def composite_collection(
    geom: shapely.geometry.Polygon,
    tile_name: str = "AF",
    tile_id: int = 0,
    product: str = "MODIS_061_MCD43A4",
    subproduct: Optional[str] = None,
    start_date: str = None,
    end_date: str = None,
    dst_folder: Optional[pathlib.Path] = pathlib.Path("./output"),
):
    """
    Composite the collection.

    Args:
    geom (shapely.geometry.Polygon): The geometry to filter the collection.
    tile_name (str): The Equi7Grid tile name.
    tile_id (int): The tile id for grid.
    product (str): The product name.
    subproducts (str): The subproduct name of selected product. Default is None.
    time (str): The time unit to composite by unit time. Can be "year", "month", "week", "day".
    start_date (str): The start date.
    end_date (str): The end date.
    dst_folder (pathlib.Path): The destination folder to save the composite. Default is "./output".
    """

    # Convert shapely.geometry.Polygon to ee.Geometry.Polygon
    shape = list(geom.exterior.coords.xy)
    coords = [[x, y] for x, y in zip(shape[0], shape[1])]
    polygon = ee.Geometry.Polygon(coords)

    # Get AOI
    aoi = polygon.getInfo()

    # Get the sensor metadata
    sensor_metadata = get_sensor_products()
    snippet = sensor_metadata.sensors[product].snippet
    resolution = sensor_metadata.sensors[product].resolution

    # Get the bands
    bands = sensor_metadata.sensors[product].bands

    # Get start and end date
    if start_date is None and end_date is None:
        start_date = sensor_metadata.sensors[product].start_date
        end_date = sensor_metadata.sensors[product].end_date

    # Get the CRS
    equi7 = get_equi7_tiles()
    tiles = equi7.tiles[tile_name]
    with open(tiles.prj_path, "r") as file:
        crs = file.read()

    # Get the folder to save the composite
    tile_id = str(tile_id).zfill(4)
    folder = dst_folder / f"{tile_name}/{product}"
    folder = pathlib.Path(folder)
    folder.mkdir(parents=True, exist_ok=True)

    # Download the composite
    if "MODIS" in product:
        # Filtering by diferent MODIS products
        if "MODIS_061_MCD43A4" in product:
            # Filter by quality bands
            quality_bands = bands[7:]
            data_bands = bands[:7]

        elif "MODIS_061_MOD11A1" in product:
            # Filter by quality bands
            quality_bands = bands[2:]
            data_bands = bands[:2]

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

        for i, months in tqdm.tqdm(enumerate(filter_months)):
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

            # Generate the out path
            month_name = datetime(2000, i + 1, 1).strftime("%B")
            out_path = folder / f"{tile_id}_{month_name}.tif"

            # Save the composite
            save_composite(
                composite, aoi, crs, resolution, out_path, data_bands
            )

    elif product in ["s1gbm", "isric", "ai0", "gwa"]:
        snippet2 = snippet[subproduct]
        collection = ee.ImageCollection(snippet2).filterBounds(aoi)
        composite = collection.median()

        # Generate the out path
        out_path = folder / subproduct
        pathlib.Path(out_path).mkdir(parents=True, exist_ok=True)

        out_path = out_path / f"{tile_id}.tif"

        # Save the composite
        save_composite(composite, aoi, crs, resolution, out_path, bands)

    else:
        collection = ee.ImageCollection(snippet).filterBounds(aoi)
        # Composite the collection
        composite = collection.median()

        out_path = folder / f"{tile_id}.tif"

        # Save the composite
        save_composite(composite, aoi, crs, resolution, out_path, bands)
