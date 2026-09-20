"""
DEM-based watershed and river network delineation for the South Bali
study area, using WhiteboxTools' D8 hydrology toolchain.

Pipeline: breach depressions -> D8 flow direction -> automatic basin
delineation -> raster-to-vector conversion -> merge micro-catchments
(area < 10 km^2) into their largest neighbour.

Output: a cleaned catchment-boundary shapefile, one polygon per
sub-catchment, ready to feed into hydrograph_model.py (see the
accompanying README).

Designed to run in Google Colab against a Google Drive working
directory.
"""

# ==========================================
# 0. INSTALLATION & ENVIRONMENT SETUP
# ==========================================
# Run this pip line if the libraries are not yet installed in your Colab runtime
!pip install whitebox rasterio geopandas shapely scipy -q

import os
import whitebox
import rasterio
from rasterio.features import shapes
import geopandas as gpd
import warnings
from shapely.errors import ShapelyDeprecationWarning

# Suppress shapely deprecation warnings to keep the console output clean
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)

def setup_working_environment():
    """Initialise WhiteboxTools and the Google Drive working directory."""
    wbt = whitebox.WhiteboxTools()
    wbt.set_verbose_mode(False)

    # Point WhiteboxTools at the Google Drive working directory
    work_dir = "/content/drive/MyDrive"
    wbt.set_working_dir(work_dir)
    return wbt, work_dir

def run_hydrological_pipeline(wbt, work_dir, dem_input):
    """Raster hydrology stages: breach depressions -> D8 flow direction -> basins."""
    print("Stage 1: Breach depressions (conditioning the DEM for hydrology)...")
    dem_breached = "01_dem_breached.tif"
    wbt.breach_depressions(dem=dem_input, output=dem_breached)

    print("Stage 2: Computing D8 flow direction...")
    fdir_raster = "02_flow_direction.tif"
    wbt.d8_pointer(dem=dem_breached, output=fdir_raster)

    print("Stage 3: Automatically delineating all basins (watersheds)...")
    basins_raster = "03_basins_raw.tif"
    wbt.basins(d8_pntr=fdir_raster, output=basins_raster)

    return os.path.join(work_dir, basins_raster)

def raster_to_geodataframe(raster_path, nodata_val=None):
    """Convert the delineated basin raster into vector polygons."""
    print("Stage 4: Converting the basin raster into vector polygons...")
    with rasterio.open(raster_path) as src:
        image = src.read(1)
        transform = src.transform
        crs = src.crs

        # Mask out nodata and background (value <= 0) pixels
        if nodata_val is not None:
            mask = (image != nodata_val) & (image > 0)
        else:
            mask = image > 0

        results = (
            {'properties': {'gridcode': v}, 'geometry': s}
            for i, (s, v)
            in enumerate(shapes(image, mask=mask, transform=transform))
        )

        gdf = gpd.GeoDataFrame.from_features(list(results), crs=crs)
        return gdf

def eliminate_micro_catchments(gdf, min_area_km2=10.0, target_epsg=32750):
    """Merge catchment polygons smaller than min_area_km2 into their largest neighbour."""
    print(f"Stage 5: Cleaning topology and eliminating micro-catchments (< {min_area_km2} km2)...")

    # Reproject to UTM (EPSG:32750, central Indonesia / Bali zone) for accurate area in square metres
    if gdf.crs is None or gdf.crs.to_epsg() != target_epsg:
        gdf = gdf.to_crs(epsg=target_epsg)

    # Compute initial area in km2
    gdf['area_km2'] = gdf.geometry.area / 10**6
    initial_count = len(gdf)
    print(f"Initial number of catchment polygons: {initial_count}")

    iteration = 1
    while (gdf['area_km2'] < min_area_km2).any():
        idx_smallest = gdf[gdf['area_km2'] < min_area_km2]['area_km2'].idxmin()
        smallest_geom = gdf.loc[idx_smallest, 'geometry']

        # Find neighbours using a 1 m buffer to catch small gaps/slivers between polygons
        neighbors = gdf[gdf.geometry.intersects(smallest_geom.buffer(1))]
        neighbors_idx = [i for i in neighbors.index if i != idx_smallest]

        if neighbors_idx:
            # Identify the neighbour with the largest area
            largest_neighbor_idx = gdf.loc[neighbors_idx, 'area_km2'].idxmax()

            # Merge the micro-catchment geometry into its dominant neighbour
            merged_geom = smallest_geom.union(gdf.loc[largest_neighbor_idx, 'geometry'])
            gdf.loc[largest_neighbor_idx, 'geometry'] = merged_geom

            # Drop the now-merged micro-catchment
            gdf = gdf.drop(index=idx_smallest)
            gdf['area_km2'] = gdf.geometry.area / 10**6
        else:
            # If the polygon is genuinely isolated, skip it to avoid an infinite loop
            gdf.loc[idx_smallest, 'area_km2'] = min_area_km2 + 0.1

        if iteration % 100 == 0:
            print(f"Iteration {iteration}: {len(gdf)} catchments remaining...")
        iteration += 1

    print(f"Elimination complete. Final catchment count: {len(gdf)} (reduced by {initial_count - len(gdf)} polygons)")
    return gdf

# ==========================================
# MAIN EXECUTION
# ==========================================
if __name__ == "__main__":
    # 1. Initial configuration
    wbt, work_dir = setup_working_environment()

    # Raw DEM filename in your Google Drive
    dem_input_file = "DEM BALI BARU.tif"
    output_shapefile_name = "Batas_DAS_Bali_Cleaned_10km2.shp"

    # Analysis parameters
    MIN_AREA_THRESHOLD_KM2 = 10.0
    UTM_ZONE_EPSG = 32750  # UTM Zone 50S for the Bali / central Indonesia region

    print(f"Starting automated delineation for: {dem_input_file}")

    # 2. Run the raster hydrology pipeline (WhiteboxTools)
    basins_raster_path = run_hydrological_pipeline(wbt, work_dir, dem_input_file)

    # 3. Convert the basin raster to vector polygons
    gdf_raw_vectors = raster_to_geodataframe(basins_raster_path)

    # 4. Eliminate micro-catchments below the area threshold
    gdf_final_cleaned = eliminate_micro_catchments(
        gdf_raw_vectors,
        min_area_km2=MIN_AREA_THRESHOLD_KM2,
        target_epsg=UTM_ZONE_EPSG
    )

    # 5. Export the final result to Google Drive
    final_output_path = os.path.join(work_dir, output_shapefile_name)
    gdf_final_cleaned.to_file(final_output_path)

    print(f"\nPROCESSING COMPLETE!")
    print(f"Final shapefile saved to: {final_output_path}")
