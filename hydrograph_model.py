"""
Master hydrological modelling script for 14 sub-catchments in South Bali
(Denpasar, Badung, Gianyar), reconstructing the design flood hydrograph
for the 9-10 September 2025 extreme rainfall event.

Pipeline: continuous SCS-CN infiltration (37-hour storm, Mononobe +
Alternating Block Method disaggregation) -> Snyder unit hydrograph
(Kirpich time of concentration, true centroidal distance) -> hourly
convolution -> per-catchment design flood hydrograph, exported to a
multi-sheet Excel workbook.

Designed to run in Google Colab against a Google Drive folder structure;
see the accompanying README for expected input file names and formats.

Companion script: catchment_and_river_delineation.py (DEM -> watershed
and river network delineation, upstream of this script).
"""

import os
import glob
import warnings
import io
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import networkx as nx
from google.colab import drive
from scipy.interpolate import griddata
from scipy.spatial import cKDTree
from shapely.geometry import LineString, MultiLineString, Point
from shapely.ops import snap, unary_union, linemerge
import openpyxl
from openpyxl.drawing.image import Image as OpenpyxlImage
from matplotlib.ticker import FuncFormatter

try:
  import rasterio
  import rasterio.mask
except ImportError:
  !pip install rasterio
  import rasterio
  import rasterio.mask

try:
  import contextily as cx
except ImportError:
  !pip install contextily
  import contextily as cx

warnings.filterwarnings('ignore')

# =====================================================================
# 1. MOUNT DRIVE & PATH SETUP
# =====================================================================
drive.mount('/content/drive')

DRIVE_ROOT = '/content/drive/MyDrive/'
RAINFALL_CSV = os.path.join(DRIVE_ROOT, 'Curah Hujan Bali Banjir September 2025.csv')
LANDCOVER_CSV = os.path.join(DRIVE_ROOT, 'Tata Guna Lahan & Manning 2025.csv')  # UPDATE FILENAME HERE IF CHANGED
CATCHMENT_SHP_DIR = os.path.join(DRIVE_ROOT, 'DAS Phyton - Tahun 1')
RIVER_SHP_DIR = os.path.join(DRIVE_ROOT, 'Sungai Phyton - Tahun 1')
DEM_DIR = os.path.join(DRIVE_ROOT, 'Bali Elevation - Tahun 1')
SOIL_RASTER_DIR = os.path.join(DRIVE_ROOT, 'Jenis Tanah')  # HYSOGs250m hydrologic soil group raster
EXCEL_OUTPUT_PATH = os.path.join(DRIVE_ROOT, 'Hasil_Analisis_Hidrologi_14_DAS_Master_Final_Manning2025.xlsx')

dem_candidates = glob.glob(os.path.join(DEM_DIR, '**/*.tif'), recursive=True) + \
                 glob.glob(os.path.join(DEM_DIR, '**/*.dem'), recursive=True)
dem_raster_path = dem_candidates[0] if dem_candidates else None

hysogs_candidates = glob.glob(os.path.join(SOIL_RASTER_DIR, '**/*.tif'), recursive=True)
hysogs_raster_path = hysogs_candidates[0] if hysogs_candidates else None
if hysogs_raster_path is None:
  print("WARNING: HYSOGs250m raster not found in the 'Jenis Tanah' folder. "
        "Hydrologic Soil Group will fall back to a Group B assumption for every catchment.")

rainfall_df = pd.read_csv(RAINFALL_CSV)

# --- Parse the date column into real datetime objects (not string matching).
#     pd.to_datetime auto-detects common formats ('9/9/2025', '2025-09-09',
#     '09/09/2025', etc.) without needing to know the format in advance --
#     far more robust than str.startswith() against a fixed pattern. ---
rainfall_df['date'] = pd.to_datetime(rainfall_df['tgl1'], errors='coerce')

n_before = len(rainfall_df)
rainfall_df = rainfall_df.dropna(subset=['lon', 'lat', 'ch', 'date'])
if len(rainfall_df) < n_before:
  print(f"WARNING: {n_before - len(rainfall_df)} rows dropped from "
        f"{RAINFALL_CSV} due to missing lon/lat/rainfall/date values.")

EVENT_DAY1 = pd.Timestamp('2025-09-09')
EVENT_DAY2 = pd.Timestamp('2025-09-10')

rainfall_day1 = rainfall_df[rainfall_df['date'].dt.normalize() == EVENT_DAY1].copy()
rainfall_day2 = rainfall_df[rainfall_df['date'].dt.normalize() == EVENT_DAY2].copy()

if rainfall_day1.empty or rainfall_day2.empty:
  available_dates = sorted(rainfall_df['date'].dt.strftime('%Y-%m-%d').unique())[:15]
  raise ValueError(
    "FAILED TO LOAD RAINFALL DATA: no rows found dated "
    f"{EVENT_DAY1.date()} or {EVENT_DAY2.date()} in "
    f"{RAINFALL_CSV}.\n"
    f"  Total rows in file : {len(rainfall_df)}\n"
    f"  Rows matching day 1: {len(rainfall_day1)}\n"
    f"  Rows matching day 2: {len(rainfall_day2)}\n"
    f"  Dates actually present in file: {available_dates}\n"
    "  -> Compare the dates above with EVENT_DAY1/EVENT_DAY2 in this "
    "script; update them if your flood event is not 9-10 September.")

rain_gauge_coords = rainfall_day1[['lon', 'lat']].values
rain_gauge_tree = cKDTree(rain_gauge_coords)

# =====================================================================
# 1b. PARAMETER CONFIGURATION (values most likely to need tuning)
# =====================================================================
# --- Curve Number by land-cover class x Hydrologic Soil Group (NRCS
#     TR-55 / NEH-630 Ch.9, USDA 1984). Column B of this table matches
#     the old fixed CN_WATER/BUILTUP/OPEN/FOREST constants (100/85/61/55)
#     exactly, so results are unchanged for any catchment whose dominant
#     HSG happens to be B. Source rows: Forest="Woods, good condition";
#     Open="Pasture, good condition"; Built-up="Residential 1/8 acre lot,
#     65% impervious"; Water=100 for every HSG (no infiltration). ---
CN_TABLE_BY_HSG = {
  'water':   {'A': 100, 'B': 100, 'C': 100, 'D': 100},
  'forest':  {'A': 30,  'B': 55,  'C': 70,  'D': 77},
  'open':    {'A': 39,  'B': 61,  'C': 74,  'D': 80},
  'builtup': {'A': 77,  'B': 85,  'C': 90,  'D': 92},
}

# --- Reclassification of HYSOGs250m raster pixel values -> HSG class.
#     Combined/dual (drained-dependent) classes 11=A/D, 12=B/D, 13=C/D,
#     14=D/D are all mapped to D (conservative, undrained scenario). ---
HSG_PIXEL_MAP = {1: 'A', 2: 'B', 3: 'C', 4: 'D', 11: 'D', 12: 'D', 13: 'D', 14: 'D'}

# --- Snyder unit hydrograph basin coefficient: TL = Ct * (L * Lc)^0.3
#     Ct: literature range 1.8-2.2 for ungauged basins (Snyder, 1938);
#     formula used without the SI unit-conversion factor C1 (see
#     methodology notes in the accompanying manuscript). ---
SNYDER_CT = 2.0

# --- Kirpich (1940) time of concentration -> used to derive Tr (unit
#     storm duration), NOT the Snyder TL itself. Valid for small,
#     well-defined-channel catchments (Roussel et al., 2005: 0.25-150
#     sq mi, S=0.002-0.1); the 14 Bali catchments (1.1-102 km^2) fall
#     within this range, though accuracy is lower for the flattest
#     coastal catchments.
#     Tc [hours] = 0.06628 * L[km]^0.77 / S[m/m]^0.385  (SI-coefficient
#     version; L is used directly in km, no conversion to miles needed)
#     Tr = TR_FRACTION_OF_TC * Tc (literature range 0.1-0.3; 0.1, the
#     lower bound, is used here because it closely matches the Snyder
#     convention tr = TL/5.5 with TL ~ 0.6*Tc, giving ~0.109*Tc). ---
TR_FRACTION_OF_TC = 0.1

# --- Seasonal column used for AMC classification (Table 6.6, standard
#     SCS/NEH-4 antecedent-moisture table). DECIDED: 'dry', based on
#     BMKG (2025), "Prakiraan Musim Kemarau 2025 di Indonesia"
#     (Directorate of Climate Change, Deputy for Climatology), zone
#     BALI_19 -- whose boundary (southern Badung, southern Gianyar,
#     southern Tabanan, city of Denpasar) matches this study area
#     exactly. That zone's 2025 dry season runs from the third dasarian
#     of April (one dasarian later than normal) for 19 dasarians (one
#     shorter than normal), i.e. roughly 21 April - 31 October 2025.
#     The 9-10 September event falls in dasarian 14-15 of 19 -- still
#     within the official dry season, not the wet season. ---
AMC_SEASON = 'dry'  # 'wet' or 'dry'

# Table 6.6 thresholds (5-day antecedent rainfall total, cm -> mm)
AMC_THRESHOLD_MM = {
  'dry': {'I': 13.0, 'III': 28.0},   # <1.3 cm ; >2.8 cm
  'wet': {'I': 36.0, 'III': 53.0},   # <3.6 cm ; >5.3 cm
}

def classify_amc(antecedent_r5_mm, season=AMC_SEASON):
  """Classify Antecedent Moisture Condition from the 5-day antecedent
  rainfall total [mm], using the Table 6.6 thresholds for the given
  season column (see AMC_SEASON citation note above)."""
  if antecedent_r5_mm is None:
    return 'II'  # fallback: no antecedent data available -> assume normal condition
  threshold = AMC_THRESHOLD_MM[season]
  if antecedent_r5_mm < threshold['I']:
    return 'I'
  if antecedent_r5_mm <= threshold['III']:
    return 'II'
  return 'III'

def adjust_cn_for_amc(cn_ii, amc_class):
  """Convert CN(II) to the CN for the given AMC class, standard SCS/NEH-4
  formulas (Chow, Maidment & Mays, 1988, Applied Hydrology, Ch. 5)."""
  if amc_class == 'I':
    return (4.2 * cn_ii) / (10 - 0.058 * cn_ii)
  elif amc_class == 'III':
    return (23 * cn_ii) / (10 + 0.13 * cn_ii)
  return cn_ii  # AMC II -> unchanged

def compute_antecedent_r5(catchment_union_wgs84, full_rainfall_df, event_start=EVENT_DAY1):
  """Total rainfall over the 5 days preceding the event [mm], using the
  'date' column parsed into real datetime objects in Section 1.

  IMPORTANT (consistency fix): this function used to check whether the
  ORIGINAL, sparse rain-gauge grid points (~11 km apart) fell exactly
  inside the catchment polygon; for small catchments or unlucky
  positioning, no original point would fall inside, silently returning
  'n/a' even though the data existed in the CSV. Day-1/Day-2 rainfall
  never had this problem because both are interpolated onto a dense
  400x400 mesh before being checked against each catchment. This
  function now uses the SAME approach: interpolate onto the global
  mesh (grid_mesh_gdf/mesh_x/mesh_y, built in Section 1), then check
  which mesh points fall inside the catchment."""
  window_start = event_start - pd.Timedelta(days=5)
  window_end = event_start - pd.Timedelta(days=1)

  window_mask = (full_rainfall_df['date'].dt.normalize() >= window_start) & \
                (full_rainfall_df['date'].dt.normalize() <= window_end)
  if window_mask.sum() == 0:
    return None  # the 5 days preceding the event are not yet in the CSV

  # Sum rainfall over the 5 days PER original gauge point first
  antecedent_df = full_rainfall_df[window_mask]
  point_totals = antecedent_df.groupby(['lon', 'lat'])['ch'].sum().reset_index()
  if len(point_totals) < 3:
    return None  # too few points for a meaningful linear interpolation

  antecedent_coords = point_totals[['lon', 'lat']].values
  antecedent_values = point_totals['ch'].values

  # Interpolate onto the SAME global mesh used for day-1/day-2 rainfall
  # (mesh_x, mesh_y, grid_mesh_gdf from Section 1) -- keeps the method consistent
  interp_linear = griddata(antecedent_coords, antecedent_values, (mesh_x, mesh_y), method='linear')
  interp_nearest = griddata(antecedent_coords, antecedent_values, (mesh_x, mesh_y), method='nearest')
  mesh_r5 = grid_mesh_gdf.copy()
  mesh_r5['ch'] = np.where(np.isnan(interp_linear), interp_nearest, interp_linear).ravel()

  within_catchment = mesh_r5[mesh_r5.geometry.within(catchment_union_wgs84)]
  if within_catchment.empty:
    # Last-resort fallback: nearest mesh point to the catchment centroid
    # (should almost never trigger given how dense the 400x400 mesh is)
    catchment_centroid = catchment_union_wgs84.centroid
    mesh_coords = np.vstack([mesh_r5.geometry.x, mesh_r5.geometry.y]).T
    mesh_tree = cKDTree(mesh_coords)
    _, nearest_idx = mesh_tree.query([catchment_centroid.x, catchment_centroid.y])
    return float(mesh_r5['ch'].iloc[nearest_idx])

  return float(within_catchment['ch'].mean())


def compute_kirpich_tc(main_channel_km, channel_slope_frac):
  """Kirpich (1940) time of concentration, SI-coefficient version.
  Tc [HOURS] = 0.06628 * L[km]^0.77 / S[m/m]^0.385
  IMPORTANT (bug fix note): with coefficient 0.06628 and L in km, this
  formula gives Tc directly in HOURS, not minutes. An earlier version
  of this script incorrectly divided the result by 60, making Tc 60x
  too small (e.g. a 43 km channel came out as Tc ~ 5 minutes instead of
  ~5 hours). Sanity check: for L=1 km, S=1% (0.01), this formula gives
  Tc=0.39 hours (~23 minutes), a sensible value for a small, steep
  catchment; interpreted as minutes (0.39 min = 23 seconds) it would be
  physically absurd."""
  if main_channel_km <= 0 or channel_slope_frac <= 0:
    return None
  tc_hours = 0.06628 * (main_channel_km ** 0.77) / (channel_slope_frac ** 0.385)
  return tc_hours

# --- Initial abstraction ratio (Ia/S = lambda). Historical SCS 1956
#     default: 0.2. Hawkins et al. (2002) recommend 0.05 based on 307
#     gauged catchments / 28,301 storm events. Set to 0.2 to revert to
#     the classic ratio. ---
LAMBDA_IA = 0.05

def convert_cn_to_lambda005(cn_lambda020):
  """Convert a Curve Number calibrated on the conventional lambda=0.2
  basis (the standard NRCS tables, including CN_TABLE_BY_HSG above)
  into its equivalent under a lambda=0.05 basis.
  Equation 9, Hawkins, R.H., Jiang, R., Woodward, D.E., Hjelmfelt, A.T.,
  Van Mullem, J.A. (2002). "Runoff Curve Number Method: Examination of
  the Initial Abstraction Ratio." Full text (verified directly):
  https://ponce.sdsu.edu/hawkins_initial_abstraction.pdf
  CN0.05 = 100 / (1.879*(100/CN0.20 - 1)^1.15 + 1)
  Valid up to CN0.20 ~ 98.5; above that CN0.05 ~ CN0.20 (S -> 0)."""
  if cn_lambda020 >= 100:
    return 100.0
  return 100.0 / (1.879 * ((100.0 / cn_lambda020 - 1) ** 1.15) + 1)

# =====================================================================
# 2. LAND COVER PROCESSING (area percentage per class)
#    CN(II) and S(II) are NOT computed here any more -- they now depend
#    on the dominant HSG per catchment (from the HYSOGs250m raster), so
#    that calculation has moved inside the per-catchment loop (Section 5).
# =====================================================================
CATCHMENT_GRIDCODES = [5931, 6310, 6553, 6653, 6799, 7028, 7278, 7498, 7514, 8489, 10247, 11445, 11905, 13106]

landcover_df = pd.read_csv(LANDCOVER_CSV)
if 'Gridcode' not in landcover_df.columns:
  landcover_df['Gridcode'] = CATCHMENT_GRIDCODES

# Ensure area columns are float, if present
for col in ['Air_m2', 'Bangunan_m2', 'Hutan_m2', 'Terbuka_m2']:
  if col in landcover_df.columns:
    landcover_df[col] = landcover_df[col].astype(np.float64)

landcover_df['Total_m2'] = landcover_df[['Air_m2', 'Bangunan_m2', 'Hutan_m2', 'Terbuka_m2']].sum(axis=1)

landcover_df['water_pct'] = (landcover_df['Air_m2'] / landcover_df['Total_m2']) * 100
landcover_df['builtup_pct'] = (landcover_df['Bangunan_m2'] / landcover_df['Total_m2']) * 100
landcover_df['forest_pct'] = (landcover_df['Hutan_m2'] / landcover_df['Total_m2']) * 100
landcover_df['open_pct'] = (landcover_df['Terbuka_m2'] / landcover_df['Total_m2']) * 100
landcover_pct_by_catchment = landcover_df.set_index('Gridcode')[['water_pct', 'builtup_pct', 'forest_pct', 'open_pct']].to_dict('index')

# =====================================================================
# 3. GEOSPATIAL DATA: SHAPEFILES & RAINFALL GRID
# =====================================================================
catchment_shp_list = glob.glob(os.path.join(CATCHMENT_SHP_DIR, '**/*.shp'), recursive=True)
river_shp_list = glob.glob(os.path.join(RIVER_SHP_DIR, '**/*.shp'), recursive=True)

catchment_gdfs = [gpd.read_file(p) for p in catchment_shp_list]
all_catchments_gdf = gpd.GeoDataFrame(pd.concat(catchment_gdfs, ignore_index=True), crs=catchment_gdfs[0].crs)

river_gdfs = [gpd.read_file(p).to_crs(epsg=32750) for p in river_shp_list]
all_rivers_gdf = gpd.GeoDataFrame(pd.concat(river_gdfs, ignore_index=True), crs='EPSG:32750')

all_catchments_wgs84 = all_catchments_gdf.to_crs(epsg=4326)
xmin, ymin, xmax, ymax = all_catchments_wgs84.total_bounds
mesh_x_1d = np.linspace(xmin, xmax, 400)
mesh_y_1d = np.linspace(ymin, ymax, 400)
mesh_x, mesh_y = np.meshgrid(mesh_x_1d, mesh_y_1d)
grid_coords_df = pd.DataFrame(np.vstack([mesh_x.ravel(), mesh_y.ravel()]).T, columns=['lon', 'lat'])
grid_mesh_gdf = gpd.GeoDataFrame(grid_coords_df, geometry=[Point(xy) for xy in zip(grid_coords_df['lon'], grid_coords_df['lat'])], crs='EPSG:4326')

day1_values = rainfall_day1['ch'].values
day2_values = rainfall_day2['ch'].values

day1_interp_linear = griddata(rain_gauge_coords, day1_values, (mesh_x, mesh_y), method='linear')
day1_interp_nearest = griddata(rain_gauge_coords, day1_values, (mesh_x, mesh_y), method='nearest')
day1_grid_result = grid_mesh_gdf.copy()
day1_grid_result['ch'] = np.where(np.isnan(day1_interp_linear), day1_interp_nearest, day1_interp_linear).ravel()

day2_interp_linear = griddata(rain_gauge_coords, day2_values, (mesh_x, mesh_y), method='linear')
day2_interp_nearest = griddata(rain_gauge_coords, day2_values, (mesh_x, mesh_y), method='nearest')
day2_grid_result = grid_mesh_gdf.copy()
day2_grid_result['ch'] = np.where(np.isnan(day2_interp_linear), day2_interp_nearest, day2_interp_linear).ravel()

# =====================================================================
# 4. DEM ELEVATION, MAIN-CHANNEL EXTRACTION & MAPPING HELPERS
# =====================================================================
def get_dem_elevation(point_utm, dem_path, label=''):
  """Read elevation from the DEM as the mean of a 3x3 pixel window
  around the point (rather than a single pixel), to dampen DEM noise
  (e.g. spurious negative-elevation artifacts near the coast). Falls
  back to the single centre pixel if the whole 3x3 window is nodata.
  IMPORTANT: any failure (corrupt file, CRS mismatch, point outside the
  raster extent) is now printed as an explicit warning -- an earlier
  version silently swallowed all exceptions and returned 0.0, which
  could make channel_slope wrong with no indication in the output."""
  if dem_path is None or not os.path.exists(dem_path):
    print(f"  [WARNING] DEM file not found for point {label} -> elevation set to 0.0")
    return 0.0
  try:
    with rasterio.open(dem_path) as src:
      pt_gdf = gpd.GeoDataFrame(geometry=[point_utm], crs='EPSG:32750').to_crs(src.crs)
      pt_geom = pt_gdf.geometry.iloc[0]
      row, col = src.index(pt_geom.x, pt_geom.y)
      data = src.read(1)
      nodata_val = src.nodata

      r0, r1 = max(0, row - 1), min(data.shape[0], row + 2)
      c0, c1 = max(0, col - 1), min(data.shape[1], col + 2)
      window = data[r0:r1, c0:c1].astype(np.float64)

      valid_mask = (window != nodata_val) if nodata_val is not None else np.ones_like(window, dtype=bool)
      valid_mask &= (window > -9999)
      window_valid = window[valid_mask]

      if window_valid.size > 0:
        return float(np.mean(window_valid))

      val = data[row, col]
      if val > -9999:
        return float(val)

      print(f"  [WARNING] All 3x3 window pixels at point {label} are nodata -> elevation set to 0.0")
      return 0.0
  except Exception as e:
    print(f"  [WARNING] Failed to read DEM at point {label}: {e} -> elevation set to 0.0")
    return 0.0

def extract_main_channel_standard(clipped_river_gdf, snap_tolerance_m=0.0):
  """Main-channel extraction for 13 of the 14 catchments.
  NOTE: catchment 7028 uses extract_main_channel_das7028() below instead
  of a variant of this function. The coordinate-rounding precision (0
  decimals), geometry-merging method (snap, not linemerge), number of
  candidate nodes (5), and lack of a longest-line fallback are all
  DELIBERATELY different from the 7028-specific version below -- do not
  merge/simplify these two without re-validating the resulting L and
  main_channel_geom for all 14 catchments."""
  if clipped_river_gdf.empty:
    return 0.0, None
  merged_union = unary_union(clipped_river_gdf.geometry)
  snapped_geom = merged_union if snap_tolerance_m == 0 else snap(merged_union, merged_union, tolerance=snap_tolerance_m)

  lines = [snapped_geom] if isinstance(snapped_geom, LineString) else (list(snapped_geom.geoms) if isinstance(snapped_geom, MultiLineString) else [])
  if not lines:
    return 0.0, None

  graph = nx.Graph()
  for line in lines:
    coords = list(line.coords)
    for i in range(len(coords) - 1):
      p1 = (round(coords[i][0], 0), round(coords[i][1], 0))
      p2 = (round(coords[i + 1][0], 0), round(coords[i + 1][1], 0))
      dist = Point(p1).distance(Point(p2))
      if dist > 0:
        graph.add_edge(p1, p2, weight=dist)

  if graph.number_of_nodes() == 0:
    return 0.0, None

  best_length, best_path = 0.0, None
  for component in nx.connected_components(graph):
    subgraph = graph.subgraph(component)
    if len(subgraph.nodes()) < 2:
      continue
    nodes_sorted = sorted(list(subgraph.nodes()), key=lambda n: n[1], reverse=True)
    top_nodes = nodes_sorted[:min(5, len(nodes_sorted))]
    bottom_nodes = nodes_sorted[-min(5, len(nodes_sorted)):]
    for upstream_node in top_nodes:
      for downstream_node in bottom_nodes:
        if upstream_node == downstream_node or upstream_node[1] <= downstream_node[1]:
          continue
        try:
          path = nx.shortest_path(subgraph, source=upstream_node, target=downstream_node, weight='weight')
          length = sum(subgraph[u][v]['weight'] for u, v in zip(path[:-1], path[1:]))
          if length > best_length:
            best_length, best_path = length, path
        except nx.NetworkXNoPath:
          continue

  return (best_length / 1000.0, LineString(best_path)) if best_path and len(best_path) > 1 else (0.0, None)

def extract_main_channel_das7028(clipped_river_gdf):
  """Main-channel extraction specific to catchment 7028 (see docstring
  of extract_main_channel_standard for why this is kept separate)."""
  if clipped_river_gdf.empty:
    return 0.0, None
  lines_only = []
  for geom in clipped_river_gdf.geometry:
    if isinstance(geom, LineString) and len(geom.coords) > 1:
      lines_only.append(geom)
    elif hasattr(geom, 'geoms'):
      for sub_geom in geom.geoms:
        if isinstance(sub_geom, LineString) and len(sub_geom.coords) > 1:
          lines_only.append(sub_geom)

  if not lines_only:
    return 0.0, None
  merged = linemerge(unary_union(lines_only))
  lines = [merged] if isinstance(merged, LineString) else (list(merged.geoms) if isinstance(merged, MultiLineString) else lines_only)

  graph = nx.Graph()
  for line in lines:
    coords = list(line.coords)
    for i in range(len(coords) - 1):
      p1 = (round(coords[i][0], 1), round(coords[i][1], 1))
      p2 = (round(coords[i + 1][0], 1), round(coords[i + 1][1], 1))
      dist = Point(p1).distance(Point(p2))
      if dist > 0:
        graph.add_edge(p1, p2, weight=dist)

  best_length, best_path = 0.0, None
  for component in nx.connected_components(graph):
    subgraph = graph.subgraph(component)
    if len(subgraph.nodes()) < 2:
      continue
    nodes_sorted = sorted(list(subgraph.nodes()), key=lambda n: n[1], reverse=True)
    top_nodes = nodes_sorted[:min(10, len(nodes_sorted))]
    bottom_nodes = nodes_sorted[-min(10, len(nodes_sorted)):]
    for upstream_node in top_nodes:
      for downstream_node in bottom_nodes:
        if upstream_node == downstream_node or upstream_node[1] <= downstream_node[1]:
          continue
        try:
          path = nx.shortest_path(subgraph, source=upstream_node, target=downstream_node, weight='weight')
          length_m = sum(subgraph[u][v]['weight'] for u, v in zip(path[:-1], path[1:]))
          if length_m > best_length:
            best_length, best_path = length_m, LineString(path)
        except nx.NetworkXNoPath:
          continue

  if best_path is None or best_length == 0:
    lines_sorted = sorted(lines, key=lambda l: l.length, reverse=True)
    best_path, best_length = lines_sorted[0], lines_sorted[0].length

  return best_length / 1000.0, best_path

def fmt_dms_lon(x, pos):
  """Matplotlib tick formatter: Web Mercator x-coordinate -> DMS longitude label."""
  pt = gpd.GeoSeries([Point(x, -967000)], crs='EPSG:3857').to_crs(epsg=4326).iloc[0]
  deg_val = abs(pt.x)
  d = int(deg_val)
  m = int((deg_val - d) * 60)
  s = int(round(((deg_val - d) * 60 - m) * 60))
  if s == 60:
    m += 1; s = 0
  if m == 60:
    d += 1; m = 0
  return f"{d}\u00b0{m}'{s}\"E"

def fmt_dms_lat(y, pos):
  """Matplotlib tick formatter: Web Mercator y-coordinate -> DMS latitude label."""
  pt = gpd.GeoSeries([Point(12820000, y)], crs='EPSG:3857').to_crs(epsg=4326).iloc[0]
  deg_val = abs(pt.y)
  d = int(deg_val)
  m = int((deg_val - d) * 60)
  s = int(round(((deg_val - d) * 60 - m) * 60))
  if s == 60:
    m += 1; s = 0
  if m == 60:
    d += 1; m = 0
  return f"{d}\u00b0{m}'{s}\"S"

def compute_abm_order(duration_hours):
  """Alternating Block Method: rank hourly rainfall blocks by placing the
  largest at the centre of the storm and alternating the remaining
  blocks in decreasing order on either side. Returns {rank: order_index}."""
  center = int(np.ceil(duration_hours / 2.0))
  abm_order = []
  left, right = center, center + 1
  place_left = True
  for _ in range(duration_hours):
    if place_left and left >= 1:
      abm_order.append(left); left -= 1; place_left = False
    elif not place_left and right <= duration_hours:
      abm_order.append(right); right += 1; place_left = True
    elif left >= 1:
      abm_order.append(left); left -= 1
    elif right <= duration_hours:
      abm_order.append(right); right += 1
  return {rank: idx + 1 for idx, rank in enumerate(abm_order)}

def compute_centroidal_distance(catchment_gdf_utm, main_channel_utm):
  """Lc = distance, measured along the main channel, from the point on
  the channel nearest the catchment centroid to the outlet (downstream
  end of the channel). Replaces the earlier approximation Lc = 0.5 * L.
  Returns None if no main-channel geometry is available (caller falls
  back to 0.5 * L in that case)."""
  if main_channel_utm is None or main_channel_utm.length == 0:
    return None
  centroid = catchment_gdf_utm.geometry.unary_union.centroid
  distance_from_upstream_end = main_channel_utm.project(centroid)
  total_length_m = main_channel_utm.length
  lc_m = total_length_m - distance_from_upstream_end
  return lc_m / 1000.0

def compute_dominant_hsg(catchment_gdf_utm, raster_path):
  """Clip the HYSOGs250m raster (pixel values 1-4=A-D, 11-14=combined
  classes) to the catchment polygon, reclassify 11-14 -> D, and return
  the dominant (majority) HSG class, its dominance percentage, and the
  full class breakdown. Falls back to ('B', 0.0, {}) if the raster is
  unavailable or does not overlap the catchment."""
  if raster_path is None:
    return 'B', 0.0, {}

  catchment_wgs84 = catchment_gdf_utm.to_crs(epsg=4326).geometry.unary_union

  with rasterio.open(raster_path) as src:
    geom_in_raster_crs = gpd.GeoSeries([catchment_wgs84], crs='EPSG:4326').to_crs(src.crs).iloc[0]
    try:
      out_image, _ = rasterio.mask.mask(src, [geom_in_raster_crs], crop=True, nodata=src.nodata)
    except ValueError:
      return 'B', 0.0, {}  # catchment polygon lies outside the raster extent
    nodata_val = src.nodata

  pixels = out_image[0]
  valid_mask = (pixels != nodata_val) if nodata_val is not None else np.ones_like(pixels, dtype=bool)
  valid_mask &= (pixels > 0)
  valid_pixels = pixels[valid_mask]

  if valid_pixels.size == 0:
    return 'B', 0.0, {}

  hsg_classes = np.array([HSG_PIXEL_MAP.get(int(v)) for v in valid_pixels])
  hsg_classes = hsg_classes[hsg_classes != None]

  if hsg_classes.size == 0:
    return 'B', 0.0, {}

  unique, counts = np.unique(hsg_classes, return_counts=True)
  total = hsg_classes.size
  class_percentages = {k: round(v / total * 100, 1) for k, v in zip(unique, counts)}
  dominant_class = max(class_percentages, key=class_percentages.get)
  dominant_percent = class_percentages[dominant_class]

  return dominant_class, dominant_percent, class_percentages

def compute_cn_from_hsg(landcover_pct, hsg_class):
  """Compute the composite CN(II) for one catchment from its land-cover
  area percentages (a dict from landcover_pct_by_catchment) and its
  dominant HSG class (A/B/C/D), using CN_TABLE_BY_HSG."""
  cn_water = (landcover_pct['water_pct'] / 100.0) * CN_TABLE_BY_HSG['water'][hsg_class]
  cn_builtup = (landcover_pct['builtup_pct'] / 100.0) * CN_TABLE_BY_HSG['builtup'][hsg_class]
  cn_forest = (landcover_pct['forest_pct'] / 100.0) * CN_TABLE_BY_HSG['forest'][hsg_class]
  cn_open = (landcover_pct['open_pct'] / 100.0) * CN_TABLE_BY_HSG['open'][hsg_class]
  cn_ii = cn_water + cn_builtup + cn_forest + cn_open
  s_ii_mm = (25400 / cn_ii) - 254 if cn_ii > 0 else 50.0
  return cn_ii, s_ii_mm

abm_rank_19 = compute_abm_order(19)
abm_rank_18 = compute_abm_order(18)

# =====================================================================
# 5. MAIN LOOP: PROCESS ALL 14 CATCHMENTS
# =====================================================================
print('--- STAGE 3: RUNNING THE 37-HOUR CONTINUOUS SIMULATION ---')

unique_gridcodes = all_catchments_gdf['gridcode'].unique()

with pd.ExcelWriter(EXCEL_OUTPUT_PATH, engine='openpyxl') as writer:
  landcover_df.to_excel(writer, sheet_name='Rekap_Tata_Guna_Lahan_CN', index=False)

  for gridcode in unique_gridcodes:
    catchment_id_str = str(int(gridcode)) if isinstance(gridcode, float) and gridcode.is_integer() else str(gridcode)
    catchment_id_val = int(gridcode) if isinstance(gridcode, float) and gridcode.is_integer() else gridcode
    print(f'-> Processing catchment: {catchment_id_str}')

    single_catchment_gdf = all_catchments_gdf[all_catchments_gdf['gridcode'] == gridcode]
    catchment_wgs84 = single_catchment_gdf.to_crs(epsg=4326)
    catchment_utm = single_catchment_gdf.to_crs(epsg=32750)

    area_km2 = catchment_utm.geometry.area.sum() / 1_000_000.0

    catchment_union_wgs84 = catchment_wgs84.geometry.unary_union
    day1_within_mask = day1_grid_result.geometry.within(catchment_union_wgs84)
    day2_within_mask = day2_grid_result.geometry.within(catchment_union_wgs84)

    if day1_within_mask.sum() > 0:
      rainfall_day1_mm = round(day1_grid_result[day1_within_mask]['ch'].mean(), 2)
      rainfall_day2_mm = round(day2_grid_result[day2_within_mask]['ch'].mean(), 2)
    else:
      catchment_centroid = catchment_wgs84.geometry.centroid.iloc[0]
      _, nearest_idx = rain_gauge_tree.query([catchment_centroid.x, catchment_centroid.y])
      rainfall_day1_mm = round(day1_values[nearest_idx], 2)
      rainfall_day2_mm = round(day2_values[nearest_idx], 2)

    catchment_union_utm = catchment_utm.geometry.unary_union
    filtered_rivers = all_rivers_gdf[all_rivers_gdf.geometry.intersects(catchment_union_utm)].copy()
    if not filtered_rivers.empty:
      filtered_rivers['geometry'] = filtered_rivers.geometry.intersection(catchment_union_utm)
      catchment_rivers = filtered_rivers[~filtered_rivers.geometry.is_empty]
    else:
      catchment_rivers = gpd.GeoDataFrame(columns=all_rivers_gdf.columns, crs='EPSG:32750')

    if "7028" in catchment_id_str:
      main_channel_km, main_channel_geom = extract_main_channel_das7028(catchment_rivers)
    else:
      main_channel_km, main_channel_geom = extract_main_channel_standard(catchment_rivers, snap_tolerance_m=0.0)

    centroidal_distance_km = compute_centroidal_distance(catchment_utm, main_channel_geom)
    if centroidal_distance_km is None or centroidal_distance_km <= 0:
      centroidal_distance_km = 0.5 * main_channel_km  # fallback if no main channel was detected

    dominant_hsg, hsg_percent, hsg_breakdown = compute_dominant_hsg(catchment_utm, hysogs_raster_path)
    hsg_label = f"Dominant HSG {dominant_hsg} ({hsg_percent:.1f}%)" if hysogs_raster_path else "Group B (assumed, HYSOGs raster not found)"

    if main_channel_geom is not None:
      main_channel_coords = list(main_channel_geom.coords)
      upstream_point = Point(main_channel_coords[0])
      outlet_point = Point(main_channel_coords[-1])
      upstream_elev = get_dem_elevation(upstream_point, dem_raster_path, label=f'upstream end, catchment {catchment_id_str}')
      outlet_elev = get_dem_elevation(outlet_point, dem_raster_path, label=f'outlet, catchment {catchment_id_str}')
      length_m = main_channel_km * 1000.0
      channel_slope = ((upstream_elev - outlet_elev) / length_m) if length_m > 0 else 0.0
      if channel_slope < 0:
        channel_slope = abs(channel_slope)
    else:
      upstream_elev, outlet_elev, channel_slope = 0.0, 0.0, 0.0

    ct = SNYDER_CT
    kirpich_tc = compute_kirpich_tc(main_channel_km, channel_slope)
    if kirpich_tc is not None and kirpich_tc > 0:
      unit_duration_tr = TR_FRACTION_OF_TC * kirpich_tc
    else:
      unit_duration_tr = 1.0  # fallback if L or slope is invalid (e.g. main channel not detected)
    lag_time_tl = ct * ((main_channel_km * centroidal_distance_km) ** 0.3) if (main_channel_km > 0 and centroidal_distance_km > 0) else 1.0
    time_to_peak_tp = lag_time_tl + 0.5 * unit_duration_tr
    time_base_tb = 5.0 * time_to_peak_tp
    peak_discharge_qp = (0.2083 * area_km2) / time_to_peak_tp if time_to_peak_tp > 0 else 0.0

    # --- Day 1 Mononobe intensity + Alternating Block Method ---
    t1 = np.arange(1, 20)
    intensity1 = (rainfall_day1_mm / 24.0) * ((24.0 / t1) ** (2.0 / 3.0))
    cumulative_depth1 = t1 * intensity1
    block_depth1 = np.zeros_like(cumulative_depth1)
    block_depth1[0] = cumulative_depth1[0]
    block_depth1[1:] = np.diff(cumulative_depth1)
    sorted_fraction1 = np.sort(block_depth1 / np.sum(block_depth1))[::-1]
    abm_fraction_day1 = np.array([sorted_fraction1[abm_rank_19[r] - 1] for r in range(1, 20)])
    hourly_rain_day1 = rainfall_day1_mm * abm_fraction_day1

    # --- Day 2 Mononobe intensity + Alternating Block Method ---
    t2 = np.arange(1, 19)
    intensity2 = (rainfall_day2_mm / 24.0) * ((24.0 / t2) ** (2.0 / 3.0))
    cumulative_depth2 = t2 * intensity2
    block_depth2 = np.zeros_like(cumulative_depth2)
    block_depth2[0] = cumulative_depth2[0]
    block_depth2[1:] = np.diff(cumulative_depth2)
    sorted_fraction2 = np.sort(block_depth2 / np.sum(block_depth2))[::-1]
    abm_fraction_day2 = np.array([sorted_fraction2[abm_rank_18[r] - 1] for r in range(1, 19)])
    hourly_rain_day2 = rainfall_day2_mm * abm_fraction_day2

    # --- Combine into the continuous 37-hour storm sequence ---
    t_continuous_37h = np.concatenate([t1, t2])
    intensity_37h = np.concatenate([intensity1, intensity2])
    cumulative_depth_37h = np.concatenate([cumulative_depth1, cumulative_depth2])
    block_depth_37h = np.concatenate([block_depth1, block_depth2])
    abm_fraction_37h = np.concatenate([abm_fraction_day1, abm_fraction_day2])
    hourly_rain_37h = np.concatenate([hourly_rain_day1, hourly_rain_day2])

    cumulative_rain_37h = np.cumsum(hourly_rain_37h)

    cn_ii_base, _ = compute_cn_from_hsg(landcover_pct_by_catchment[catchment_id_val], dominant_hsg)

    antecedent_r5_mm = compute_antecedent_r5(catchment_union_wgs84, rainfall_df)
    amc_class = classify_amc(antecedent_r5_mm)
    if antecedent_r5_mm is None:
      print(f"  [WARNING] 5-day antecedent rainfall data not found for catchment {catchment_id_str} "
            f"-> AMC assumed to be II (no CN adjustment applied)")

    # Required order: (1) base CN from the HSG table (lambda=0.2 basis) ->
    # (2) apply the AMC adjustment (the NEH-4 formulas also assume a
    # lambda=0.2 basis) -> (3) ONLY THEN convert the result to a
    # lambda=0.05 basis. Reversing this order (converting lambda first,
    # then applying AMC) would be incorrect, because the AMC formulas
    # above are not valid for a CN already on a lambda=0.05 basis.
    cn_amc_adjusted = adjust_cn_for_amc(cn_ii_base, amc_class)
    cn_final = convert_cn_to_lambda005(cn_amc_adjusted) if LAMBDA_IA == 0.05 else cn_amc_adjusted
    s_retention_mm = (25400 / cn_final) - 254 if cn_final > 0 else 50.0
    initial_abstraction_mm = LAMBDA_IA * s_retention_mm

    cumulative_peff_37h = np.where(
      cumulative_rain_37h > initial_abstraction_mm,
      ((cumulative_rain_37h - initial_abstraction_mm) ** 2) / (cumulative_rain_37h - initial_abstraction_mm + s_retention_mm),
      0.0
    )
    hourly_peff_37h = np.zeros_like(cumulative_peff_37h)
    hourly_peff_37h[0] = cumulative_peff_37h[0]
    hourly_peff_37h[1:] = np.diff(cumulative_peff_37h)
    hourly_infiltration_37h = hourly_rain_37h - hourly_peff_37h

    clock_day1 = [f'{h:02d}:00' for h in range(5, 24)]
    clock_day2 = [f'{h:02d}:00' for h in range(0, 18)]
    clock_time_list = clock_day1 + clock_day2

    hyetograph_table = pd.DataFrame({
        't Kontinu (jam)': np.arange(1, 38),
        'Tanggal Event': ['9 Sept 2025'] * 19 + ['10 Sept 2025'] * 18,
        'Jam Kejadian': clock_time_list,
        't (relatif)': t_continuous_37h,
        'Ir (mm/jam)': np.round(intensity_37h, 2),
        'Ir * Td': np.round(cumulative_depth_37h, 2),
        '\u0394p (mm)': np.round(block_depth_37h, 2),
        'pt (ABM)': [f'{v*100:.2f}%' for v in abm_fraction_37h],
        'P (mm)': np.round(hourly_rain_37h, 2),
        '\u03a3P (mm)': np.round(cumulative_rain_37h, 2),
        'S (mm)': np.round(s_retention_mm, 2),
        'Ia (mm)': np.round(initial_abstraction_mm, 2),
        '\u03a3 Peff (mm)': np.round(cumulative_peff_37h, 2),
        'Peff (mm)': np.round(hourly_peff_37h, 2),
        'Infil (mm)': np.round(hourly_infiltration_37h, 2)
    })

    # SCS dimensionless unit hydrograph (standard curvilinear shape)
    dimensionless_t_ratio = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0])
    dimensionless_q_ratio = np.array([0.0, 0.03, 0.1, 0.19, 0.31, 0.47, 0.66, 0.82, 0.93, 0.99, 1.0, 0.99, 0.93, 0.86, 0.78, 0.68, 0.56, 0.46, 0.39, 0.33, 0.28, 0.2435, 0.207, 0.177, 0.147, 0.127, 0.107, 0.092, 0.077, 0.066, 0.055, 0.0475, 0.04, 0.0345, 0.029, 0.025, 0.021, 0.018, 0.015, 0.013, 0.011, 0.0098, 0.0086, 0.0074, 0.0062, 0.005, 0.004, 0.003, 0.002, 0.001, 0.0])

    unit_hydrograph_t = dimensionless_t_ratio * time_to_peak_tp
    unit_hydrograph_q = dimensionless_q_ratio * peak_discharge_qp
    max_interp_hour = int(np.ceil(time_base_tb)) + 37
    hourly_time_axis = np.arange(0, max_interp_hour + 1)
    unit_hydrograph_hourly = np.interp(hourly_time_axis, unit_hydrograph_t, unit_hydrograph_q, right=0.0)

    num_rows = len(unit_hydrograph_hourly)
    num_rain_hours = 37
    convolution_matrix = np.zeros((num_rows, num_rain_hours))

    for j in range(num_rain_hours):
      peff_value = hourly_peff_37h[j]
      if j < num_rows:
        length_to_copy = min(len(unit_hydrograph_hourly), num_rows - j)
        convolution_matrix[j:j + length_to_copy, j] = unit_hydrograph_hourly[:length_to_copy] * peff_value

    total_flood_discharge = np.sum(convolution_matrix, axis=1)
    convolution_dict = {'t (jam)': hourly_time_axis}
    for j in range(num_rain_hours):
      convolution_dict[f'Jam-{j+1} ({hourly_peff_37h[j]:.2f}mm)'] = np.round(convolution_matrix[:, j], 3)
    convolution_dict['Debit Total (m\u00b3/s)'] = np.round(total_flood_discharge, 2)
    convolution_table = pd.DataFrame(convolution_dict)

    # --- Location map for this catchment (saved into the Excel sheet) ---
    fig, ax = plt.subplots(figsize=(6, 8), dpi=200)

    all_catchments_web = all_catchments_gdf.to_crs(epsg=3857)
    focus_catchment_web = all_catchments_web[all_catchments_web['gridcode'] == gridcode]
    other_catchments_web = all_catchments_web[all_catchments_web['gridcode'] != gridcode]

    other_catchments_web.plot(ax=ax, facecolor='#D1C4E9', edgecolor='black', linewidth=0.6, alpha=0.5, label='CA')
    focus_catchment_web.plot(ax=ax, facecolor='#FFD700', edgecolor='black', linewidth=1.2, alpha=0.85, label=f'CA - {catchment_id_str}')

    if not catchment_rivers.empty:
      rivers_web = catchment_rivers.to_crs(epsg=3857)
      rivers_web.plot(ax=ax, color='#1E88E5', linewidth=0.8, alpha=0.7, label='River')

    if main_channel_geom is not None:
      main_channel_web = gpd.GeoSeries([main_channel_geom], crs='EPSG:32750').to_crs(epsg=3857).iloc[0]
      gpd.GeoSeries([main_channel_web], crs='EPSG:3857').plot(ax=ax, color='red', linewidth=2.2, label='Main River')
      web_coords = list(main_channel_web.coords)
      ax.plot(web_coords[0][0], web_coords[0][1], 'go', markersize=5)
      ax.plot(web_coords[-1][0], web_coords[-1][1], 'mo', markersize=5)

    try:
      cx.add_basemap(ax, source=cx.providers.Esri.WorldImagery, zoom=12, attribution="")
    except Exception:
      pass

    focus_bounds = focus_catchment_web.total_bounds
    margin_x = (focus_bounds[2] - focus_bounds[0]) * 0.45
    margin_y = (focus_bounds[3] - focus_bounds[1]) * 0.45
    ax.set_xlim(focus_bounds[0] - margin_x, focus_bounds[2] + margin_x)
    ax.set_ylim(focus_bounds[1] - margin_y, focus_bounds[3] + margin_y)

    ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=2))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=4))

    ax.xaxis.set_major_formatter(FuncFormatter(fmt_dms_lon))
    ax.yaxis.set_major_formatter(FuncFormatter(fmt_dms_lat))

    ax.tick_params(axis='both', which='major', labelsize=6, labeltop=True, labelbottom=True, labelleft=True, labelright=True)
    ax.grid(True, linestyle=':', color='gray', alpha=0.5)

    ax.legend(loc='lower right', frameon=True, facecolor='white', framealpha=0.9, fontsize=7)

    esri_attribution = "Tiles (C) Esri \u2014 Source: Esri, i-cubed, USDA, USGS, AEX, GeoEye, Getmapping, Aerogrid, IGN, IGP, UPR-EGP, and the GIS User Community"
    ax.text(0.02, 0.02, esri_attribution, transform=ax.transAxes, fontsize=3.5, fontweight='bold',
            color='white', bbox=dict(boxstyle="square,pad=0.2", facecolor="black", alpha=0.5, edgecolor="none"),
            wrap=True, verticalalignment='bottom')

    img_buffer = io.BytesIO()
    plt.savefig(img_buffer, format='png', bbox_inches='tight', dpi=150)
    img_buffer.seek(0)
    plt.close(fig)

    # --- Write the summary sheet for this catchment ---
    sheet_name = f"DAS_{catchment_id_str}"
    summary_table = pd.DataFrame({
        'Parameter Hidrologi & Geometri': [
            'Gridcode DAS', 'Hydrologic Soil Group (HYSOGs250m)',
            'Curve Number (CN II, basis lambda=0.2)',
            'Hujan 5 Hari Antesenden (R5) [mm]', 'Kelas AMC',
            f'Curve Number Final (basis lambda={LAMBDA_IA})',
            'Luas DAS (A) [km\u00b2]',
            'Panjang Sungai Utama (L) [km]', 'Jarak Centroid-Outlet (Lc) [km]',
            'Elevasi Hulu Sungai (Z_hulu) [m]', 'Elevasi Hilir Sungai (Z_hilir) [m]',
            'Kemiringan Saluran (S0) [m/m]', 'Kemiringan Saluran (S0) [%]',
            'Hujan Harian 9 Sept (R1) [mm]', 'Hujan Harian 10 Sept (R2) [mm]',
            'Total Durasi Kejadian [jam]',
            'Time of Concentration Kirpich (Tc) [jam]', 'Durasi Hujan Satuan (Tr) [jam]',
            'Time Lag (TL) [jam]', 'Time Peak (TP) [jam]',
            'Time Base (TB) [jam]', 'Debit Puncak Hidrograf Satuan (QP) [m\u00b3/s/mm]',
            'DEBIT PUNCAK BANJIR TOTAL (KONTINU) [m\u00b3/s]'
        ],
        'Nilai': [
            catchment_id_str, hsg_label, round(cn_ii_base, 2),
            'n/a' if antecedent_r5_mm is None else round(antecedent_r5_mm, 2), amc_class,
            round(cn_final, 2),
            round(area_km2, 2), round(main_channel_km, 2), round(centroidal_distance_km, 2),
            round(upstream_elev, 2), round(outlet_elev, 2), round(channel_slope, 5), round(channel_slope * 100, 3),
            round(rainfall_day1_mm, 2), round(rainfall_day2_mm, 2), 37,
            'n/a' if kirpich_tc is None else round(kirpich_tc, 3), round(unit_duration_tr, 3),
            round(lag_time_tl, 3), round(time_to_peak_tp, 3), round(time_base_tb, 3), round(peak_discharge_qp, 3),
            round(np.max(total_flood_discharge), 2)
        ]
    })

    summary_table.to_excel(writer, sheet_name=sheet_name, index=False, startrow=0)

    hyetograph_start_row = len(summary_table) + 3
    worksheet = writer.sheets[sheet_name]
    worksheet.cell(row=hyetograph_start_row, column=1, value="TABEL INTENSITAS HUJAN 37 JAM (MONONOBE, ABM & EKSTENSI SCS-CN)")
    hyetograph_table.to_excel(writer, sheet_name=sheet_name, index=False, startrow=hyetograph_start_row)

    convolution_start_row = hyetograph_start_row + len(hyetograph_table) + 3
    worksheet.cell(row=convolution_start_row, column=1, value="TABEL SUPERPOSISI KONVOLUSI HIDROGRAF DEBIT BANJIR (37 JAM KONTINU)")
    convolution_table.to_excel(writer, sheet_name=sheet_name, index=False, startrow=convolution_start_row)

    temp_image_path = f"/tmp/plot_das_{catchment_id_str}.png"
    with open(temp_image_path, "wb") as f:
      f.write(img_buffer.getbuffer())

    excel_image = OpenpyxlImage(temp_image_path)
    excel_image.width, excel_image.height = 340, 440
    worksheet.add_image(excel_image, "AU2")

print('\n=====================================================================')
print(f'PROCESSING COMPLETE. OUTPUT FILE:\n{EXCEL_OUTPUT_PATH}')
print('=====================================================================')
