================================================================================
SOUTH BALI FLOOD HYDROGRAPH MODELLING
================================================================================

Python / Google Colab pipeline for reconstructing the design flood hydrograph
of the 9-10 September 2025 extreme rainfall event across 14 sub-catchments in
Denpasar, Badung, and Gianyar, South Bali, Indonesia. Developed for a
companion manuscript on comparative hydrodynamic simulation and field-based
flood validation.


--------------------------------------------------------------------------------
WHAT'S IN THIS REPOSITORY
--------------------------------------------------------------------------------

  catchment_and_river_delineation.py
      DEM -> watershed and river network delineation (WhiteboxTools D8
      pipeline). Run FIRST.

  hydrograph_model.py
      Full hydrologic model: rainfall analysis, curve number, Snyder unit
      hydrograph, hourly convolution -> design flood hydrograph per
      catchment. Run SECOND, using the output of the first script.

Both scripts are written to run in Google Colab against a Google Drive
folder structure (see INPUT DATA below). They are not currently packaged
as a pip-installable library; copy the cell contents into a Colab notebook,
or adapt the Drive paths at the top of each file to run elsewhere.


--------------------------------------------------------------------------------
METHODOLOGY SUMMARY
--------------------------------------------------------------------------------

  1. Terrain processing (catchment_and_river_delineation.py)
     Breach depressions, D8 flow direction, automatic basin delineation
     from a national DEM, then merge any catchment smaller than 10 km^2
     into its largest neighbour.

  2. Land cover & soil
     Land-cover class shares per catchment (expected as a pre-computed
     CSV, typically from a Google Earth Engine Random Forest
     classification) combined with the dominant Hydrologic Soil Group
     per catchment, extracted directly from the HYSOGs250m global raster
     (Ross et al., 2018). DOI: 10.1038/sdata.2018.91

  3. Curve Number
     NRCS TR-55/NEH-630 table by land-cover class x HSG, adjusted for
     Antecedent Moisture Condition (5-day antecedent rainfall vs. a
     seasonal threshold table) and converted from the conventional
     initial-abstraction ratio (lambda = Ia/S = 0.2) to lambda = 0.05
     following Hawkins et al. (2002).
     Source: https://ponce.sdsu.edu/hawkins_initial_abstraction.pdf

  4. Design storm
     Mononobe intensity-duration disaggregation of daily rainfall depth,
     reordered by the Alternating Block Method, into a continuous
     37-hour hyetograph spanning two calendar days without resetting
     cumulative infiltration at midnight.

  5. Unit hydrograph
     Snyder method (Ct = 2.0, within Snyder 1938's typical range for
     ungauged basins), with the centroidal distance (Lc) computed as the
     true along-channel distance from centroid to outlet, and the unit
     storm duration (Tr) derived from a Kirpich (1940) time-of-
     concentration estimate.

  6. Convolution
     Hourly effective rainfall convolved with the SCS dimensionless unit
     hydrograph to produce the continuous design flood hydrograph per
     catchment.

Full citations and the underlying reasoning for each methodological choice
are documented in the companion manuscript and inline in the code comments.


--------------------------------------------------------------------------------
REQUIREMENTS
--------------------------------------------------------------------------------

  geopandas
  rasterio
  shapely
  scipy
  numpy
  pandas
  networkx
  matplotlib
  openpyxl
  contextily
  whitebox          (catchment_and_river_delineation.py only)

Both scripts auto-install rasterio / contextily / whitebox if missing when
run in Colab. A Google account with Drive access is required
(google.colab.drive.mount).


--------------------------------------------------------------------------------
INPUT DATA
--------------------------------------------------------------------------------

Both scripts expect a fixed Google Drive folder layout (folder names are
matched literally; adjust the path constants near the top of each script
if yours differ):

  DEM BALI BARU.tif
      Raw DEM covering the study area (input to the delineation script).

  DAS Phyton - Tahun 1/
      Catchment boundary shapefile(s) (output of the delineation script,
      or your own) -- must include a 'gridcode' field.

  Sungai Phyton - Tahun 1/
      River network line shapefile(s).

  Bali Elevation - Tahun 1/
      DEM raster (.tif / .dem) used for elevation/slope extraction.

  Jenis Tanah/
      HYSOGs250m raster tile(s) covering the study area.

  Curah Hujan Bali Banjir September 2025.csv
      Rainfall grid. Columns: lon, lat, tgl1 (date, any common format),
      ch (rainfall depth, mm). Must include at least the two event days
      and, for the AMC calculation, the 5 preceding days.

  Tata Guna Lahan & Manning 2025.csv
      Land-cover area per catchment. Columns: Gridcode, Air_m2,
      Bangunan_m2, Hutan_m2, Terbuka_m2 (water / built-up / forest /
      open, in square metres).


--------------------------------------------------------------------------------
RUNNING THE PIPELINE
--------------------------------------------------------------------------------

  1. Open catchment_and_river_delineation.py in Colab, run all cells.
     -> produces Batas_DAS_Bali_Cleaned_10km2.shp in your Drive root.

  2. Move/rename that output (or your own catchment shapefile) into the
     'DAS Phyton - Tahun 1' folder referenced by hydrograph_model.py.

  3. Open hydrograph_model.py in Colab, run all cells.
     -> produces Hasil_Analisis_Hidrologi_14_DAS_Master_Final_Manning2025.xlsx,
        one sheet per catchment (hourly hyetograph, infiltration, and full
        37-hour convolution table) plus a land-cover/CN summary sheet.

hydrograph_model.py is defensive about common data issues: it validates
the rainfall CSV's date range before proceeding (raising a clear error
listing the dates actually found if the expected event dates are
missing), and prints an explicit warning -- rather than failing silently
-- whenever a DEM read, HSG raster lookup, or antecedent-rainfall lookup
falls back to a default value.


--------------------------------------------------------------------------------
OUTPUT
--------------------------------------------------------------------------------

hydrograph_model.py writes a single Excel workbook containing:

  - Rekap_Tata_Guna_Lahan_CN
      Land-cover percentages per catchment.

  - One DAS_<gridcode> sheet per catchment, containing:
      * A parameter summary (area, main channel length, centroidal
        distance, slope, HSG, CN before/after AMC and lambda adjustment,
        Snyder/Kirpich parameters, peak discharge).
      * The full 37-hour hyetograph and infiltration table.
      * The full hourly convolution (superposition) table and resulting
        design flood hydrograph.
      * A location map of the catchment within the study area.


--------------------------------------------------------------------------------
KNOWN LIMITATIONS
--------------------------------------------------------------------------------

  - Parameter values (Snyder Ct, Kirpich coefficients, Manning's n, NRCS
    curve numbers) are drawn from published literature developed outside
    Indonesia's tropical, volcanic-derived terrain, in the absence of
    local calibration data (see the companion manuscript's limitations
    section).

  - The main-channel extraction uses two separate algorithms for one
    geometrically irregular catchment (7028) versus the other thirteen;
    see the docstrings in hydrograph_model.py before modifying either.

  - This repository produces the HYDROLOGIC design hydrograph only.
    Hydrodynamic (2D inundation) simulation and field validation are
    handled downstream in HEC-RAS and are not part of this repository.


--------------------------------------------------------------------------------
CITATION
--------------------------------------------------------------------------------

If you use this code, please cite the companion manuscript (citation
details to be added upon publication) and the primary data/method
sources listed in the Methodology Summary section above.


--------------------------------------------------------------------------------
LICENSE
--------------------------------------------------------------------------------

Specify a license (e.g. MIT) before making this repository public.
