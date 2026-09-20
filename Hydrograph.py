# =====================================================================
# MASTER SCRIPT: PEMODELAN BANJIR 14 DAS BALI (SEPTEMBER 2025)
# UPDATE FILE: 'Tata Guna Lahan & Manning 2025.csv'
# =====================================================================
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
# 1. MOUNT DRIVE & SETUP PATH
# =====================================================================
drive.mount('/content/drive')

PATH_DRIVE = '/content/drive/MyDrive/'
FILE_CSV_HUJAN = os.path.join(PATH_DRIVE, 'Curah Hujan Bali Banjir September 2025.csv')
FILE_CSV_TATA_GUNA = os.path.join(PATH_DRIVE, 'Tata Guna Lahan & Manning 2025.csv') # UPDATE NAMA FILE
PATH_DAS = os.path.join(PATH_DRIVE, 'DAS Phyton - Tahun 1')
PATH_SUNGAI = os.path.join(PATH_DRIVE, 'Sungai Phyton - Tahun 1')
PATH_DEM_FOLDER = os.path.join(PATH_DRIVE, 'Bali Elevation - Tahun 1')
PATH_JENIS_TANAH = os.path.join(PATH_DRIVE, 'Jenis Tanah')  # folder raster HYSOGs250m
EXCEL_OUTPUT_PATH = os.path.join(PATH_DRIVE, 'Hasil_Analisis_Hidrologi_14_DAS_Master_Final_Manning2025.xlsx')

list_dem_files = glob.glob(os.path.join(PATH_DEM_FOLDER, '**/*.tif'), recursive=True) + \
                 glob.glob(os.path.join(PATH_DEM_FOLDER, '**/*.dem'), recursive=True)

dem_raster_path = list_dem_files[0] if list_dem_files else None

list_hysogs_files = glob.glob(os.path.join(PATH_JENIS_TANAH, '**/*.tif'), recursive=True)
hysogs_raster_path = list_hysogs_files[0] if list_hysogs_files else None
if hysogs_raster_path is None:
  print("PERINGATAN: raster HYSOGs250m tidak ditemukan di folder 'Jenis Tanah'. "
        "HSG akan fallback ke asumsi Group B untuk semua DAS.")

df_hujan_grid = pd.read_csv(FILE_CSV_HUJAN)

# --- Parse kolom tgl1 jadi objek tanggal ASLI (bukan cocokkan teks).
#     pd.to_datetime otomatis mengenali berbagai format ('9/9/2025',
#     '2025-09-09', '09/09/2025', dst) tanpa perlu tahu formatnya lebih
#     dulu -- jauh lebih tahan banting daripada str.startswith(). ---
df_hujan_grid['tgl1'] = pd.to_datetime(df_hujan_grid['tgl1'], errors='coerce')

# Buang baris tanpa koordinat/nilai hujan/tanggal yang valid
n_sebelum = len(df_hujan_grid)
df_hujan_grid = df_hujan_grid.dropna(subset=['lon', 'lat', 'ch', 'tgl1'])
if len(df_hujan_grid) < n_sebelum:
  print(f"PERINGATAN: {n_sebelum - len(df_hujan_grid)} baris dibuang dari "
        f"{FILE_CSV_HUJAN} karena lon/lat/ch/tgl1 kosong atau tanggal "
        f"tidak terbaca.")

TANGGAL_EVENT_HARI1 = pd.Timestamp('2025-09-09')
TANGGAL_EVENT_HARI2 = pd.Timestamp('2025-09-10')

df_hujan_day1 = df_hujan_grid[df_hujan_grid['tgl1'].dt.normalize() == TANGGAL_EVENT_HARI1].copy()
df_hujan_day2 = df_hujan_grid[df_hujan_grid['tgl1'].dt.normalize() == TANGGAL_EVENT_HARI2].copy()

if df_hujan_day1.empty or df_hujan_day2.empty:
  contoh_tgl1 = sorted(df_hujan_grid['tgl1'].dt.strftime('%Y-%m-%d').unique())[:15]
  raise ValueError(
    "GAGAL MEMUAT DATA HUJAN: tidak ada baris bertanggal "
    f"{TANGGAL_EVENT_HARI1.date()} atau {TANGGAL_EVENT_HARI2.date()} di "
    f"{FILE_CSV_HUJAN}.\n"
    f"  Total baris di file      : {len(df_hujan_grid)}\n"
    f"  Baris cocok hari 1       : {len(df_hujan_day1)}\n"
    f"  Baris cocok hari 2       : {len(df_hujan_day2)}\n"
    f"  Tanggal yang ADA di file : {contoh_tgl1}\n"
    "  -> Bandingkan tanggal di atas dengan TANGGAL_EVENT_HARI1/2 di kode "
    "ini; kalau event flood kamu bukan 9-10 September, kabari Claude "
    "tanggal yang benar.")

coords_grid = df_hujan_day1[['lon', 'lat']].values
tree_hujan = cKDTree(coords_grid)

# =====================================================================
# 1b. KONFIGURASI PARAMETER (kumpulan nilai yang paling sering di-tuning)
# =====================================================================
# --- Tabel CN per kelas tata guna lahan x HSG (NRCS TR-55 / NEH-630
#     Ch.9, USDA 1984). Kolom B dari tabel ini identik dengan konstanta
#     CN_AIR/BANGUNAN/TERBUKA/HUTAN versi lama (100/85/61/55) -- jadi
#     baseline lama otomatis konsisten begitu HSG dominan = 'B'.
#     Baris sumber: Hutan="Woods, good condition"; Terbuka="Pasture,
#     good condition"; Bangunan="Residential 1/8 acre lot, 65% impervious";
#     Air=100 untuk semua HSG (badan air tidak berinfiltrasi). ---
TABEL_CN_HSG = {
  'Air':      {'A': 100, 'B': 100, 'C': 100, 'D': 100},
  'Hutan':    {'A': 30,  'B': 55,  'C': 70,  'D': 77},
  'Terbuka':  {'A': 39,  'B': 61,  'C': 74,  'D': 80},
  'Bangunan': {'A': 77,  'B': 85,  'C': 90,  'D': 92},
}

# --- Reklasifikasi pixel raster HYSOGs250m -> kelas HSG A-D.
#     Kelas gabungan (dual/undrained) 11=A/D, 12=B/D, 13=C/D, 14=D/D
#     seluruhnya dikonversi jadi D (skenario konservatif/undrained). ---
MAP_HSG_PIXEL = {1: 'A', 2: 'B', 3: 'C', 4: 'D', 11: 'D', 12: 'D', 13: 'D', 14: 'D'}

# --- Parameter Snyder Unit Hydrograph: TL = Ct * (L * Lc)^0.3 ---
# Ct  : koefisien basin, rentang literatur 1.8-2.2 (rumus tanpa faktor
#       konversi C1 SI; lihat catatan diskusi metodologi)
CT_SNYDER = 2.0

# --- Kirpich (1940) time of concentration -> dipakai untuk turunkan Tr
#     (durasi hujan satuan), BUKAN untuk TL Snyder. Valid untuk DAS kecil
#     berlereng jelas (Roussel et al. 2005: 0.25-150 mil^2, S=0.002-0.1);
#     14 DAS Bali (1.1-102 km^2) masuk rentang ini, meski akurasi lebih
#     rendah untuk sub-DAS pesisir berlereng sangat landai.
#     Tc [menit] = 0.06628 * L[km]^0.77 / S[m/m]^0.385  (versi SI Kirpich,
#     persis konstanta di slide kamu; L langsung km, tidak perlu ke mil)
#     Tr = FRAKSI_TR_DARI_TC * Tc  (rentang literatur umum 0.1-0.3;
#     dipilih 0.1 -- batas bawah -- karena hampir sama dengan konvensi
#     Snyder tr=TL/5.5 dengan TL~0.6*Tc, memberi ~0.109*Tc)
FRAKSI_TR_DARI_TC = 0.1

# --- Pilihan kolom musim untuk klasifikasi AMC (Tabel 6.6). DIPUTUSKAN:
#     'kering', berdasarkan BMKG (2025), "Prakiraan Musim Kemarau 2025 di
#     Indonesia", Direktorat Perubahan Iklim, Deputi Bidang Klimatologi,
#     Jakarta, Maret 2025 -- zona BALI_19 (cakupan: Badung bagian selatan,
#     Gianyar bagian selatan, Tabanan bagian selatan, Kota Denpasar --
#     PERSIS area studi ini). Musim kemarau zona ini: awal April III
#     (mundur 1 dasarian dari normal), panjang 19 dasarian (lebih pendek 1
#     dasarian dari normal) -> berlangsung ±21 April s/d ±31 Oktober 2025.
#     Event 9-10 September jatuh di dasarian ke-14/15 dari 19 -- MASIH DI
#     DALAM musim kemarau resmi BMKG untuk zona ini, bukan musim hujan. ---
MUSIM_AMC = 'kering'  # 'semi' atau 'kering'

# Ambang Tabel 6.6 (jumlah hujan 5 hari terdahulu, cm -> mm)
AMBANG_AMC_MM = {
  'kering': {'I': 13.0, 'III': 28.0},   # <1.3 cm ; >2.8 cm
  'semi':   {'I': 36.0, 'III': 53.0},   # <3.6 cm ; >5.3 cm
}

def klasifikasi_amc(R5_mm, musim=MUSIM_AMC):
  """Klasifikasi AMC dari total hujan 5 hari antesenden [mm], memakai
  ambang Tabel 6.6 (textbook hidrologi terapan, hal. 157) sesuai kolom
  musim yang dipilih lewat MUSIM_AMC (lihat catatan sitasi BMKG di atas)."""
  if R5_mm is None:
    return 'II'  # fallback: tidak ada data antesenden -> asumsi kondisi normal
  batas = AMBANG_AMC_MM[musim]
  if R5_mm < batas['I']: return 'I'
  if R5_mm <= batas['III']: return 'II'
  return 'III'

def adjust_cn_amc(CN_II, kelas_amc):
  """Konversi CN(II) -> CN sesuai kelas AMC, formula standar SCS/NEH-4
  (Chow, Maidment & Mays, 1988, Applied Hydrology, Ch.5)."""
  if kelas_amc == 'I':
    return (4.2 * CN_II) / (10 - 0.058 * CN_II)
  elif kelas_amc == 'III':
    return (23 * CN_II) / (10 + 0.13 * CN_II)
  return CN_II  # AMC II -> tidak berubah

def hitung_r5_antesenden(geom_das_wgs84_union, df_hujan_grid_full, tanggal_mulai_event=TANGGAL_EVENT_HARI1):
  """Total hujan 5 hari sebelum event (memakai kolom 'tgl1' yang sudah
  di-parse jadi objek tanggal asli di Tahap 1).

  PENTING (perbaikan konsistensi metode): dulu fungsi ini mengecek
  apakah titik grid ASLI (jarang, ~11 km antar titik) jatuh persis di
  dalam poligon DAS -- untuk DAS kecil/posisi tertentu, kebetulan tidak
  ada satu pun titik asli yang jatuh di dalamnya, sehingga hasilnya
  'n/a' padahal datanya ADA di CSV. R_day1/R_day2 tidak punya masalah
  ini karena keduanya diinterpolasi dulu ke mesh rapat 400x400 sebelum
  dicek per DAS. Fungsi ini sekarang memakai pendekatan yang SAMA:
  interpolasi ke mesh global gdf_grid_coords/xx/yy (sudah dibuat di
  Tahap 1), baru dicek titik mana yang masuk DAS."""
  tgl_awal = tanggal_mulai_event - pd.Timedelta(days=5)
  tgl_akhir = tanggal_mulai_event - pd.Timedelta(days=1)

  mask_ada = (df_hujan_grid_full['tgl1'].dt.normalize() >= tgl_awal) & \
             (df_hujan_grid_full['tgl1'].dt.normalize() <= tgl_akhir)
  if mask_ada.sum() == 0:
    return None  # data 5 hari sebelum event belum ada di CSV

  # Jumlahkan hujan 5 hari PER titik grid ASLI dulu (total R5 per idgrid)
  df_r5 = df_hujan_grid_full[mask_ada]
  total_per_titik_asli = df_r5.groupby(['lon', 'lat'])['ch'].sum().reset_index()
  if len(total_per_titik_asli) < 3:
    return None  # terlalu sedikit titik untuk interpolasi linear yang berarti

  coords_r5 = total_per_titik_asli[['lon', 'lat']].values
  val_r5 = total_per_titik_asli['ch'].values

  # Interpolasi ke mesh global yang SAMA dipakai R_day1/R_day2 (variabel
  # gdf_grid_coords, xx, yy dari Tahap 1) -- konsisten metodenya
  grid_r5_lin = griddata(coords_r5, val_r5, (xx, yy), method='linear')
  grid_r5_near = griddata(coords_r5, val_r5, (xx, yy), method='nearest')
  gdf_grid_res_r5 = gdf_grid_coords.copy()
  gdf_grid_res_r5['ch'] = np.where(np.isnan(grid_r5_lin), grid_r5_near, grid_r5_lin).ravel()

  dalam_das = gdf_grid_res_r5[gdf_grid_res_r5.geometry.within(geom_das_wgs84_union)]
  if dalam_das.empty:
    # fallback terakhir: titik mesh terdekat ke centroid DAS (harusnya
    # nyaris tidak pernah kepakai, mesh 400x400 sangat rapat)
    centroid_das = geom_das_wgs84_union.centroid
    coords_mesh = np.vstack([gdf_grid_res_r5.geometry.x, gdf_grid_res_r5.geometry.y]).T
    tree_mesh = cKDTree(coords_mesh)
    _, idx_terdekat = tree_mesh.query([centroid_das.x, centroid_das.y])
    return float(gdf_grid_res_r5['ch'].iloc[idx_terdekat])

  return float(dalam_das['ch'].mean())


def hitung_tc_kirpich(L_km, slope_frac):
  """Kirpich (1940) time of concentration, versi SI dari slide referensi.
  Tc [JAM] = 0.06628 * L[km]^0.77 / S[m/m]^0.385
  PENTING (koreksi bug): dengan koefisien 0.06628 dan L dalam km, rumus
  ini menghasilkan Tc LANGSUNG DALAM JAM -- bukan menit. Versi sebelumnya
  salah membagi hasil dengan 60, membuat Tc 60x lebih kecil dari
  seharusnya (mis. sungai 43 km jadi Tc~5 menit, padahal seharusnya ~5
  jam). Cross-check: untuk L=1km, S=1% (0,01), formula ini memberi
  Tc=0,39 jam (~23 menit) -- masuk akal untuk DAS kecil berlereng
  curam; kalau dianggap menit (0,39 menit = 23 detik) jelas tidak
  masuk akal secara fisik."""
  if L_km <= 0 or slope_frac <= 0:
    return None
  Tc_jam = 0.06628 * (L_km ** 0.77) / (slope_frac ** 0.385)
  return Tc_jam

# --- Rasio abstraksi awal Ia/S (lambda). Default historis SCS 1956: 0.2.
#     Hawkins et al. (2002) merekomendasikan 0.05 berdasarkan 307 DAS/28.301
#     event. Set ke 0.2 untuk kembali ke perilaku lama. ---
LAMBDA_IA = 0.05

def cn_ke_lambda_005(CN_020):
  """Konversi CN yang dikalibrasi pada basis lambda=0.2 (standar tabel NRCS,
  termasuk TABEL_CN_HSG di atas) menjadi CN ekuivalen pada basis lambda=0.05.
  Persamaan 9, Hawkins, R.H., Jiang, R., Woodward, D.E., Hjelmfelt, A.T.,
  Van Mullem, J.A. (2002). "Runoff Curve Number Method: Examination of the
  Initial Abstraction Ratio." Full text (diverifikasi langsung dari PDF):
  https://ponce.sdsu.edu/hawkins_initial_abstraction.pdf
  CN0.05 = 100 / (1.879*(100/CN0.20 - 1)^1.15 + 1)
  Berlaku hingga CN0.20~98.5; di atas itu CN0.05~CN0.20 (S->0)."""
  if CN_020 >= 100:
    return 100.0
  return 100.0 / (1.879 * ((100.0 / CN_020 - 1) ** 1.15) + 1)

# =====================================================================
# 2. PROSES TATA GUNA LAHAN (persentase area per kelas)
#    CN_II & S_II_mm TIDAK dihitung di sini lagi -- sekarang tergantung
#    HSG dominan per-DAS (hasil raster HYSOGs250m), jadi dipindah ke
#    dalam loop iterasi DAS (Tahap 5).
# =====================================================================
list_gridcode = [5931, 6310, 6553, 6653, 6799, 7028, 7278, 7498, 7514, 8489, 10247, 11445, 11905, 13106]

df_luas = pd.read_csv(FILE_CSV_TATA_GUNA)
if 'Gridcode' not in df_luas.columns:
  df_luas['Gridcode'] = list_gridcode

# Konversi kolom m2 ke float jika ada
for col in ['Air_m2', 'Bangunan_m2', 'Hutan_m2', 'Terbuka_m2']:
  if col in df_luas.columns:
    df_luas[col] = df_luas[col].astype(np.float64)

df_luas['Total_m2'] = df_luas[['Air_m2', 'Bangunan_m2', 'Hutan_m2', 'Terbuka_m2']].sum(axis=1)

df_luas['Air_%'] = (df_luas['Air_m2'] / df_luas['Total_m2']) * 100
df_luas['Bangunan_%'] = (df_luas['Bangunan_m2'] / df_luas['Total_m2']) * 100
df_luas['Hutan_%'] = (df_luas['Hutan_m2'] / df_luas['Total_m2']) * 100
df_luas['Terbuka_%'] = (df_luas['Terbuka_m2'] / df_luas['Total_m2']) * 100
dict_persen_lahan = df_luas.set_index('Gridcode')[['Air_%', 'Bangunan_%', 'Hutan_%', 'Terbuka_%']].to_dict('index')

# =====================================================================
# 3. GEOSPASIAL DATA SHAPEFILE & GRID HUJAN
# =====================================================================
list_shp_das = glob.glob(os.path.join(PATH_DAS, '**/*.shp'), recursive=True)
list_shp_sungai = glob.glob(os.path.join(PATH_SUNGAI, '**/*.shp'), recursive=True)

gdfs_das = [gpd.read_file(p) for p in list_shp_das]
gdf_das_all = gpd.GeoDataFrame(pd.concat(gdfs_das, ignore_index=True), crs=gdfs_das[0].crs)

gdfs_s = [gpd.read_file(p).to_crs(epsg=32750) for p in list_shp_sungai]
gdf_sungai_all = gpd.GeoDataFrame(pd.concat(gdfs_s, ignore_index=True), crs='EPSG:32750')

gdf_das_wgs84_total = gdf_das_all.to_crs(epsg=4326)
xmin, ymin, xmax, ymax = gdf_das_wgs84_total.total_bounds
x = np.linspace(xmin, xmax, 400)
y = np.linspace(ymin, ymax, 400)
xx, yy = np.meshgrid(x, y)
df_grid_coords = pd.DataFrame(np.vstack([xx.ravel(), yy.ravel()]).T, columns=['lon', 'lat'])
gdf_grid_coords = gpd.GeoDataFrame(df_grid_coords, geometry=[Point(xyz) for xyz in zip(df_grid_coords['lon'], df_grid_coords['lat'])], crs='EPSG:4326')

val_day1 = df_hujan_day1['ch'].values
val_day2 = df_hujan_day2['ch'].values

grid_day1_lin = griddata(coords_grid, val_day1, (xx, yy), method='linear')
grid_day1_near = griddata(coords_grid, val_day1, (xx, yy), method='nearest')
gdf_grid_res_day1 = gdf_grid_coords.copy()
gdf_grid_res_day1['ch'] = np.where(np.isnan(grid_day1_lin), grid_day1_near, grid_day1_lin).ravel()

grid_day2_lin = griddata(coords_grid, val_day2, (xx, yy), method='linear')
grid_day2_near = griddata(coords_grid, val_day2, (xx, yy), method='nearest')
gdf_grid_res_day2 = gdf_grid_coords.copy()
gdf_grid_res_day2['ch'] = np.where(np.isnan(grid_day2_lin), grid_day2_near, grid_day2_lin).ravel()

# =====================================================================
# 4. FUNGSI ELEVASI DEM, EKSTRAKSI SUNGAI UTAMA & KARTOGRAFI
# =====================================================================
def get_elevation_from_dem(point_utm32750, dem_path, label=''):
  """Ambil elevasi dari DEM sebagai rata-rata window 3x3 pixel di sekitar
  titik (bukan 1 pixel tunggal) untuk meredam noise DEM (mis. artefak
  elevasi negatif di zona pesisir). Kalau window 3x3 semuanya nodata,
  fallback ke 1 pixel pusat. PENTING: kegagalan apa pun (file rusak, CRS
  mismatch, titik di luar extent) SEKARANG dicetak sebagai peringatan
  eksplisit -- versi lama menelan semua error diam-diam dan return 0.0,
  sehingga slope_S0 bisa salah tanpa ada tanda apa pun di output."""
  if dem_path is None or not os.path.exists(dem_path):
    print(f"  [PERINGATAN] DEM tidak ditemukan untuk titik {label} -> elevasi diisi 0.0")
    return 0.0
  try:
    with rasterio.open(dem_path) as src:
      pt_gdf = gpd.GeoDataFrame(geometry=[point_utm32750], crs='EPSG:32750').to_crs(src.crs)
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

      print(f"  [PERINGATAN] Semua pixel window 3x3 di titik {label} nodata -> elevasi diisi 0.0")
      return 0.0
  except Exception as e:
    print(f"  [PERINGATAN] Gagal baca DEM di titik {label}: {e} -> elevasi diisi 0.0")
    return 0.0

def hitung_sungai_utama_standar(sungai_clipped_gdf, tolerance_m=0.0):
  """Delineasi sungai utama untuk 13 dari 14 DAS.
  CATATAN: DAS 7028 pakai hitung_sungai_utama_khusus_7028() di bawah,
  BUKAN varian dari fungsi ini. Presisi pembulatan koordinat (0 desimal),
  metode gabung geometri (snap, bukan linemerge), jumlah node kandidat (5),
  dan ketiadaan fallback garis-terpanjang semuanya sengaja beda dari versi
  khusus 7028 -- jangan digabung/disederhanakan tanpa validasi ulang hasil
  L & geom_utama untuk seluruh 14 DAS."""
  if sungai_clipped_gdf.empty: return 0.0, None
  merged_union = unary_union(sungai_clipped_gdf.geometry)
  snapped_geom = merged_union if tolerance_m == 0 else snap(merged_union, merged_union, tolerance=tolerance_m)

  lines = [snapped_geom] if isinstance(snapped_geom, LineString) else (list(snapped_geom.geoms) if isinstance(snapped_geom, MultiLineString) else [])
  if not lines: return 0.0, None

  G = nx.Graph()
  for line in lines:
    coords = list(line.coords)
    for i in range(len(coords) - 1):
      p1 = (round(coords[i][0], 0), round(coords[i][1], 0))
      p2 = (round(coords[i + 1][0], 0), round(coords[i + 1][1], 0))
      dist = Point(p1).distance(Point(p2))
      if dist > 0: G.add_edge(p1, p2, weight=dist)

  if G.number_of_nodes() == 0: return 0.0, None

  best_len, best_path = 0.0, None
  for comp in nx.connected_components(G):
    subG = G.subgraph(comp)
    if len(subG.nodes()) < 2: continue
    nodes_sorted = sorted(list(subG.nodes()), key=lambda n: n[1], reverse=True)
    top_nodes, bottom_nodes = nodes_sorted[:min(5, len(nodes_sorted))], nodes_sorted[-min(5, len(nodes_sorted)):]
    for hulu in top_nodes:
      for hilir in bottom_nodes:
        if hulu == hilir or hulu[1] <= hilir[1]: continue
        try:
          path = nx.shortest_path(subG, source=hulu, target=hilir, weight='weight')
          length = sum(subG[u][v]['weight'] for u, v in zip(path[:-1], path[1:]))
          if length > best_len: best_len, best_path = length, path
        except nx.NetworkXNoPath: continue

  return (best_len / 1000.0, LineString(best_path)) if best_path and len(best_path) > 1 else (0.0, None)

def hitung_sungai_utama_khusus_7028(sungai_clipped_gdf):
  if sungai_clipped_gdf.empty: return 0.0, None
  lines_only = []
  for geom in sungai_clipped_gdf.geometry:
    if isinstance(geom, LineString) and len(geom.coords) > 1: lines_only.append(geom)
    elif hasattr(geom, 'geoms'):
      for sub_geom in geom.geoms:
        if isinstance(sub_geom, LineString) and len(sub_geom.coords) > 1: lines_only.append(sub_geom)

  if not lines_only: return 0.0, None
  merged = linemerge(unary_union(lines_only))
  lines = [merged] if isinstance(merged, LineString) else (list(merged.geoms) if isinstance(merged, MultiLineString) else lines_only)

  G = nx.Graph()
  for line in lines:
    coords = list(line.coords)
    for i in range(len(coords) - 1):
      p1 = (round(coords[i][0], 1), round(coords[i][1], 1))
      p2 = (round(coords[i + 1][0], 1), round(coords[i + 1][1], 1))
      dist = Point(p1).distance(Point(p2))
      if dist > 0: G.add_edge(p1, p2, weight=dist)

  best_len, best_path = 0.0, None
  for comp in nx.connected_components(G):
    subG = G.subgraph(comp)
    if len(subG.nodes()) < 2: continue
    nodes_sorted = sorted(list(subG.nodes()), key=lambda n: n[1], reverse=True)
    top_nodes, bottom_nodes = nodes_sorted[:min(10, len(nodes_sorted))], nodes_sorted[-min(10, len(nodes_sorted)):]
    for hulu in top_nodes:
      for hilir in bottom_nodes:
        if hulu == hilir or hulu[1] <= hilir[1]: continue
        try:
          path = nx.shortest_path(subG, source=hulu, target=hilir, weight='weight')
          length_m = sum(subG[u][v]['weight'] for u, v in zip(path[:-1], path[1:]))
          if length_m > best_len: best_len, best_path = length_m, LineString(path)
        except nx.NetworkXNoPath: continue

  if best_path is None or best_len == 0:
    lines_sorted = sorted(lines, key=lambda l: l.length, reverse=True)
    best_path, best_len = lines_sorted[0], lines_sorted[0].length

  return best_len / 1000.0, best_path

def fmt_dms_lon(x, pos):
  pt = gpd.GeoSeries([Point(x, -967000)], crs='EPSG:3857').to_crs(epsg=4326).iloc[0]
  deg_val = abs(pt.x); d = int(deg_val); m = int((deg_val - d) * 60); s = int(round(((deg_val - d) * 60 - m) * 60))
  if s == 60: m += 1; s = 0
  if m == 60: d += 1; m = 0
  return f"{d}°{m}'{s}\"E"

def fmt_dms_lat(y, pos):
  pt = gpd.GeoSeries([Point(12820000, y)], crs='EPSG:3857').to_crs(epsg=4326).iloc[0]
  deg_val = abs(pt.y); d = int(deg_val); m = int((deg_val - d) * 60); s = int(round(((deg_val - d) * 60 - m) * 60))
  if s == 60: m += 1; s = 0
  if m == 60: d += 1; m = 0
  return f"{d}°{m}'{s}\"S"

def get_abm_order(duration_hours):
  center = int(np.ceil(duration_hours / 2.0))
  abm_order = []
  left, right = center, center + 1
  toggle = True
  for _ in range(duration_hours):
    if toggle and left >= 1:
      abm_order.append(left); left -= 1; toggle = False
    elif not toggle and right <= duration_hours:
      abm_order.append(right); right += 1; toggle = True
    elif left >= 1:
      abm_order.append(left); left -= 1
    elif right <= duration_hours:
      abm_order.append(right); right += 1
  return {r: idx + 1 for idx, r in enumerate(abm_order)}

def hitung_lc_centroid(gdf_das_utm, geom_utama_utm):
  """Lc = jarak sepanjang sungai utama dari titik pada sungai yang
  terdekat dengan centroid DAS, hingga outlet (ujung hilir sungai).
  Menggantikan pendekatan lama Lc = 0.5 * L.
  Return None jika geometri sungai tidak tersedia (caller fallback ke 0.5*L)."""
  if geom_utama_utm is None or geom_utama_utm.length == 0:
    return None
  centroid = gdf_das_utm.geometry.unary_union.centroid
  jarak_dari_hulu = geom_utama_utm.project(centroid)
  panjang_total_m = geom_utama_utm.length
  Lc_m = panjang_total_m - jarak_dari_hulu
  return Lc_m / 1000.0

def hitung_hsg_dominan(gdf_das_utm, raster_path):
  """Clip raster HYSOGs250m (pixel 1-4=A-D, 11-14=kelas gabungan) ke
  poligon DAS, reklasifikasi 11-14 -> D, lalu kembalikan kelas HSG
  paling dominan (majority) + persentase dominasi + rincian per kelas.
  Fallback ke ('B', 0.0, {}) kalau raster tidak ada / tidak overlap DAS."""
  if raster_path is None:
    return 'B', 0.0, {}

  geom_wgs84 = gdf_das_utm.to_crs(epsg=4326).geometry.unary_union

  with rasterio.open(raster_path) as src:
    geom_raster_crs = gpd.GeoSeries([geom_wgs84], crs='EPSG:4326').to_crs(src.crs).iloc[0]
    try:
      out_image, _ = rasterio.mask.mask(src, [geom_raster_crs], crop=True, nodata=src.nodata)
    except ValueError:
      return 'B', 0.0, {}  # poligon DAS di luar extent raster
    nodata_val = src.nodata

  pixels = out_image[0]
  valid_mask = (pixels != nodata_val) if nodata_val is not None else np.ones_like(pixels, dtype=bool)
  valid_mask &= (pixels > 0)
  pixel_valid = pixels[valid_mask]

  if pixel_valid.size == 0:
    return 'B', 0.0, {}

  kelas_hsg = np.array([MAP_HSG_PIXEL.get(int(v)) for v in pixel_valid])
  kelas_hsg = kelas_hsg[kelas_hsg != None]

  if kelas_hsg.size == 0:
    return 'B', 0.0, {}

  unique, counts = np.unique(kelas_hsg, return_counts=True)
  total = kelas_hsg.size
  rincian_persen = {k: round(v / total * 100, 1) for k, v in zip(unique, counts)}
  kelas_dominan = max(rincian_persen, key=rincian_persen.get)
  persen_dominan = rincian_persen[kelas_dominan]

  return kelas_dominan, persen_dominan, rincian_persen

def hitung_cn_dari_hsg(persen_lahan_dict, hsg_kelas):
  """Hitung CN(II) komposit untuk satu DAS berdasarkan persentase area
  per kelas tata guna lahan (dict dari dict_persen_lahan) dan kelas HSG
  dominan (A/B/C/D), memakai TABEL_CN_HSG."""
  cn_air = (persen_lahan_dict['Air_%'] / 100.0) * TABEL_CN_HSG['Air'][hsg_kelas]
  cn_bangunan = (persen_lahan_dict['Bangunan_%'] / 100.0) * TABEL_CN_HSG['Bangunan'][hsg_kelas]
  cn_hutan = (persen_lahan_dict['Hutan_%'] / 100.0) * TABEL_CN_HSG['Hutan'][hsg_kelas]
  cn_terbuka = (persen_lahan_dict['Terbuka_%'] / 100.0) * TABEL_CN_HSG['Terbuka'][hsg_kelas]
  cn_ii = cn_air + cn_bangunan + cn_hutan + cn_terbuka
  s_ii_mm = (25400 / cn_ii) - 254 if cn_ii > 0 else 50.0
  return cn_ii, s_ii_mm

row_to_rank_19 = get_abm_order(19)
row_to_rank_18 = get_abm_order(18)

# =====================================================================
# 5. ITERASI PEMROSESAN 14 DAS
# =====================================================================
print('--- TAHAP 3: EKSEKUSI PEMODELAN KONTINU 37 JAM ---')

unique_gridcodes = gdf_das_all['gridcode'].unique()

with pd.ExcelWriter(EXCEL_OUTPUT_PATH, engine='openpyxl') as writer:
  df_luas.to_excel(writer, sheet_name='Rekap_Tata_Guna_Lahan_CN', index=False)
  
  for gc in unique_gridcodes:
    gc_str = str(int(gc)) if isinstance(gc, float) and gc.is_integer() else str(gc)
    gc_val = int(gc) if isinstance(gc, float) and gc.is_integer() else gc
    print(f'-> Memproses DAS Gridcode: {gc_str}')
    
    gdf_single_das = gdf_das_all[gdf_das_all['gridcode'] == gc]
    gdf_das_wgs84 = gdf_single_das.to_crs(epsg=4326)
    gdf_das_utm = gdf_single_das.to_crs(epsg=32750)

    A = gdf_das_utm.geometry.area.sum() / 1_000_000.0

    geom_das_union = gdf_das_wgs84.geometry.unary_union
    mask_day1 = gdf_grid_res_day1.geometry.within(geom_das_union)
    mask_day2 = gdf_grid_res_day2.geometry.within(geom_das_union)

    if mask_day1.sum() > 0:
      R_day1 = round(gdf_grid_res_day1[mask_day1]['ch'].mean(), 2)
      R_day2 = round(gdf_grid_res_day2[mask_day2]['ch'].mean(), 2)
    else:
      centroid_das = gdf_das_wgs84.geometry.centroid.iloc[0]
      _, idx_terdekat = tree_hujan.query([centroid_das.x, centroid_das.y])
      R_day1 = round(val_day1[idx_terdekat], 2)
      R_day2 = round(val_day2[idx_terdekat], 2)

    das_geom_utm_union = gdf_das_utm.geometry.unary_union
    sungai_filtered = gdf_sungai_all[gdf_sungai_all.geometry.intersects(das_geom_utm_union)].copy()
    if not sungai_filtered.empty:
      sungai_filtered['geometry'] = sungai_filtered.geometry.intersection(das_geom_utm_union)
      sungai_das = sungai_filtered[~sungai_filtered.geometry.is_empty]
    else:
      sungai_das = gpd.GeoDataFrame(columns=gdf_sungai_all.columns, crs='EPSG:32750')

    if "7028" in gc_str:
      L, geom_utama = hitung_sungai_utama_khusus_7028(sungai_das)
    else:
      L, geom_utama = hitung_sungai_utama_standar(sungai_das, tolerance_m=0.0)

    Lc = hitung_lc_centroid(gdf_das_utm, geom_utama)
    if Lc is None or Lc <= 0:
      Lc = 0.5 * L  # fallback kalau sungai utama tidak terdeteksi

    hsg_dominan, hsg_persen, hsg_rincian = hitung_hsg_dominan(gdf_das_utm, hysogs_raster_path)
    label_hsg = f"Dominan HSG {hsg_dominan} ({hsg_persen:.1f}%)" if hysogs_raster_path else "Group B (asumsi, raster HYSOGs tidak ditemukan)"

    if geom_utama is not None:
      coords_main = list(geom_utama.coords)
      pt_hulu = Point(coords_main[0])
      pt_hilir = Point(coords_main[-1])
      z_hulu = get_elevation_from_dem(pt_hulu, dem_raster_path, label=f'hulu DAS {gc_str}')
      z_hilir = get_elevation_from_dem(pt_hilir, dem_raster_path, label=f'hilir DAS {gc_str}')
      panjang_m = L * 1000.0
      slope_S0 = ((z_hulu - z_hilir) / panjang_m) if panjang_m > 0 else 0.0
      if slope_S0 < 0: slope_S0 = abs(slope_S0)
    else:
      z_hulu, z_hilir, slope_S0 = 0.0, 0.0, 0.0

    Ct = CT_SNYDER
    Tc_kirpich = hitung_tc_kirpich(L, slope_S0)
    if Tc_kirpich is not None and Tc_kirpich > 0:
      Tr = FRAKSI_TR_DARI_TC * Tc_kirpich
    else:
      Tr = 1.0  # fallback kalau L atau slope_S0 tidak valid (mis. sungai utama gagal terdeteksi)
    TL = Ct * ((L * Lc) ** 0.3) if (L > 0 and Lc > 0) else 1.0
    TP = TL + 0.5 * Tr
    TB = 5.0 * TP
    QP = (0.2083 * A) / TP if TP > 0 else 0.0

    t1 = np.arange(1, 20)
    Ir1 = (R_day1 / 24.0) * ((24.0 / t1) ** (2.0 / 3.0))
    Ir_Td1 = t1 * Ir1
    delta_p1 = np.zeros_like(Ir_Td1); delta_p1[0] = Ir_Td1[0]; delta_p1[1:] = np.diff(Ir_Td1)
    pt1_sorted = np.sort(delta_p1 / np.sum(delta_p1))[::-1]
    pt2_day1 = np.array([pt1_sorted[row_to_rank_19[r] - 1] for r in range(1, 20)])
    p_hourly1 = R_day1 * pt2_day1

    t2 = np.arange(1, 19)
    Ir2 = (R_day2 / 24.0) * ((24.0 / t2) ** (2.0 / 3.0))
    Ir_Td2 = t2 * Ir2
    delta_p2 = np.zeros_like(Ir_Td2); delta_p2[0] = Ir_Td2[0]; delta_p2[1:] = np.diff(Ir_Td2)
    pt2_sorted = np.sort(delta_p2 / np.sum(delta_p2))[::-1]
    pt2_day2 = np.array([pt2_sorted[row_to_rank_18[r] - 1] for r in range(1, 19)])
    p_hourly2 = R_day2 * pt2_day2

    t_37_rel = np.concatenate([t1, t2])
    Ir_37 = np.concatenate([Ir1, Ir2])
    Ir_Td_37 = np.concatenate([Ir_Td1, Ir_Td2])
    delta_p_37 = np.concatenate([delta_p1, delta_p2])
    pt_37 = np.concatenate([pt2_day1, pt2_day2])
    p_37 = np.concatenate([p_hourly1, p_hourly2])

    sigma_P37 = np.cumsum(p_37)

    CN_II_base, _ = hitung_cn_dari_hsg(dict_persen_lahan[gc_val], hsg_dominan)

    R5_das = hitung_r5_antesenden(geom_das_union, df_hujan_grid)
    kelas_amc = klasifikasi_amc(R5_das)
    if R5_das is None:
      print(f"  [PERINGATAN] Data hujan 5 hari antesenden tidak ditemukan untuk DAS {gc_str} "
            f"-> AMC diasumsikan II (tidak ada penyesuaian CN)")

    # Urutan wajib: (1) CN dasar dari tabel HSG (basis lambda=0.2) ->
    # (2) sesuaikan AMC (formula NEH-4 juga berbasis lambda=0.2) ->
    # (3) BARU konversi seluruh hasil ke basis lambda=0.05. Membalik urutan
    # ini (konversi dulu baru AMC) akan salah karena formula AMC di atas
    # tidak berlaku untuk CN berbasis lambda=0.05.
    CN_amc = adjust_cn_amc(CN_II_base, kelas_amc)
    CN_final = cn_ke_lambda_005(CN_amc) if LAMBDA_IA == 0.05 else CN_amc
    S_base = (25400 / CN_final) - 254 if CN_final > 0 else 50.0
    Ia_base = LAMBDA_IA * S_base

    sigma_Peff37 = np.where(sigma_P37 > Ia_base, ((sigma_P37 - Ia_base)**2) / (sigma_P37 - Ia_base + S_base), 0.0)
    peff_37 = np.zeros_like(sigma_Peff37)
    peff_37[0] = sigma_Peff37[0]
    peff_37[1:] = np.diff(sigma_Peff37)
    infil_37 = p_37 - peff_37

    jam_day1 = [f'{h:02d}:00' for h in range(5, 24)]
    jam_day2 = [f'{h:02d}:00' for h in range(0, 18)]
    list_jam_kejadian = jam_day1 + jam_day2

    df_hujan_ext = pd.DataFrame({
        't Kontinu (jam)': np.arange(1, 38),
        'Tanggal Event': ['9 Sept 2025']*19 + ['10 Sept 2025']*18,
        'Jam Kejadian': list_jam_kejadian,
        't (relatif)': t_37_rel,
        'Ir (mm/jam)': np.round(Ir_37, 2),
        'Ir * Td': np.round(Ir_Td_37, 2),
        'Δp (mm)': np.round(delta_p_37, 2),
        'pt (ABM)': [f'{v*100:.2f}%' for v in pt_37],
        'P (mm)': np.round(p_37, 2),
        'ΣP (mm)': np.round(sigma_P37, 2),
        'S (mm)': np.round(S_base, 2),
        'Ia (mm)': np.round(Ia_base, 2),
        'Σ Peff (mm)': np.round(sigma_Peff37, 2),
        'Peff (mm)': np.round(peff_37, 2),
        'Infil (mm)': np.round(infil_37, 2)
    })

    t_tp_std = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0])
    q_qp_std = np.array([0.0, 0.03, 0.1, 0.19, 0.31, 0.47, 0.66, 0.82, 0.93, 0.99, 1.0, 0.99, 0.93, 0.86, 0.78, 0.68, 0.56, 0.46, 0.39, 0.33, 0.28, 0.2435, 0.207, 0.177, 0.147, 0.127, 0.107, 0.092, 0.077, 0.066, 0.055, 0.0475, 0.04, 0.0345, 0.029, 0.025, 0.021, 0.018, 0.015, 0.013, 0.011, 0.0098, 0.0086, 0.0074, 0.0062, 0.005, 0.004, 0.003, 0.002, 0.001, 0.0])

    t_hidrograf = t_tp_std * TP
    q_hidrograf = q_qp_std * QP
    max_jam_interpolasi = int(np.ceil(TB)) + 37
    jam_1jam = np.arange(0, max_jam_interpolasi + 1)
    q_1jam = np.interp(jam_1jam, t_hidrograf, q_hidrograf, right=0.0)

    num_rows = len(q_1jam)
    num_rain = 37
    matrix_hidrograf = np.zeros((num_rows, num_rain))

    for j in range(num_rain):
      p_val = peff_37[j]
      if j < num_rows:
        length_to_copy = min(len(q_1jam), num_rows - j)
        matrix_hidrograf[j : j + length_to_copy, j] = q_1jam[:length_to_copy] * p_val

    debit_banjir_total = np.sum(matrix_hidrograf, axis=1)
    dict_konvolusi = {'t (jam)': jam_1jam}
    for j in range(num_rain):
      dict_konvolusi[f'Jam-{j+1} ({peff_37[j]:.2f}mm)'] = np.round(matrix_hidrograf[:, j], 3)
    dict_konvolusi['Debit Total (m³/s)'] = np.round(debit_banjir_total, 2)
    df_konvolusi_final = pd.DataFrame(dict_konvolusi)

    fig, ax = plt.subplots(figsize=(6, 8), dpi=200)

    gdf_das_web = gdf_das_all.to_crs(epsg=3857)
    gdf_fokus_web = gdf_das_web[gdf_das_web['gridcode'] == gc]
    gdf_around_web = gdf_das_web[gdf_das_web['gridcode'] != gc]

    gdf_around_web.plot(ax=ax, facecolor='#D1C4E9', edgecolor='black', linewidth=0.6, alpha=0.5, label='CA')
    gdf_fokus_web.plot(ax=ax, facecolor='#FFD700', edgecolor='black', linewidth=1.2, alpha=0.85, label=f'CA - {gc_str}')

    if not sungai_das.empty:
      sungai_web = sungai_das.to_crs(epsg=3857)
      sungai_web.plot(ax=ax, color='#1E88E5', linewidth=0.8, alpha=0.7, label='River')

    if geom_utama is not None:
      geom_utama_web = gpd.GeoSeries([geom_utama], crs='EPSG:32750').to_crs(epsg=3857).iloc[0]
      gpd.GeoSeries([geom_utama_web], crs='EPSG:3857').plot(ax=ax, color='red', linewidth=2.2, label='Main River')
      coords_w = list(geom_utama_web.coords)
      ax.plot(coords_w[0][0], coords_w[0][1], 'go', markersize=5)
      ax.plot(coords_w[-1][0], coords_w[-1][1], 'mo', markersize=5)

    try:
      cx.add_basemap(ax, source=cx.providers.Esri.WorldImagery, zoom=12, attribution="")
    except Exception:
      pass

    bounds_fokus = gdf_fokus_web.total_bounds
    margin_x = (bounds_fokus[2] - bounds_fokus[0]) * 0.45
    margin_y = (bounds_fokus[3] - bounds_fokus[1]) * 0.45
    ax.set_xlim(bounds_fokus[0] - margin_x, bounds_fokus[2] + margin_x)
    ax.set_ylim(bounds_fokus[1] - margin_y, bounds_fokus[3] + margin_y)

    ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=2))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=4))

    ax.xaxis.set_major_formatter(FuncFormatter(fmt_dms_lon))
    ax.yaxis.set_major_formatter(FuncFormatter(fmt_dms_lat))
    
    ax.tick_params(axis='both', which='major', labelsize=6, labeltop=True, labelbottom=True, labelleft=True, labelright=True)
    ax.grid(True, linestyle=':', color='gray', alpha=0.5)

    ax.legend(loc='lower right', frameon=True, facecolor='white', framealpha=0.9, fontsize=7)

    esri_attr_text = "Tiles (C) Esri — Source: Esri, i-cubed, USDA, USGS, AEX, GeoEye, Getmapping, Aerogrid, IGN, IGP, UPR-EGP, and the GIS User Community"
    ax.text(0.02, 0.02, esri_attr_text, transform=ax.transAxes, fontsize=3.5, fontweight='bold',
            color='white', bbox=dict(boxstyle="square,pad=0.2", facecolor="black", alpha=0.5, edgecolor="none"),
            wrap=True, verticalalignment='bottom')

    img_buffer = io.BytesIO()
    plt.savefig(img_buffer, format='png', bbox_inches='tight', dpi=150)
    img_buffer.seek(0)
    plt.close(fig)

    sheet_name = f"DAS_{gc_str}"
    df_rekap = pd.DataFrame({
        'Parameter Hidrologi & Geometri': [
            'Gridcode DAS', 'Hydrologic Soil Group (HYSOGs250m)',
            'Curve Number (CN II, basis lambda=0.2)',
            'Hujan 5 Hari Antesenden (R5) [mm]', 'Kelas AMC',
            f'Curve Number Final (basis lambda={LAMBDA_IA})',
            'Luas DAS (A) [km²]',
            'Panjang Sungai Utama (L) [km]', 'Jarak Centroid-Outlet (Lc) [km]',
            'Elevasi Hulu Sungai (Z_hulu) [m]', 'Elevasi Hilir Sungai (Z_hilir) [m]',
            'Kemiringan Saluran (S0) [m/m]', 'Kemiringan Saluran (S0) [%]',
            'Hujan Harian 9 Sept (R1) [mm]', 'Hujan Harian 10 Sept (R2) [mm]',
            'Total Durasi Kejadian [jam]',
            'Time of Concentration Kirpich (Tc) [jam]', 'Durasi Hujan Satuan (Tr) [jam]',
            'Time Lag (TL) [jam]', 'Time Peak (TP) [jam]', 
            'Time Base (TB) [jam]', 'Debit Puncak Hidrograf Satuan (QP) [m³/s/mm]',
            'DEBIT PUNCAK BANJIR TOTAL (KONTINU) [m³/s]'
        ],
        'Nilai': [
            gc_str, label_hsg, round(CN_II_base, 2),
            'n/a' if R5_das is None else round(R5_das, 2), kelas_amc,
            round(CN_final, 2),
            round(A, 2), round(L, 2), round(Lc, 2),
            round(z_hulu, 2), round(z_hilir, 2), round(slope_S0, 5), round(slope_S0 * 100, 3),
            round(R_day1, 2), round(R_day2, 2), 37,
            'n/a' if Tc_kirpich is None else round(Tc_kirpich, 3), round(Tr, 3),
            round(TL, 3), round(TP, 3), round(TB, 3), round(QP, 3),
            round(np.max(debit_banjir_total), 2)
        ]
    })

    df_rekap.to_excel(writer, sheet_name=sheet_name, index=False, startrow=0)
    
    start_row_hujan = len(df_rekap) + 3
    worksheet = writer.sheets[sheet_name]
    worksheet.cell(row=start_row_hujan, column=1, value="TABEL INTENSITAS HUJAN 37 JAM (MONONOBE, ABM & EKSTENSI SCS-CN)")
    df_hujan_ext.to_excel(writer, sheet_name=sheet_name, index=False, startrow=start_row_hujan)

    start_row_konv = start_row_hujan + len(df_hujan_ext) + 3
    worksheet.cell(row=start_row_konv, column=1, value="TABEL SUPERPOSISI KONVOLUSI HIDROGRAF DEBIT BANJIR (37 JAM KONTINU)")
    df_konvolusi_final.to_excel(writer, sheet_name=sheet_name, index=False, startrow=start_row_konv)

    img_temp_path = f"/tmp/plot_das_{gc_str}.png"
    with open(img_temp_path, "wb") as f: f.write(img_buffer.getbuffer())
    
    xl_img = OpenpyxlImage(img_temp_path)
    xl_img.width, xl_img.height = 340, 440
    worksheet.add_image(xl_img, "AU2")

print(f'\n=====================================================================')
print(f'PEMROSESAN SELESAI DENGAN FILE TATA GUNA LAHAN BARU:\n{EXCEL_OUTPUT_PATH}')
print(f'=====================================================================')