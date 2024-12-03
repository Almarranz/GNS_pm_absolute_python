        
# GNS_pm_absolute_offsets using just python scripts (and no astroalign)

## Obtaining Absolute Proper Motions from GNS1 and GNS2 Using GAIA as the Reference Frame

### GRF (ICRS): Offsets (GNS xy coordinates and Deprojected Gaia Coordinates)

0. **`gns1_B_lists.py`**  
   Combines the H and Ks lists from GNS1.

1. **`gaia_refstars_offsets_gns_pixels.py`**  
   Generates Gaia offset coordinates, shifting them to the corresponding epoch. This script uses `astroalign` to align the *GNS xy coordinates* with the Gaia offsets.

2. **`IDL_PY_gaia_offset_gns_xy_align.pro`**  
   Aligns GNS with Gaia using a polynomial of degree 1 or 2. This script uses the **Python** script `compare_lists` and then calls the script **`PY_GaiaRf_align_gns_offsets.py`**.

   - 2.1 **(A) `IDL_BSgaiaRF_offsets_align_gns.pro`**  
     Generates lists using bootstrapping to estimate alignment uncertainties.

   - 2.1 **(B) `IDL_JKgaiaRF_offsets_align_gns.pro`**  
     Generates lists using jackknifing to estimate alignment uncertainties.

   - 2.2 **(A) `bs_uncert_offsets.py`**  
     Provides the RA and Dec bootstrapping uncertainties for each star.

   - 2.2 **(B) `jk_uncert_offsets.py`**  
     Provides the RA and Dec jackknife uncertainties for each star.

3. **(Old) `IDL_GaiaRf_align_gns_offsets.pro`**  
   Finds matching stars and computes the offset for proper motions.

3. **`PY_GaiaRf_align_gns_offsets.py`** (replaces `IDL_GaiaRf_align_gns_offsets.pro`, faster execution)  
   Calculates proper motions from GNSI and GNSII.

   - 3.1 **`gaia_stars_offsets_comp.py`**  
     Computes residuals with the Gaia catalog. Identifies and excludes 3-sigma outliers, then repeats the alignment (from step 2) without them.

4. **`pm_analysis_offsets_radec.py`**  
   Identifies and visualizes clusters using DBSCAN.

#### Quality Check

5. **`Hosek_cluster_com.py`**  
   Calculates proper motion and position residuals between Hosek and GNS.

6. **`gaia_hosek_gaia_gns.py`**  
   Computes residuals between Gaia stars in Hosek's catalog and those in the GNS catalog.

   - 6.1 **`gaia_libralato_gaia_gns.py`**  
     Computes residuals between Gaia stars in Hosek's catalog and GNS's Gaia stars.

7. **`Libralato_comp.py`**  
   Compares GNS proper motions with Libralato's measurements.

8. **`Hosek_dbscan.pro`**  
   Runs DBSCAN on Hosek's data.

9. **`NSD_icrs_disk.py`**  
   Uses Hosek's data to visualize the proper motion distribution of the disk (Galactic and relative), providing insights into expected values.

> **Note**: `NSD_icrs_disk.py` is located in `/Users/amartinez/Desktop/PhD/HAWK/GNS_2/GNS_pm`.
