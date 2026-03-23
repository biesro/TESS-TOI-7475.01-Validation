"""
01_detection_BLS.py
=============================================================================
HYBRID ANALYSIS PIPELINE: NATURAL DATA + ELITE PHYSICS
=============================================================================

DESCRIPTION:
This script merges two analytical philosophies to validate Exoplanet TOI 7475.01.

1. DATA PROCESSING (The "Natural" Approach):
   - We apply minimal cleaning (removing NaNs/Outliers) but STRICTLY AVOID 
     flattening or aggressive detrending.
   - Goal: To preserve the natural stellar variability and the true geometric 
     shape of the transit without distortion.

2. DETECTION ENGINE (The "Elite" Approach):
   - Multi-Duration BLS: Searches for signals across various transit durations 
     to find the best fit.
   - Robust SNR: Mathematically calculates Signal-to-Noise Ratio to gauge 
     detection confidence.
   - Odd/Even Veto: Splits transits to check for binary star false positives.

TARGET:
- ID:     TOI 7475.01 / TIC 376866659
- Output: Professional plots and CSV data export.
"""

import lightkurve as lk
import numpy as np
import matplotlib.pyplot as plt
import gc
import warnings

warnings.simplefilter('ignore')

# --- CONFIGURATION ---
TIC_ID = "TIC 376866659"
TOI_ID = "TOI 7475.01"

SEARCH_PERIOD_MIN = 3.15    
SEARCH_PERIOD_MAX = 3.35
DURATIONS_TO_TEST = np.linspace(0.10, 0.25, 10)

print(f"🔬 FUSION ANALYSIS FOR {TOI_ID}")
print("---------------------------------------------------------")

# --- DATA PIPELINE (BASIC CLEANING WITHOUT FLATTEN) ---
print("📡 Connecting to MAST servers...")

try:
    search = lk.search_lightcurve(TIC_ID, author="SPOC", exptime=120)
    if len(search) == 0:
        search = lk.search_lightcurve(TIC_ID, author="TESS-SPOC")
    if len(search) == 0:
        search = lk.search_lightcurve(TIC_ID, author="QLP")

    if len(search) == 0:
        raise ValueError("❌ No data found.")

    print(f"✅ Data found! Sectors: {len(search)}")
    
    lc_collection = []

    for i, lc_item in enumerate(search):
        try:
            print(f"   -> Processing Sector {lc_item.mission[0]}...")
            temp_lc = lc_item.download()
            
            # ✅ BASIC CLEANING (LIKE VISUALIZER)
            temp_lc = temp_lc.remove_nans().normalize()
            temp_lc = temp_lc.remove_outliers(sigma_upper=4, sigma_lower=15)
            
            # ❌ NO FLATTEN - We leave data natural!
            
            lc_collection.append(temp_lc)
            del temp_lc
            gc.collect()
            
        except Exception as e:
            print(f"   ⚠️ Error reading sector: {e}")
            continue

    if not lc_collection:
        raise ValueError("No clean data extracted.")

    lc_combined = lk.LightCurveCollection(lc_collection).stitch()
    lc_combined.flux = lc_combined.flux.astype(np.float32)
    
    print(f"✅ Pipeline complete. Clean points: {len(lc_combined)}")

except Exception as e:
    print(f"❌ CRITICAL ERROR: {e}")
    exit()

# --- BLS SEARCH (MULTI-DURATION) ---
print("\n🔍 Running High-Res BLS (Multi-Duration)...")

period_grid = np.linspace(SEARCH_PERIOD_MIN, SEARCH_PERIOD_MAX, 100000)

bls = lc_combined.to_periodogram(
    method='bls', 
    period=period_grid, 
    duration=DURATIONS_TO_TEST
)

best_period = bls.period_at_max_power.value
best_t0 = bls.transit_time_at_max_power.value
best_depth = bls.depth_at_max_power.value
best_duration = bls.duration_at_max_power.value

print(f"🎯 MATCH FOUND:")
print(f"   Period:   {best_period:.6f} days")
print(f"   T0:       {best_t0:.4f} BTJD")
print(f"   Duration: {best_duration*24:.2f} hours")
print(f"   Depth:    {best_depth:.4f} ({best_depth*1e6:.0f} ppm)")

# --- ROBUST SNR CALCULATION ---
print("\n🧮 Calculating Signal-to-Noise Ratio (SNR)...")

lc_folded = lc_combined.fold(period=best_period, epoch_time=best_t0)

phase_mask_transit = (np.abs(lc_folded.phase.value) < (best_duration * 0.55))
phase_mask_out = (np.abs(lc_folded.phase.value) > (best_duration * 2.0))

flux_in = lc_folded.flux[phase_mask_transit]
flux_out = lc_folded.flux[phase_mask_out]

if len(flux_in) > 5 and len(flux_out) > 5:
    median_flux_out = np.nanmedian(flux_out)
    median_flux_in = np.nanmedian(flux_in)
    
    depth_calc = median_flux_out - median_flux_in
    noise = np.nanstd(flux_out)
    snr_final = (depth_calc / noise) * np.sqrt(len(flux_in))
    
    depth_ppm = depth_calc * 1e6
    noise_ppm = noise * 1e6
else:
    print("⚠️ Warning: Not enough points to calculate reliable SNR.")
    snr_final = 0.0
    depth_ppm = 0.0

# --- PROFESSIONAL PLOTS ---
fig = plt.figure(figsize=(14, 8))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.2])

# A: BLS Periodogram
ax1 = fig.add_subplot(gs[0, :])
ax1.plot(bls.period, bls.power, color='black', lw=1)
ax1.set_title(f"BLS Periodogram - Peak at {best_period:.4f} d", fontsize=13, fontweight='bold')
ax1.set_xlabel("Period [days]", fontsize=11)
ax1.set_ylabel("BLS Power", fontsize=11)
ax1.axvline(best_period, color='red', alpha=0.5, lw=6, label=f'Detected: {best_period:.4f} d')
ax1.legend()
ax1.grid(True, alpha=0.3)

# B: Folded Transit (NATURAL - LIKE VISUALIZER)
ax2 = fig.add_subplot(gs[1, 0])
# Gray points
lc_folded.scatter(ax=ax2, color='gray', alpha=0.3, s=2, label='Raw Data')
# 10-minute binning (0.007 days)
bin_size = 10 / (24*60)
lc_binned = lc_folded.bin(time_bin_size=bin_size)
lc_binned.scatter(ax=ax2, color='red', s=25, marker='o', label='Binned (10min)', zorder=5)

ax2.set_title(f"Folded Transit - Natural View (T0={best_t0:.2f})", fontsize=13, fontweight='bold')
ax2.set_xlabel("Phase [days]", fontsize=11)
ax2.set_ylabel("Normalized Flux", fontsize=11)
ax2.set_xlim(-0.12, 0.12)
ax2.legend(loc='lower right')
ax2.grid(True, alpha=0.3)

# C: Odd/Even Check
ax3 = fig.add_subplot(gs[1, 1])
transit_number = np.round((lc_combined.time.value - best_t0) / best_period)
is_even = (transit_number % 2 == 0)

lc_even = lc_combined[is_even].fold(period=best_period, epoch_time=best_t0)
lc_odd = lc_combined[~is_even].fold(period=best_period, epoch_time=best_t0)

if len(lc_even) > 10:
    lc_even.bin(time_bin_size=bin_size).plot(ax=ax3, color='blue', lw=2.5, label='Even Transits')
if len(lc_odd) > 10:
    lc_odd.bin(time_bin_size=bin_size).plot(ax=ax3, color='orange', lw=2.5, linestyle='--', label='Odd Transits')

ax3.set_title("Odd vs Even Transit Check", fontsize=13, fontweight='bold')
ax3.set_xlabel("Phase [days]", fontsize=11)
ax3.set_ylabel("Normalized Flux", fontsize=11)
ax3.set_xlim(-0.1, 0.1)
ax3.legend()
ax3.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# --- FINAL REPORT ---
print("\n" + "="*60)
print(f"📊 FUSION REPORT FOR {TOI_ID}")
print("="*60)
print(f"Period Detected:    {best_period:.6f} days")
print(f"Transit Depth:      {depth_ppm:.0f} ppm")
print(f"Duration:           {best_duration*24:.2f} hours")
print(f"Epoch (T0):         {best_t0:.4f} BTJD")
print("-" * 60)
print(f"SNR (Signal/Noise): {snr_final:.2f}")
print(f"Noise Level:        {noise_ppm:.0f} ppm")
print("-" * 60)

if snr_final > 10:
    print("✅ SIGNAL STATUS: EXCELLENT. High confidence detection.")
elif snr_final > 7:
    print("✅ SIGNAL STATUS: GOOD. Likely real transit.")
elif snr_final > 5:
    print("⚠️ SIGNAL STATUS: MODERATE. Needs validation.")
else:
    print("⚠️ SIGNAL STATUS: WEAK. Requires caution.")
    
print("="*60)
print("\n🎯 NOTES:")
print("   - No detrending applied (natural stellar variations preserved)")
print("   - Multi-duration BLS search ensures optimal fit")
print("   - Odd/Even check helps identify potential false positives")
print("="*60)
# -----------------------------------------------------------
# --- DATA EXPORT (CSV) ---
# Add this to the end of your script
# -----------------------------------------------------------
import pandas as pd

print("\n💾 EXPORTING PLOT DATA...")

# 1. RED Points (Binned / 10 min average)
# These are best for analyzing transit shape
df_binned = pd.DataFrame({
    'time_days_from_transit': lc_binned.time.value, # Time from center (0.0)
    'flux': lc_binned.flux.value,
    'flux_err': lc_binned.flux_err.value
})
# Sort by time so it appears ordered in CSV
df_binned = df_binned.sort_values(by='time_days_from_transit')

# Clean filename
filename_red = f"{TOI_ID.replace(' ', '')}_Folded_Binned_RedPoints.csv"
df_binned.to_csv(filename_red, index=False, header=False)
print(f"✅ Saved (Red Points): {filename_red}")


# 2. GRAY Points (Raw Folded)
# This is the background point cloud
df_raw = pd.DataFrame({
    'time_days_from_transit': lc_folded.time.value,
    'flux': lc_folded.flux.value,
    'flux_err': lc_folded.flux_err.value
})
df_raw = df_raw.sort_values(by='time_days_from_transit')

# Raw filename
filename_gray = f"{TOI_ID.replace(' ', '')}_Folded_Raw_GrayPoints.csv"
df_raw.to_csv(filename_gray, index=False, header=False)
print(f"✅ Saved (Gray Points):   {filename_gray}")
