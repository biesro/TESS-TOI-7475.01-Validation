"""
02_centroid_FINAL.py
--------------------
DEFINITIVE Centroid Test for TOI 7475.01.
Objective: Determine if the star's position shifts during the eclipse.

INTERPRETATION:
- If the trend lines are FLAT = SUCCESS (Likely a Planet).
- If the trend lines look like MOUNTAINS/VALLEYS = FALSE POSITIVE (Background Binary).
"""

import lightkurve as lk
import matplotlib.pyplot as plt
import numpy as np
import warnings

# Ignore technical warnings for a cleaner output
warnings.simplefilter('ignore')

# --- 1. UPDATED DATA (From your Elite Analysis) ---
TIC_ID = "TIC 376866659"
KNOWN_PERIOD = 3.253773  # <--- Precise updated value
KNOWN_T0 = 3775.5819     # <--- Precise updated value
SECTOR = 91             

print(f"🚀 STARTING CENTROID TEST: {TIC_ID}")
print(f"   Parameters: P={KNOWN_PERIOD:.4f} d, T0={KNOWN_T0:.4f}")
print("-" * 50)

try:
    # 1. DOWNLOAD THE PIXEL 'VIDEO' (TPF)
    print("Step 1: Downloading Target Pixel File (TPF)...")
    search = lk.search_targetpixelfile(TIC_ID, author="SPOC", sector=SECTOR)
    
    if len(search) == 0:
        print("⚠️ SPOC not found. Trying TESS-SPOC...")
        search = lk.search_targetpixelfile(TIC_ID, author="TESS-SPOC", sector=SECTOR)
        
    if len(search) == 0:
        # Last attempt with QLP if nothing else exists
        print("⚠️ TESS-SPOC not found. Trying QLP...")
        search = lk.search_targetpixelfile(TIC_ID, author="QLP", sector=SECTOR)

    if len(search) == 0:
        raise ValueError("❌ NO PIXEL DATA (TPF) FOUND. Cannot perform centroid test.")

    tpf = search.download()
    print("✅ TPF downloaded successfully.")

    # 2. CALCULATE THE CENTER OF LIGHT
    print("Step 2: Calculating the center of light for each frame...")
    # We use the 'moments' method, which is the standard for this analysis
    centroids = tpf.estimate_centroids(method='moments')
    
    # Separate X (Columns) and Y (Rows) coordinates
    x_vals = centroids[0].value
    y_vals = centroids[1].value
    time = tpf.time.value
    
    # 3. DATA CLEANING (Crucial to avoid plotting errors)
    # Create a mask to remove empty points (NaNs)
    mask = ~np.isnan(x_vals) & ~np.isnan(y_vals)
    
    time_clean = time[mask]
    x_clean = x_vals[mask]
    y_clean = y_vals[mask]
    
    print(f"✅ Data cleaned: {len(time_clean)} valid points.")

    # 4. FOLDING
    # The trick: We create 'fake' lightcurves where the flux is actually the X or Y position
    print("Step 3: Folding data with the planet's period...")
    
    lc_x = lk.LightCurve(time=time_clean, flux=x_clean)
    lc_y = lk.LightCurve(time=time_clean, flux=y_clean)
    
    folded_x = lc_x.fold(period=KNOWN_PERIOD, epoch_time=KNOWN_T0)
    folded_y = lc_y.fold(period=KNOWN_PERIOD, epoch_time=KNOWN_T0)

    # 5. VISUALIZATION
    print("Step 4: Generating Plots...")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # -- PLOT X (Columns) --
    # Gray points = Raw motion noise
    folded_x.scatter(ax=ax1, s=3, alpha=0.3, c='gray', label='Raw Motion')
    # Red line = Average movement (BINNING)
    # Using a larger bin (0.005) to see the clear trend
    folded_x.bin(time_bin_size=0.005).plot(ax=ax1, c='red', lw=3, label='Average (Trend)')
    
    ax1.set_title("X Motion (Columns)", fontsize=14, fontweight='bold')
    ax1.set_ylabel("Shift (Pixels)")
    ax1.set_xlabel("Phase")
    ax1.set_xlim(-0.15, 0.15) # Zoom into the transit
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # -- PLOT Y (Rows) --
    folded_y.scatter(ax=ax2, s=3, alpha=0.3, c='gray')
    folded_y.bin(time_bin_size=0.015).plot(ax=ax2, c='blue', lw=3, label='Average (Trend)')
    
    ax2.set_title("Y Motion (Rows)", fontsize=14, fontweight='bold')
    ax2.set_xlabel("Phase")
    ax2.set_xlim(-0.15, 0.15)
    ax2.grid(True, alpha=0.3)
    
    plt.suptitle(f"CENTROID TEST: TOI 7475.01 (Sector {SECTOR})", fontsize=16)
    plt.tight_layout()
    
    print("🖥️  OPENING WINDOW... CHECK YOUR TASKBAR!")
    plt.show()

except Exception as e:
    print("\n" + "!"*50)
    print(f"❌ FATAL ERROR: {e}")
    print("!"*50)