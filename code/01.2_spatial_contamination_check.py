"""
01.2_spatial_contamination_check.py
=============================================================================
SPATIAL CONTAMINATION CHECK (TPF & GAIA STAR MAP)
=============================================================================

DESCRIPTION:
This script performs a deep spatial analysis to detect potential contamination
from background stars (binaries) that might mimic a planet transit.

IT GENERATES TWO VISUALS:
1. TESS TPF IMAGE: The raw pixels seen by the telescope.
2. GAIA STAR MAP (Top-Down View): 
   - A precise reconstruction of the star field using Gaia DR3 data.
   - Shows relative positions in arcseconds (offsets from center).
   - Circle sizes are proportional to brightness (larger = brighter).
   - Reference circles show the size of a single TESS pixel (21").

TARGET:
- ID:     TOI 7475.01 / TIC 376866659
- Output: Professional plots and a console report of nearby stars.
"""

import lightkurve as lk
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.vizier import Vizier
import numpy as np
import warnings

# Ignore technical warnings for a cleaner output
warnings.simplefilter('ignore')

# 1. SET YOUR TARGET COORDINATES
# Replace these with the specific coordinates of your candidate
target_coords = SkyCoord("12:06:42.386 +00:16:26.83", unit=(u.hourangle, u.deg))

print("📡 Downloading Target Pixel Files (TPF) from NASA...")

# 2. SEARCH FOR TARGET PIXEL FILES (TPF)
search_result = lk.search_targetpixelfile(target_coords, radius=60)

if len(search_result) > 0:
    print(f"✅ Found in {len(search_result)} sectors!")
    
    # Take the first available sector
    tpf = search_result[0].download()
    
    print(f"📊 Sector: {tpf.sector if hasattr(tpf, 'sector') else 'N/A'}")
    print(f"🔭 Mission: {tpf.mission}")
    
    # 3. DRAW THE STAR FIELD
    fig, ax = plt.subplots(1, 2, figsize=(16, 7))
    
    # ===== IMAGE 1: RAW TPF PIXELS =====
    tpf.plot(ax=ax[0], aperture_mask=tpf.pipeline_mask, title="Pipeline Aperture Mask")
    
    # ===== IMAGE 2: PROFESSIONAL STAR MAP =====
    # Styling for a "Space" look
    ax[1].set_facecolor('#0a0e27')  # Dark blue background
    ax[1].grid(True, alpha=0.2, color='white', linestyle='--')
    ax[1].set_xlabel('Offset RA (arcsec)', fontsize=12, color='white')
    ax[1].set_ylabel('Offset Dec (arcsec)', fontsize=12, color='white')
    ax[1].set_title('🌌 STAR MAP - Top-Down View', fontsize=14, color='white', pad=20)
    ax[1].tick_params(colors='white')
    
    # ===== QUERY GAIA CATALOG =====
    try:
        print("\n⏳ Querying Gaia DR3 catalog via Vizier...")
        
        # Configure Vizier to return all columns and limit rows
        v = Vizier(columns=["*", "+_r"], row_limit=50)
        v.ROW_LIMIT = 50
        
        # Query region: 1 arcminute radius around target
        result = v.query_region(target_coords, radius=1.0*u.arcmin, catalog="I/355/gaiadr3")
        
        if len(result) > 0:
            gaia_table = result[0]
            print(f"✅ Found {len(gaia_table)} Gaia sources within 1 arcmin")
            
            # Arrays to store plotting data
            offsets_ra = []
            offsets_dec = []
            magnitudes = []
            ruwes = []
            
            # Calculate offsets for each star relative to the target
            for star in gaia_table:
                star_coord = SkyCoord(ra=star['RA_ICRS']*u.deg, dec=star['DE_ICRS']*u.deg)
                
                # Calculate offset in arcseconds
                # Correct RA for declination (cosine correction)
                delta_ra = (star_coord.ra.deg - target_coords.ra.deg) * 3600 * np.cos(np.radians(target_coords.dec.deg))
                delta_dec = (star_coord.dec.deg - target_coords.dec.deg) * 3600
                
                offsets_ra.append(delta_ra)
                offsets_dec.append(delta_dec)
                # Use Gmag if available, else assume faint (15)
                magnitudes.append(star['Gmag'] if 'Gmag' in star.colnames else 15)
                # RUWE (Renormalized Unit Weight Error) checks for binary instability
                ruwes.append(star['RUWE'] if 'RUWE' in star.colnames else 1.0)
            
            # Convert to numpy arrays for vector operations
            offsets_ra = np.array(offsets_ra)
            offsets_dec = np.array(offsets_dec)
            magnitudes = np.array(magnitudes)
            ruwes = np.array(ruwes)
            
            # Convert magnitudes to circle sizes (Logarithmic scale)
            # Smaller magnitude = Brighter star = Larger circle
            sizes = 1000 * 10**(-magnitudes/5)
            
            # Color map based on brightness
            colors_stars = plt.cm.plasma((magnitudes - magnitudes.min()) / (magnitudes.max() - magnitudes.min() + 0.01))
            
            # Draw the stars
            for i, (ra, dec, mag, size, color, ruwe) in enumerate(zip(offsets_ra, offsets_dec, magnitudes, sizes, colors_stars, ruwes)):
                
                # CASE A: The Target Star (Brightest one, usually index 0 but checked via min mag)
                if i == np.argmin(magnitudes):
                    ax[1].scatter(ra, dec, s=size*2, c='gold', marker='*', 
                                 edgecolors='yellow', linewidths=3, zorder=100,
                                 label=f'Target (Mag {mag:.1f})')
                    # Label "OUR STAR"
                    ax[1].text(ra, dec-3, 'OUR STAR', ha='center', va='top', 
                              color='yellow', fontsize=10, fontweight='bold')
                
                # CASE B: Background Stars
                else:
                    # Highlight high RUWE (possible binaries) with red edge
                    edge_color = 'red' if ruwe > 1.4 else 'white'
                    
                    ax[1].scatter(ra, dec, s=size, c=[color], marker='o', 
                                 edgecolors=edge_color, linewidths=2, alpha=0.8, zorder=50)
                    
                    # Label magnitude for bright stars only
                    if mag < 16:
                        ax[1].text(ra, dec+2, f'{mag:.1f}', ha='center', va='bottom',
                                  color='white', fontsize=8, alpha=0.7)
            
            # ===== REFERENCE CIRCLE: 1 TESS PIXEL =====
            # TESS pixels are huge (21 arcseconds per side)
            circle_tess = plt.Circle((0, 0), 21, color='red', fill=False, 
                                    linewidth=3, linestyle='--', alpha=0.7,
                                    label='1 TESS Pixel (21")')
            ax[1].add_patch(circle_tess)
            
            # ===== REFERENCE CIRCLE: TYPICAL APERTURE =====
            # The photometric aperture is often smaller
            circle_aperture = plt.Circle((0, 0), 10, color='lime', fill=False,
                                        linewidth=2, linestyle=':', alpha=0.5,
                                        label='Typical Aperture (~10")')
            ax[1].add_patch(circle_aperture)
            
            # Adjust plot limits to see everything clearly
            max_offset = max(abs(offsets_ra).max(), abs(offsets_dec).max()) * 1.2
            max_offset = max(max_offset, 30)  # Minimum view of 30 arcsec
            ax[1].set_xlim(-max_offset, max_offset)
            ax[1].set_ylim(-max_offset, max_offset)
            ax[1].set_aspect('equal')
            
            # Legend
            ax[1].legend(loc='upper right', fontsize=9, facecolor='#0a0e27', 
                        edgecolor='white', labelcolor='white')
            
            print("✅ Star map generated successfully!")
            
            # ===== CONSOLE REPORT =====
            print("\n--- 🌟 NEARBY STARS (Gaia G Magnitude) ---")
            
            # Sort stars by brightness for the report
            gaia_sorted = gaia_table[magnitudes.argsort()][:5]
            ra_sorted = offsets_ra[magnitudes.argsort()][:5]
            dec_sorted = offsets_dec[magnitudes.argsort()][:5]
            
            for i, (star, offset_ra, offset_dec) in enumerate(zip(gaia_sorted, ra_sorted, dec_sorted)):
                mag = star['Gmag'] if 'Gmag' in star.colnames else 0
                separation = np.sqrt(offset_ra**2 + offset_dec**2)
                ruwe = star['RUWE'] if 'RUWE' in star.colnames else 0
                
                symbol = "🌟" if i == 0 else "⭐"
                warning = ""
                
                # Warning if a bright star is inside the TESS pixel
                if separation < 21 and mag < 15 and i > 0:
                    warning = " ⚠️  INSIDE 1 TESS PIXEL!"
                
                print(f"{symbol} {i+1}. Mag G = {mag:.2f} | Sep = {separation:.1f}\" | "
                      f"RUWE = {ruwe:.2f}{warning}")
        
        else:
            ax[1].text(0.5, 0.5, '❌ No Gaia sources found', 
                      ha='center', va='center', transform=ax[1].transAxes,
                      color='white', fontsize=14)
    
    except Exception as e:
        print(f"❌ Error querying Gaia: {e}")
        ax[1].text(0.5, 0.5, f'❌ Error: {str(e)}', 
                  ha='center', va='center', transform=ax[1].transAxes,
                  color='red', fontsize=12)
        import traceback
        traceback.print_exc()

    plt.tight_layout()
    plt.show()
    
    # 4. FINAL DIAGNOSTIC GUIDE
    print("\n" + "="*70)
    print("--- 🔍 HOW TO INTERPRET THE STAR MAP ---")
    print("="*70)
    print("🌟 GOLD STAR   = Your Target (Center)")
    print("⭐ WHITE DOTS  = Nearby Background Stars")
    print("🔴 RED DASHED  = Size of 1 TESS Pixel (21 arcsec)")
    print("🟢 GREEN DOT   = Typical Aperture (~10 arcsec)")
    print("")
    print("✅ CLEAN FIELD: Only the Gold Star is inside the Green Circle.")
    print("⚠️  CONTAMINATION: Bright White Stars exist inside the Red Circle.")
    print("")
    print("💡 Dot Size = Brightness (Larger dot = Brighter star)")
    print("="*70)

else:
    print("❌ No Target Pixel File (TPF) data found for these coordinates at MAST.")