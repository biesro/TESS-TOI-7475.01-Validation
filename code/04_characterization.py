"""
EXOPLANET FULL CHARACTERIZATION PIPELINE
=========================================

Performs a complete photometric characterization of a TESS hot Jupiter candidate:
  1. Downloads and preprocesses SPOC light curves via lightkurve
  2. Pre-fit diagnostics: full light curve + phase-folded transit check
  3. Bayesian transit fit using Juliet + Dynesty nested sampling
  4. Monte Carlo error propagation combining fit posteriors and stellar uncertainties
  5. Derived physical parameters: Rp, Mp (Chen-Kipping), Teq, a, rho_p, g_p, H_atm, TSM
  6. Output plots and summary report

Features:
  - Fit cache system: if posteriors.hdf5 exists, skips re-running Dynesty (~70 min)
  - MAP (Maximum A Posteriori) estimates via KDE alongside median/percentiles
  - Classification with uncertainty ranges for both radius and mass
  - Population comparison against NASA Exoplanet Archive

Dependencies:
  numpy, matplotlib, corner, juliet, lightkurve, astropy, mr_forecast, scipy

Author : Biel Escolà Rodrigo
Date   : February 2025
"""

# ============================================================================
# WINDOWS FIX
# ============================================================================
import os
if 'HOME' not in os.environ:
    os.environ['HOME'] = os.environ.get('USERPROFILE', os.path.expanduser("~"))

# ============================================================================
# IMPORTS
# ============================================================================
import numpy as np
import matplotlib.pyplot as plt
import corner
from pathlib import Path
import shutil
import warnings
warnings.filterwarnings('ignore')

import juliet
import lightkurve as lk
import astropy.units as u
from astropy import constants as const
import mr_forecast as mr

# ============================================================================
# CONFIGURATION
# ============================================================================

class Config:
    # TARGET
    tic_id = "376866659"
    toi    = "7475.01"

    # STELLAR PARAMETERS  (from spectroscopy / TIC)
    teff      = 6467        # Effective temperature [K]
    teff_err  = 52
    logg      = 4.136       # Surface gravity [cgs]
    logg_err  = 0.1
    r_star     = 1.7257     # Stellar radius [R_sun]
    r_star_err = 0.10
    m_star     = 1.279      # Stellar mass [M_sun]
    m_star_err = 0.20
    feh     = -0.099        # Metallicity [Fe/H]
    feh_err =  0.08
    mag_J   = 7.676         # 2MASS J-band magnitude (verified via astroquery)

    # ORBITAL PARAMETERS  (initial estimates from BLS / ExoFOP)
    period     = 3.253773   # Orbital period [days]
    period_err = 0.0001
    epoch_bjd  = 2460775.5819  # Transit epoch [BJD]
    epoch_err  = 0.005
    depth_ppm      = 4601.0    # Transit depth [ppm]
    duration_hours = 4.08      # Transit duration [hours]

    # PREPROCESSING
    binning_min   = 3       # Light curve binning [minutes]
    use_windows   = True    # Use transit windows only (reduces data volume)
    window_factor = 1.5     # Window half-width = window_factor × duration

    # DYNESTY SETTINGS
    nlive    = 1000         # Number of live points
    dlogz    = 0.1          # Convergence criterion
    nthreads = 4            # Parallel threads

    # MONTE CARLO
    n_samples = 10000       # Number of MC samples for error propagation

    # CACHE — if True and posteriors.hdf5 exists, skips re-running Dynesty
    # Set to False or delete the output folder to force a new fit
    use_cached_fit = True
    out_folder     = 'results_characterization'


# ============================================================================
# PHYSICAL FUNCTIONS
# ============================================================================

def calc_tsm(Rp_Rearth, Mp_Mearth, R_star_Rsol, Teq, mag_J):
    """
    Transmission Spectroscopy Metric (TSM).
    Kempton et al. (2018), PASP 130, 114401.

    Parameters
    ----------
    Rp_Rearth   : float — planet radius [R_Earth]
    Mp_Mearth   : float — planet mass [M_Earth]
    R_star_Rsol : float — stellar radius [R_sun]
    Teq         : float — equilibrium temperature [K]
    mag_J       : float — 2MASS J-band magnitude

    Returns
    -------
    float — TSM value
    """
    if   Rp_Rearth < 1.5:  scale = 0.190
    elif Rp_Rearth < 2.75: scale = 1.260
    elif Rp_Rearth < 4.0:  scale = 1.280
    else:                   scale = 0.167
    tsm = scale * (Rp_Rearth**3 * Teq) / (Mp_Mearth * R_star_Rsol**2)
    return tsm * 10.0**(-mag_J / 5.0)


def classify_planet(Rp_Rjup, Teq, Mp_Mjup=None):
    """
    Classify a planet by size, temperature regime, and mass.

    Parameters
    ----------
    Rp_Rjup  : float — planet radius [R_Jup]
    Teq      : float — equilibrium temperature [K]
    Mp_Mjup  : float or None — planet mass [M_Jup], optional

    Returns
    -------
    dict with keys: size, size_emoji, temp, temp_emoji, temp_note,
                    mass (if Mp provided), overall
    """
    result = {}

    # Size classification
    if   Rp_Rjup < 0.5: result['size'] = 'Sub-Saturn';             result['size_emoji'] = '🪐'
    elif Rp_Rjup < 0.9: result['size'] = 'Saturn-sized';           result['size_emoji'] = '🪐'
    elif Rp_Rjup < 1.3: result['size'] = 'Jupiter-sized';          result['size_emoji'] = '♃'
    else:                result['size'] = 'Super-Jupiter (inflated)'; result['size_emoji'] = '🎈'

    # Temperature classification
    if   Teq > 2200: result['temp'] = 'Ultra-Hot Jupiter'; result['temp_emoji'] = '🔥'; result['temp_note'] = 'Thermal dissociation'
    elif Teq > 1800: result['temp'] = 'Very Hot Jupiter';  result['temp_emoji'] = '🌡️'; result['temp_note'] = 'TiO/VO, thermal inversion'
    elif Teq > 1000: result['temp'] = 'Hot Jupiter';       result['temp_emoji'] = '☀️'; result['temp_note'] = 'Typical hot Jupiter'
    else:            result['temp'] = 'Warm Jupiter';      result['temp_emoji'] = '🌤️'; result['temp_note'] = 'Cooler giant'

    # Mass classification (optional)
    if Mp_Mjup is not None:
        if   Mp_Mjup < 0.5: result['mass'] = 'Low-mass giant'
        elif Mp_Mjup < 2.0: result['mass'] = 'Jupiter-mass'
        else:                result['mass'] = 'Super-Jupiter'

    result['overall'] = f"{result['temp']}, {result['size']}"
    return result


# ============================================================================
# PRE-FIT DIAGNOSTICS
# ============================================================================

def plot_full_lc_with_transits(lc_raw, config):
    """
    Plot the full light curve with expected transit times marked in red.

    Verifies that the epoch is correct BEFORE running the Bayesian fit.
    If the red lines do not coincide with flux dips → epoch is wrong.

    Parameters
    ----------
    lc_raw : lightkurve LightCurve object
    config : Config instance
    """
    epoch_btjd = config.epoch_bjd - 2457000.0
    t = lc_raw.time.value
    f = lc_raw.flux.value

    fig, ax = plt.subplots(figsize=(20, 4))
    ax.scatter(t, f, s=0.3, alpha=0.2, c='steelblue')

    # Mark each expected transit time
    n_start = int((t.min() - epoch_btjd) / config.period) - 1
    n_end   = int((t.max() - epoch_btjd) / config.period) + 2
    for n in range(n_start, n_end):
        t_tr = epoch_btjd + n * config.period
        if t.min() <= t_tr <= t.max():
            ax.axvline(t_tr, color='red', alpha=0.6, linewidth=1.2)

    ax.set_xlabel('BTJD [days]', fontsize=11)
    ax.set_ylabel('Normalized flux', fontsize=11)
    ax.set_title(f'TOI {config.toi} — Full light curve + expected transits (red)\n'
                 f'Verify that red lines coincide with flux dips!',
                 fontsize=11)
    ax.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig('diagnostic_full_lc.png', dpi=120, bbox_inches='tight')
    plt.close()
    print("  → diagnostic_full_lc.png  ← CHECK that red lines coincide with flux dips!")


def plot_phased_check(times, flux, config):
    """
    Phase-folded light curve PRE-FIT with:
      - Individual data points (blue, transparent)
      - 10-min binned data (black) to make the transit visible
      - Auto-scaled Y axis to show the real transit depth
      - Expected depth reference line

    Also prints a numerical flux diagnostic to verify epoch and normalization.

    Parameters
    ----------
    times  : array — time array [BTJD]
    flux   : array — normalized flux
    config : Config instance
    """
    epoch_btjd = config.epoch_bjd - 2457000.0
    phase = ((times - epoch_btjd) % config.period) / config.period
    phase = np.where(phase > 0.5, phase - 1, phase)
    phase_h = phase * config.period * 24   # convert to hours

    p_expected = np.sqrt(config.depth_ppm / 1e6)

    # --- Manual 10-min binning to reveal the transit ---
    bin_width_h = 10 / 60   # 10 minutes in hours
    bins = np.arange(phase_h.min(), phase_h.max() + bin_width_h, bin_width_h)
    bin_centers, bin_flux, bin_err = [], [], []
    for i in range(len(bins) - 1):
        mask = (phase_h >= bins[i]) & (phase_h < bins[i+1])
        if mask.sum() >= 3:
            bin_centers.append((bins[i] + bins[i+1]) / 2)
            bin_flux.append(np.median(flux[mask]))          # median: robust to outliers
            bin_err.append(np.std(flux[mask]) / np.sqrt(mask.sum()))  # standard error
    bin_centers = np.array(bin_centers)
    bin_flux    = np.array(bin_flux)
    bin_err     = np.array(bin_err)

    # Y-axis: 30% margin below transit depth and above baseline
    depth_val  = 1.0 - config.depth_ppm / 1e6
    flux_range = 1.0 - depth_val
    y_min = depth_val - flux_range * 0.3
    y_max = flux.max() + flux_range * 0.3

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    fig.suptitle(f'TOI {config.toi} — Pre-fit verification\n'
                 f'Expected p = {p_expected:.4f}   '
                 f'Expected Rp ≈ {p_expected * config.r_star * 9.95492:.3f} RJup',
                 fontsize=12)

    # Upper panel: zoom around the expected transit
    ax = axes[0]
    ax.scatter(phase_h, flux, s=2, alpha=0.15, c='steelblue', label='Individual points')
    if len(bin_centers) > 0:
        ax.errorbar(bin_centers, bin_flux, yerr=bin_err,
                    fmt='o', color='black', markersize=4, linewidth=1.2,
                    capsize=2, label='10-min bins', zorder=5)
    ax.axvspan(-config.duration_hours/2, config.duration_hours/2,
               alpha=0.12, color='red', label=f'Expected transit ({config.duration_hours:.1f} h)')
    ax.axhline(depth_val, color='red', linestyle='--', alpha=0.8,
               label=f'Expected depth ({config.depth_ppm:.0f} ppm)')
    ax.axhline(1.0, color='gray', linestyle=':', alpha=0.5)
    ax.set_ylim(y_min, y_max)
    ax.set_xlim(-config.duration_hours * 2.5, config.duration_hours * 2.5)
    ax.set_ylabel('Normalized flux', fontsize=11)
    ax.legend(fontsize=9, loc='lower right')
    ax.grid(alpha=0.25)
    ax.set_title('Transit zoom', fontsize=10)

    # Lower panel: full phase view
    ax2 = axes[1]
    ax2.scatter(phase_h, flux, s=1, alpha=0.1, c='steelblue')
    if len(bin_centers) > 0:
        ax2.errorbar(bin_centers, bin_flux, yerr=bin_err,
                     fmt='o', color='black', markersize=2, linewidth=0.8,
                     capsize=1, zorder=5)
    ax2.axvspan(-config.duration_hours/2, config.duration_hours/2,
                alpha=0.12, color='red')
    ax2.axhline(depth_val, color='red', linestyle='--', alpha=0.8)
    ax2.set_xlabel('Phase × Period [hours]', fontsize=11)
    ax2.set_ylabel('Normalized flux', fontsize=11)
    ax2.set_title('Full phase', fontsize=10)
    ax2.grid(alpha=0.25)

    plt.tight_layout()
    plt.savefig('phased_check.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Numerical flux diagnostic
    print(f"\n  FLUX DIAGNOSTIC:")
    print(f"    Median flux (total)    : {np.median(flux):.6f}")
    print(f"    Minimum flux           : {flux.min():.6f}")
    print(f"    Maximum flux           : {flux.max():.6f}")
    print(f"    Expected depth         : {config.depth_ppm/1e6:.6f}  ({config.depth_ppm:.0f} ppm)")
    print(f"    Expected in-transit flux: {1.0 - config.depth_ppm/1e6:.6f}")

    in_transit_mask = np.abs(phase_h) < config.duration_hours / 2
    n_in_transit = in_transit_mask.sum()
    if n_in_transit > 0:
        flux_in_transit = flux[in_transit_mask]
        print(f"    Points in-transit      : {n_in_transit}")
        print(f"    Median in-transit flux : {np.median(flux_in_transit):.6f}")
        delta = np.median(flux_in_transit) - np.median(flux[~in_transit_mask])
        print(f"    Observed Δflux         : {delta*1e6:.0f} ppm  (expected: -{config.depth_ppm:.0f} ppm)")
        if abs(delta) < config.depth_ppm / 1e6 * 0.1:
            print(f"    ⚠️  TRANSIT NOT VISIBLE — possible epoch or normalization issue!")
        else:
            print(f"    ✓ Transit detected in phase-folded curve")
    else:
        print(f"    ⚠️  No points inside transit window — epoch may be wrong!")

    print(f"  → phased_check.png")


# ============================================================================
# BAYESIAN FIT WITH CACHE
# ============================================================================

def run_juliet_fit(times, flux, flux_err, config):
    """
    Run or load the Bayesian transit fit using Juliet + Dynesty.

    If config.use_cached_fit=True and posteriors.hdf5 already exists,
    loads posteriors without re-running Dynesty (~70 min saved).
    To force a new fit: set use_cached_fit=False or delete the output folder.

    Parameters
    ----------
    times    : array — time array [BTJD]
    flux     : array — normalized flux
    flux_err : array — flux uncertainties
    config   : Config instance

    Returns
    -------
    juliet results object
    """
    print("\n" + "="*70)
    print("BAYESIAN FIT — JULIET + DYNESTY")
    print("="*70)

    epoch_btjd  = config.epoch_bjd - 2457000.0
    p_expected  = np.sqrt(config.depth_ppm / 1e6)   # depth = p² → p = sqrt(depth)
    p_lo = max(0.01,  p_expected * 0.5)
    p_hi = min(0.25,  p_expected * 1.5)

    tim = {'TESS': np.array(times)}
    fl  = {'TESS': np.array(flux)}
    fle = {'TESS': np.array(flux_err)}

    priors = {
        # Limb darkening — quadratic law, Kipping (2013) parametrization
        'ld_law_TESS': {'distribution': 'fixed',      'hyperparameters': 'quadratic'},
        'q1_TESS':     {'distribution': 'normal',     'hyperparameters': [0.348, 0.05]},  # Claret (2017)
        'q2_TESS':     {'distribution': 'normal',     'hyperparameters': [0.246, 0.05]},  # Claret (2017)

        # Orbital parameters — Gaussian priors from BLS solution
        'P_p1':        {'distribution': 'normal',     'hyperparameters': [config.period, config.period_err]},
        't0_p1':       {'distribution': 'normal',     'hyperparameters': [epoch_btjd, config.epoch_err]},

        # Transit shape — uniform priors
        'p_p1':        {'distribution': 'uniform',    'hyperparameters': [p_lo, p_hi]},    # Rp/R*
        'b_p1':        {'distribution': 'uniform',    'hyperparameters': [0.0, 0.99]},     # impact parameter

        # Stellar density — used to compute a/R* via Kepler's 3rd law
        # rho = M* / (4/3 * pi * R*^3), computed from stellar parameters
        'rho':         {'distribution': 'normal',     'hyperparameters': [350.914, 93.149]},

        # Eccentricity — fixed to 0 (tidal circularization, P < 10 days)
        'ecc_p1':      {'distribution': 'fixed',      'hyperparameters': 0.0},
        'omega_p1':    {'distribution': 'fixed',      'hyperparameters': 90.0},

        # Instrumental parameters
        'dilution_TESS':  {'distribution': 'fixed',      'hyperparameters': 1.0},
        'mdilution_TESS': {'distribution': 'fixed',      'hyperparameters': 0.0},
        'mflux_TESS':     {'distribution': 'normal',     'hyperparameters': [0, 0.1]},      # flux offset
        'sigma_w_TESS':   {'distribution': 'loguniform', 'hyperparameters': [0.1, 1000]},   # white noise [ppm]
    }

    # Check if a previous fit exists
    cached = (Path(config.out_folder) / 'posteriors.hdf5').exists()

    if config.use_cached_fit and cached:
        print(f"\n  ✓ CACHE HIT — loading existing fit from '{config.out_folder}/'")
        print(f"    (To re-run: set use_cached_fit=False or delete the folder)")
    else:
        if cached:
            print(f"\n  Deleting previous fit and re-running...")
            shutil.rmtree(config.out_folder)
        print(f"\n  Prior p_p1     : Uniform [{p_lo:.4f}, {p_hi:.4f}]")
        print(f"  Eccentricity   : FIXED at 0")
        print(f"  Live points    : {config.nlive}")
        print(f"  Running Dynesty nested sampling...")

    dataset = juliet.load(
        priors=priors,
        t_lc=tim, y_lc=fl, yerr_lc=fle,
        out_folder=config.out_folder
    )
    results = dataset.fit(
        sampler='dynesty',
        nthreads=config.nthreads,
        n_live_points=config.nlive,
        dlogz=config.dlogz
    )

    if config.use_cached_fit and cached:
        print("  ✓ Posteriors loaded successfully")
    else:
        print("  ✓ Fit completed")

    # Corner plot of fit posteriors
    try:
        post   = results.posteriors['posterior_samples']
        params = [p for p in ['p_p1', 'b_p1', 'rho', 't0_p1', 'P_p1', 'q1_TESS', 'q2_TESS']
                  if p in post]
        labels_map = {
            'p_p1':    r'$R_p/R_*$',    'b_p1':    r'$b$',
            'rho':     r'$\rho_*$',     't0_p1':   r'$T_0$',
            'P_p1':    r'$P$ [d]',      'q1_TESS': r'$q_1$',
            'q2_TESS': r'$q_2$',
        }
        data_c = np.column_stack([post[p] for p in params])
        lbls   = [labels_map.get(p, p) for p in params]
        fig = corner.corner(data_c, labels=lbls,
                            quantiles=[0.16, 0.5, 0.84],
                            show_titles=True,
                            title_kwargs={"fontsize": 10},
                            label_kwargs={"fontsize": 11})
        fig.savefig('corner_characterization.png', dpi=120, bbox_inches='tight')
        plt.close(fig)
        print("  → corner_characterization.png")
    except Exception as e:
        print(f"  ⚠️  Corner plot failed: {e}")

    return results


# ============================================================================
# MONTE CARLO ERROR PROPAGATION
# ============================================================================

def propagate_uncertainties(results, config):
    """
    Derive physical parameters and propagate uncertainties via Monte Carlo.

    Generates n_samples realizations by drawing from:
      - Gaussian distributions for stellar parameters (M*, R*, Teff)
      - Posterior samples from the Juliet fit (p, b, rho, P, T0)

    This correctly accounts for the correlation between fit posteriors
    and stellar uncertainties, which analytical error propagation cannot handle.

    Parameters
    ----------
    results : juliet results object
    config  : Config instance

    Returns
    -------
    dict of arrays, one per derived parameter, each of length n_samples
    """
    print("\n" + "="*70)
    print("MONTE CARLO ERROR PROPAGATION")
    print("="*70)
    n = config.n_samples
    print(f"\n  Drawing {n} samples...")

    # --- Stellar parameter samples ---
    M_star_s = np.random.normal(config.m_star,    config.m_star_err, n)
    R_star_s = np.random.normal(config.r_star,    config.r_star_err, n)
    Teff_s   = np.random.normal(config.teff,      config.teff_err,   n)

    # --- Fit posterior samples (drawn directly, preserving true distribution shape) ---
    post = results.posteriors['posterior_samples']
    idx  = np.random.randint(0, len(post['p_p1']), n)
    P_s  = post['P_p1'][idx]
    T0_s = post['t0_p1'][idx] + 2457000.0   # convert BTJD → BJD
    p_s   = post['p_p1'][idx]
    b_s   = post['b_p1'][idx]
    rho_s = np.clip(post['rho'][idx], 10.0, None)  # clip to avoid divide-by-zero

    # --- Semi-major axis via Kepler's 3rd law reformulated with stellar density ---
    # a/R* = (G × rho* × P² / 3π)^(1/3)
    G_SI   = const.G.value
    P_sec  = P_s * 86400.0
    a_rs_s = (G_SI * rho_s * P_sec**2 / (3.0 * np.pi))**(1.0 / 3.0)

    # --- Planet radius ---
    # p = Rp/R*  →  Rp = p × R* × conversion_factor
    Rp_Rjup   = p_s * R_star_s * 9.95492    # R_sun / R_Jup (equatorial)
    Rp_Rearth = Rp_Rjup * 11.2089           # R_Jup / R_Earth (equatorial)
    a_au      = a_rs_s * R_star_s * 0.00465047  # R_sun → AU conversion

    # --- Equilibrium temperature ---
    # Teq = Teff × (R*/(2a))^(1/2) × (1-A)^(1/4)
    # Assumes Bond albedo A=0.3 and heat redistribution factor f=0.5
    R_sun_to_au = (const.R_sun / const.au).decompose().value
    Teq = Teff_s * np.sqrt(R_star_s * R_sun_to_au / (2.0 * a_au)) * (0.5 * 0.7)**0.25

    # --- Incident flux relative to Earth ---
    # F = sigma × Teff^4 × (R*/a)²  /  F_sun
    sigma    = const.sigma_sb.cgs.value
    R_sun_cm = const.R_sun.to(u.cm).value
    au_cm    = const.au.to(u.cm).value
    F_rel    = sigma * Teff_s**4 * (R_star_s * R_sun_cm / (a_au * au_cm))**2 / 1.36e6

    # --- Planet mass via Chen-Kipping (2017) M-R relation ---
    print("  → Computing masses (mr_forecast / Chen-Kipping 2017)...")
    Rp_safe   = np.clip(Rp_Rearth, 0.11, 99.0)   # valid domain of the model
    Mp_Mearth = np.array(mr.Rpost2M(Rp_safe, unit='Earth'), dtype=float)
    Mp_Mearth = np.where(Mp_Mearth > 0, Mp_Mearth, np.nan)  # remove numerical artifacts
    Mp_Mjup   = Mp_Mearth / 317.8

    # --- Derived parameters ---
    print("  → Computing density, gravity, and atmospheric scale height...")
    M_jup_kg = const.M_jup.value
    R_jup_m  = const.R_jup.value
    k_B      = const.k_B.value
    m_H      = const.m_p.value

    rho_p  = (Mp_Mjup * M_jup_kg) / ((4.0/3.0) * np.pi * (Rp_Rjup * R_jup_m)**3)
    g_p    = G_SI * (Mp_Mjup * M_jup_kg) / (Rp_Rjup * R_jup_m)**2
    logg_p = np.log10(g_p * 100.0)       # log g in cgs [cm/s²]
    v_esc  = np.sqrt(2.0 * G_SI * (Mp_Mjup * M_jup_kg) / (Rp_Rjup * R_jup_m)) / 1000.0
    H_atm  = (k_B * Teq) / (2.3 * m_H * g_p) / 1000.0   # mu=2.3 for H2-dominated atmosphere [km]

    print("  ✓ Propagation complete")
    return {
        'Rp_Rjup':    Rp_Rjup,   'Rp_Rearth': Rp_Rearth,  'a_au':      a_au,
        'Teq':        Teq,        'Mp_Mjup':   Mp_Mjup,     'Mp_Mearth': Mp_Mearth,
        'rho_planet': rho_p,      'g_planet':  g_p,         'logg_planet': logg_p,
        'v_esc':      v_esc,      'H_atm':     H_atm,       'F_incident': F_rel,
        'b_impact':   b_s,        'P_orbital': P_s,          'T0_bjd':    T0_s,
        'eccentricity': np.zeros(n),
    }


# ============================================================================
# RESULTS — HELPER FUNCTIONS
# ============================================================================

def _map(array, name, unit, decimals=3):
    """
    Compute and print the MAP (Maximum A Posteriori) via Gaussian KDE.
    More informative than the median for asymmetric distributions (e.g. Mp).

    Parameters
    ----------
    array    : array — posterior samples
    name     : str   — parameter label for printing
    unit     : str   — unit string
    decimals : int   — decimal places

    Returns
    -------
    float — MAP value
    """
    from scipy.stats import gaussian_kde
    clean = array[np.isfinite(array)]
    kde = gaussian_kde(clean)
    x = np.linspace(np.percentile(clean, 1), np.percentile(clean, 99), 1000)
    map_val = x[np.argmax(kde(x))]
    fmt = f"{{:.{decimals}f}}"
    print(f"  {name:32s} = {fmt.format(map_val)} {unit}  [MAP]")
    return map_val


def _pct(array, name, unit, decimals=3):
    """
    Compute and print median ± 1-sigma (16th/84th percentiles).

    Parameters
    ----------
    array    : array — posterior samples
    name     : str   — parameter label for printing
    unit     : str   — unit string
    decimals : int   — decimal places

    Returns
    -------
    tuple (median, upper_error, lower_error)
    """
    p16, p50, p84 = np.nanpercentile(array, [16, 50, 84])
    ep, em = p84 - p50, p50 - p16
    fmt = f"{{:.{decimals}f}}"
    print(f"  {name:32s} = {fmt.format(p50)} +{fmt.format(ep)} -{fmt.format(em)} {unit}")
    return p50, ep, em


# ============================================================================
# RESULTS — MAIN OUTPUT
# ============================================================================

def print_results(mc, config):
    """
    Print the complete characterization summary to stdout.

    Parameters
    ----------
    mc     : dict — Monte Carlo output from propagate_uncertainties()
    config : Config instance

    Returns
    -------
    tuple (cls, Rp_med, Mp_med, Teq_med, rho_med)
    """
    print("\n" + "="*70)
    print("FULL CHARACTERIZATION RESULTS")
    print("="*70)

    print("\n[1] STELLAR PARAMETERS:")
    print(f"  TIC {config.tic_id} / TOI {config.toi}")
    print(f"  Teff   = {config.teff:.0f} ± {config.teff_err:.0f} K")
    print(f"  M*     = {config.m_star:.3f} ± {config.m_star_err:.3f} M_sun")
    print(f"  R*     = {config.r_star:.4f} ± {config.r_star_err:.4f} R_sun")
    print(f"  log g  = {config.logg:.2f} ± {config.logg_err:.2f}")
    print(f"  [Fe/H] = {config.feh:.2f} ± {config.feh_err:.2f}")

    print("\n[2] ORBITAL PARAMETERS:")
    _pct(mc['P_orbital'], "P",               "days", decimals=7)
    _pct(mc['T0_bjd'],    "T0",              "BJD",  decimals=6)
    print(f"  e  = 0 (fixed — HJ, tidally circularized)")
    a_med, *_ = _pct(mc['a_au'],     "a (median)",      "AU", 4)
    a_map     = _map(mc['a_au'],     "a (MAP)",         "AU", 4)
    _pct(mc['b_impact'],  "b (impact parameter)", "",  3)

    print("\n[3] PLANETARY PHYSICAL PARAMETERS:")
    Rp_med, Rp_ep, Rp_em = _pct(mc['Rp_Rjup'],   "Rp (median)",   "R_Jup",   3)
    Rp_map                = _map(mc['Rp_Rjup'],   "Rp (MAP)",      "R_Jup",   3)
    Re_med, *_            = _pct(mc['Rp_Rearth'], "Rp (median)",   "R_Earth", 2)
    Re_map                = _map(mc['Rp_Rearth'], "Rp (MAP)",      "R_Earth", 2)
    Mp_med, Mp_ep, Mp_em  = _pct(mc['Mp_Mjup'],   "Mp (median)",   "M_Jup",   3)
    Mp_map                = _map(mc['Mp_Mjup'],   "Mp (MAP)",      "M_Jup",   3)
    Me_med, *_            = _pct(mc['Mp_Mearth'], "Mp (median)",   "M_Earth", 1)
    Me_map                = _map(mc['Mp_Mearth'], "Mp (MAP)",      "M_Earth", 1)
    rho_med, *_           = _pct(mc['rho_planet'],"rho_p (median)","kg/m³",   1)
    rho_map               = _map(mc['rho_planet'],"rho_p (MAP)",   "kg/m³",   1)

    print("\n[4] GRAVITY AND ESCAPE:")
    g_med,    *_ = _pct(mc['g_planet'],    "g_p",     "m/s²",  1)
    logg_med, *_ = _pct(mc['logg_planet'], "log g_p", "(cgs)", 2)
    v_med,    *_ = _pct(mc['v_esc'],       "v_esc",   "km/s",  1)

    print("\n[5] TEMPERATURE AND IRRADIATION:")
    Teq_med, *_ = _pct(mc['Teq'],        "T_eq (median)", "K",       0)
    Teq_map     = _map(mc['Teq'],        "T_eq (MAP)",    "K",       0)
    F_med,   *_ = _pct(mc['F_incident'], "F_p",           "F_Earth", 0)

    Teq_A0 = Teq_med * (1.0 / 0.7)**0.25
    Teq_A5 = Teq_med * (0.5 / 0.7)**0.25
    print(f"\n  T_eq range by albedo:")
    print(f"    A = 0.0 (blackbody):        {Teq_A0:.0f} K")
    print(f"    A = 0.3 (typical HJ):       {Teq_med:.0f} K")
    print(f"    A = 0.5 (highly reflective):{Teq_A5:.0f} K")

    print("\n[6] ATMOSPHERE:")
    H_med, *_ = _pct(mc['H_atm'], "H (median)", "km", 0)
    H_map     = _map(mc['H_atm'], "H (MAP)",    "km", 0)

    # TSM with forecaster mass and with reference 1 MJup
    TSM       = calc_tsm(Re_med, Me_med, config.r_star, Teq_med, config.mag_J)
    TSM_1Mjup = calc_tsm(Re_med, 318.0,  config.r_star, Teq_med, config.mag_J)
    print(f"\n  TSM (Kempton 2018, forecaster mass) ≈ {TSM:.1f}  (mJ = {config.mag_J})")
    print(f"  TSM assuming Mp = 1 MJup            ≈ {TSM_1Mjup:.1f}  (reference case)")
    if TSM_1Mjup > 90:
        print(f"    → ✓ EXCELLENT for JWST if Mp ~ 1 MJup")
    elif TSM_1Mjup > 40:
        print(f"    → Favorable for JWST if Mp ~ 1 MJup")
    else:
        print(f"    → Low TSM even for 1 MJup")

    print("\n[7] CLASSIFICATION:")

    # Based on median
    print("  Based on MEDIAN:")
    cls = classify_planet(Rp_med, Teq_med, Mp_med)
    print(f"  {cls['size_emoji']}  Radius:      {cls['size']}  ({Rp_med:.2f} RJup)")
    print(f"  {cls['temp_emoji']}  Temperature: {cls['temp']}")
    if 'mass' in cls:
        print(f"  ⚖️   Mass:        {cls['mass']}  ({Mp_med:.2f} MJup)")

    # Based on MAP
    print("\n  Based on MAP:")
    cls_map = classify_planet(Rp_map, Teq_map, Mp_map)
    print(f"  {cls_map['size_emoji']}  Radius:      {cls_map['size']}  ({Rp_map:.2f} RJup)")
    print(f"  {cls_map['temp_emoji']}  Temperature: {cls_map['temp']}")
    if 'mass' in cls_map:
        print(f"  ⚖️   Mass:        {cls_map['mass']}  ({Mp_map:.2f} MJup)")

    # Uncertainty ranges
    Rp_lo = Rp_med - Rp_em
    Rp_hi = Rp_med + Rp_ep
    Mp_lo = max(0, Mp_med - Mp_em)
    Mp_hi = Mp_med + Mp_ep
    print(f"\n  ⚠️  Rp 1-sigma range: {Rp_lo:.2f} – {Rp_hi:.2f} RJup")
    print(f"      Size classification may vary within this range")
    print(f"\n  ⚠️  Mp 1-sigma range: {Mp_lo:.2f} – {Mp_hi:.2f} MJup")
    print(f"      Highly asymmetric distribution — MAP={Mp_map:.2f} MJup is more representative")
    print(f"      Radial velocity follow-up required to constrain mass")

    print(f"\n[8] COMPARISON WITH JUPITER:")
    print(f"  Radius:      {Rp_med:.2f}× R_Jup")
    print(f"  Mass:        {Mp_med:.2f}× M_Jup (estimated, no RV)")
    print(f"  Density:     {rho_med/1326:.2f}× rho_Jup  ({rho_med:.0f} vs 1326 kg/m³)")
    print(f"  Gravity:     {g_med/24.79:.2f}× g_Jup     ({g_med:.1f} vs 24.8 m/s²)")
    print(f"  Temperature: {Teq_med/165:.0f}× T_Jup     ({Teq_med:.0f} vs 165 K)")
    print(f"  Irradiation: {F_med:.0f}× F_Earth")

    rho_ratio = rho_med / 1326
    if   rho_ratio < 0.3: print(f"\n  ⚠️  HIGHLY INFLATED PLANET (rho < 0.3 rho_Jup)")
    elif rho_ratio < 0.7: print(f"\n  ℹ️  Slightly inflated (typical for hot Jupiters)")
    else:                  print(f"\n  ℹ️  Normal density (not significantly inflated)")

    return cls, Rp_med, Mp_med, Teq_med, rho_med


# ============================================================================
# PLOTS
# ============================================================================

def plot_population_comparison(Rp_Rjup, Teq, config):
    """
    Compare this planet against the known hot Jupiter population
    from the NASA Exoplanet Archive.

    Downloads data on first run and caches locally as nasa_exoplanet_cache.csv.

    Parameters
    ----------
    Rp_Rjup : float — planet radius [R_Jup]
    Teq     : float — equilibrium temperature [K]
    config  : Config instance
    """
    import pandas as pd

    cache_file = 'nasa_exoplanet_cache.csv'

    if os.path.exists(cache_file):
        print("  → Loading exoplanet population from local cache...")
        df = pd.read_csv(cache_file)
    else:
        print("  → Downloading hot Jupiter population from NASA Exoplanet Archive...")
        try:
            from astroquery.nasa_exoplanet_archive import NasaExoplanetArchive
            raw = NasaExoplanetArchive.query_criteria(
                table   = "pscomppars",
                select  = "pl_name, pl_rade, pl_eqt, pl_orbper, pl_bmassj",
                where   = "pl_orbper < 10 AND pl_rade > 8 AND pl_eqt > 500"
            )
            df = raw.to_pandas()
            df.to_csv(cache_file, index=False)
            print(f"  → {len(df)} planets downloaded and saved to '{cache_file}'")
        except Exception as e:
            print(f"  ⚠️  NASA Archive download failed: {e}")
            print(f"  → Using fallback population (do not cite in paper)")
            df = pd.DataFrame({
                'pl_rade': [8.9, 10.7, 11.8, 12.7, 13.5, 14.3, 15.2,
                            16.1, 11.9, 13.8, 15.6, 16.6, 12.2, 14.1],
                'pl_eqt':  [1240, 1620, 1850, 1450, 1730, 2020, 2180,
                            1910, 1380, 1650, 1820, 2100, 1540, 1750],
            })
            df['pl_rade_jup'] = df['pl_rade'] / 11.2089
        else:
            df['pl_rade_jup'] = pd.to_numeric(df['pl_rade'], errors='coerce') / 11.2089

    df['pl_eqt']      = pd.to_numeric(df['pl_eqt'], errors='coerce')
    df['pl_rade_jup'] = pd.to_numeric(df.get('pl_rade_jup',
                         df.get('pl_rade', pd.Series()) / 11.2089), errors='coerce')

    mask = (
        df['pl_eqt'].notna() &
        df['pl_rade_jup'].notna() &
        (df['pl_eqt'] > 500) &
        (df['pl_rade_jup'] > 0.5) &
        (df['pl_rade_jup'] < 3.0)
    )
    df_clean = df[mask]
    print(f"  → {len(df_clean)} hot Jupiters with complete data")

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(
        df_clean['pl_eqt'], df_clean['pl_rade_jup'],
        c='lightgray', s=60, alpha=0.6,
        edgecolors='darkgray', linewidths=0.3,
        label=f'NASA Exoplanet Archive (N={len(df_clean)})',
        zorder=2
    )
    ax.scatter(
        Teq, Rp_Rjup,
        c='red', s=300, marker='*',
        edgecolors='darkred', linewidths=1.5,
        label=f'TOI {config.toi}',
        zorder=5
    )
    ax.axhline(1.3, color='orange', linestyle='--', alpha=0.5, linewidth=1.2)
    ax.text(0.02, 0.70, 'Inflated (Rp > 1.3 RJ)',
            transform=ax.transAxes, fontsize=9, color='darkorange', alpha=0.9)
    ax.axvline(1800, color='red', linestyle='--', alpha=0.5, linewidth=1.2)
    ax.text(0.72, 0.02, 'Ultra-Hot (T > 1800 K)',
            transform=ax.transAxes, fontsize=9, color='red', alpha=0.9)
    ax.set_xlabel('Equilibrium Temperature [K]', fontsize=13)
    ax.set_ylabel(r'Planet Radius [$R_{Jup}$]', fontsize=13)
    ax.set_title('Hot Jupiter Population — NASA Exoplanet Archive',
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(alpha=0.25, linestyle='--')
    plt.tight_layout()
    plt.savefig('population_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  → population_comparison.png")


def plot_posterior_distributions(mc, config):
    """
    Plot histograms of the six main derived posterior distributions
    with median and 1-sigma lines overlaid.

    Parameters
    ----------
    mc     : dict — Monte Carlo output from propagate_uncertainties()
    config : Config instance
    """
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    fig.suptitle(f'TOI {config.toi} — Derived posteriors', fontsize=14, fontweight='bold')
    params = [
        (mc['Rp_Rjup'],    r'$R_p$ [$R_{Jup}$]',  'steelblue'),
        (mc['Mp_Mjup'],    r'$M_p$ [$M_{Jup}$]',  'darkorange'),
        (mc['rho_planet'], r'$\rho_p$ [kg/m³]',    'seagreen'),
        (mc['Teq'],        r'$T_{eq}$ [K]',         'firebrick'),
        (mc['a_au'],       r'$a$ [AU]',             'purple'),
        (mc['H_atm'],      r'$H$ [km]',             'brown'),
    ]
    for ax, (data, label, color) in zip(axes.flat, params):
        clean = data[np.isfinite(data)]
        p16, p50, p84 = np.nanpercentile(clean, [16, 50, 84])
        ax.hist(clean, bins=60, color=color, alpha=0.75, edgecolor='none')
        ax.axvline(p50, color='black', linewidth=1.5, label='median')
        ax.axvline(p16, color='black', linewidth=1.0, linestyle='--')
        ax.axvline(p84, color='black', linewidth=1.0, linestyle='--')
        ax.set_xlabel(label, fontsize=11)
        ax.set_ylabel('N', fontsize=10)
        ax.set_title(f'{p50:.3g} $^{{+{p84-p50:.3g}}}_{{-{p50-p16:.3g}}}$', fontsize=10)
        ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig('posterior_distributions.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  → posterior_distributions.png")


# ============================================================================
# SAVE SUMMARY REPORT
# ============================================================================

def save_summary(mc, config, cls):
    """
    Write a complete characterization summary to characterization_summary.txt.

    Parameters
    ----------
    mc     : dict — Monte Carlo output from propagate_uncertainties()
    config : Config instance
    cls    : dict — classification output from classify_planet()
    """
    def wp(f, array, name, unit, dec=4):
        p16, p50, p84 = np.nanpercentile(array, [16, 50, 84])
        f.write(f"  {name}: {p50:.{dec}f} +{p84-p50:.{dec}f} -{p50-p16:.{dec}f} {unit}\n")

    def wmap(f, array, name, unit, dec=4):
        from scipy.stats import gaussian_kde
        clean = array[np.isfinite(array)]
        kde = gaussian_kde(clean)
        x = np.linspace(np.percentile(clean, 1), np.percentile(clean, 99), 1000)
        map_val = x[np.argmax(kde(x))]
        f.write(f"  {name}: {map_val:.{dec}f} {unit}  [MAP]\n")

    with open('characterization_summary.txt', 'w', encoding='utf-8') as f:
        f.write("="*70 + "\n")
        f.write(f"FULL CHARACTERIZATION — TOI {config.toi}\n")
        f.write("="*70 + "\n\n")
        f.write(f"TARGET        : TIC {config.tic_id}\n")
        f.write(f"Classification: {cls['overall']}\n\n")

        f.write("STELLAR PARAMETERS:\n")
        f.write(f"  Teff   = {config.teff} ± {config.teff_err} K\n")
        f.write(f"  M*     = {config.m_star} ± {config.m_star_err} M_sun\n")
        f.write(f"  R*     = {config.r_star} ± {config.r_star_err} R_sun\n")
        f.write(f"  log g  = {config.logg} ± {config.logg_err}\n")
        f.write(f"  [Fe/H] = {config.feh} ± {config.feh_err}\n")
        f.write(f"  mJ     = {config.mag_J} (2MASS, verified via astroquery)\n")

        f.write("\nORBITAL PARAMETERS:\n")
        wp(f,  mc['P_orbital'], "P",  "days", dec=7)
        wp(f,  mc['T0_bjd'],   "T0", "BJD",  dec=6)
        f.write(f"  e  = 0 (fixed, HJ tidal circularization)\n")
        wp(f,   mc['a_au'],     "a (median)", "AU")
        wmap(f, mc['a_au'],     "a (MAP)",    "AU")
        wp(f,   mc['b_impact'], "b",          "")

        f.write("\nPLANETARY PARAMETERS:\n")
        wp(f,   mc['Rp_Rjup'],   "Rp (median)", "R_Jup")
        wmap(f, mc['Rp_Rjup'],   "Rp (MAP)",    "R_Jup")
        wp(f,   mc['Rp_Rearth'], "Rp (median)", "R_Earth", dec=2)
        wmap(f, mc['Rp_Rearth'], "Rp (MAP)",    "R_Earth", dec=2)
        wp(f,   mc['Mp_Mjup'],   "Mp (median)", "M_Jup")
        wmap(f, mc['Mp_Mjup'],   "Mp (MAP)",    "M_Jup")
        wp(f,   mc['Mp_Mearth'], "Mp (median)", "M_Earth", dec=1)
        wmap(f, mc['Mp_Mearth'], "Mp (MAP)",    "M_Earth", dec=1)
        wp(f,   mc['rho_planet'],"rho_p (median)","kg/m³", dec=1)
        wmap(f, mc['rho_planet'],"rho_p (MAP)",   "kg/m³", dec=1)

        f.write("\nTEMPERATURE AND ATMOSPHERE:\n")
        wp(f,   mc['Teq'],        "T_eq (median)", "K",       dec=0)
        wmap(f, mc['Teq'],        "T_eq (MAP)",    "K",       dec=0)
        wp(f,   mc['H_atm'],      "H (median)",    "km",      dec=0)
        wmap(f, mc['H_atm'],      "H (MAP)",       "km",      dec=0)
        wp(f,   mc['F_incident'], "F_p",           "F_Earth", dec=0)

        Re      = float(np.nanmedian(mc['Rp_Rearth']))
        Me      = float(np.nanmedian(mc['Mp_Mearth']))
        Tq      = float(np.nanmedian(mc['Teq']))
        tsm       = calc_tsm(Re, Me,    config.r_star, Tq, config.mag_J)
        tsm_1mjup = calc_tsm(Re, 318.0, config.r_star, Tq, config.mag_J)
        f.write(f"\nATMOSPHERIC OBSERVABILITY:\n")
        f.write(f"  TSM (Chen-Kipping mass)  = {tsm:.1f}\n")
        f.write(f"  TSM (Mp = 1 MJup)        = {tsm_1mjup:.1f}\n")

        f.write("\nMETHODOLOGY AND REFERENCES:\n")
        f.write("  Bayesian fit     : Juliet + Dynesty\n")
        f.write("    Ref              : Espinoza et al. (2019), MNRAS 490, 2262\n")
        f.write("  Nested sampling  : Dynesty\n")
        f.write("    Ref              : Speagle (2020), MNRAS 493, 3132\n")
        f.write("  Limb darkening   : Quadratic law, Kipping (2013) parametrization\n")
        f.write("    Coefficients     : Claret (2017), A&A 600, A30\n")
        f.write("  Mass estimate    : Chen-Kipping (2017) M-R relation via mr_forecast\n")
        f.write("    Ref              : Chen & Kipping (2017), ApJ 834, 17\n")
        f.write("  MAP estimate     : Gaussian KDE (scipy.stats.gaussian_kde)\n")
        f.write("  Error propagation: Monte Carlo (10,000 samples)\n")
        f.write("  TSM formula      : Kempton et al. (2018), PASP 130, 114401\n")
    print("  → characterization_summary.txt")


# ============================================================================
# MAIN
# ============================================================================

def main():
    config = Config()

    print("="*70)
    print("EXOPLANET FULL CHARACTERIZATION PIPELINE")
    print(f"Target: TOI {config.toi} (TIC {config.tic_id})")
    print("="*70)

    # [1] Download TESS light curve
    print("\n[1/6] Downloading TESS light curve...")
    try:
        search = lk.search_lightcurve(f"TIC {config.tic_id}", mission="TESS", author="SPOC")
        if len(search) == 0:
            search = lk.search_lightcurve(f"TIC {config.tic_id}", mission="TESS")
        if len(search) == 0:
            raise RuntimeError(f"No data found for TIC {config.tic_id}")
        lc_raw = search.download_all().stitch().remove_nans()
        print(f"  → {len(lc_raw)} data points downloaded")
    except Exception as e:
        print(f"  ✗ Error: {e}"); raise

    # [2] Diagnostic: full light curve with transit markers
    print("\n[2/6] Full light curve diagnostic...")
    plot_full_lc_with_transits(lc_raw, config)

    # [3] Preprocessing and pre-fit verification
    print("\n[3/6] Preprocessing and pre-fit check...")
    epoch_btjd = config.epoch_bjd - 2457000.0

    if config.use_windows:
        # Keep only data within window_factor × duration around each transit
        phase = ((lc_raw.time.value - epoch_btjd) % config.period) / config.period
        phase = np.where(phase > 0.5, phase - 1, phase)
        dur_d = config.duration_hours / 24.0
        mask  = np.abs(phase) < (dur_d * config.window_factor / config.period)
        lc_windowed = lc_raw[mask]
    else:
        lc_windowed = lc_raw

    lc_final = lc_windowed.bin(time_bin_size=config.binning_min * u.min).remove_nans()
    times    = lc_final.time.value
    flux     = lc_final.flux.value
    flux_err = lc_final.flux_err.value
    print(f"  → {len(lc_final)} data points after preprocessing")

    plot_phased_check(times, flux, config)

    # [4] Bayesian fit (or load from cache)
    print("\n[4/6] Bayesian transit fit...")
    results = run_juliet_fit(times, flux, flux_err, config)

    # [5] Monte Carlo error propagation
    print("\n[5/6] Monte Carlo error propagation...")
    mc = propagate_uncertainties(results, config)

    # [6] Results and output files
    print("\n[6/6] Generating results...")
    cls, Rp, Mp, Teq, rho = print_results(mc, config)
    plot_population_comparison(Rp, Teq, config)
    plot_posterior_distributions(mc, config)
    save_summary(mc, config, cls)

    print("\n" + "="*70)
    print("✓✓✓ CHARACTERIZATION COMPLETE")
    print("="*70)
    print("\nOUTPUT FILES:")
    print("  📄 diagnostic_full_lc.png        ← Full light curve + transit markers")
    print("  📄 phased_check.png              ← Phase-folded binned curve (transit visible?)")
    print("  📄 corner_characterization.png   ← Fit posterior corner plot")
    print("  📄 posterior_distributions.png   ← Derived parameter posteriors")
    print("  📄 population_comparison.png     ← Hot Jupiter population context")
    print("  📄 characterization_summary.txt  ← Full numerical summary")
    print(f"  📁 {config.out_folder}/          ← Juliet posteriors (HDF5)")
    print("\n  ⚠️  PRIORITY CHECK: diagnostic_full_lc.png")
    print("      Red lines must coincide with flux dips!")
    print("="*70)


if __name__ == '__main__':
    main()