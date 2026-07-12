#!/usr/bin/env python3
"""
Astrometrica Report Parser & Publication-Grade Plotting Pipeline
Author: Aditya Rajiv Ratnam
"""

import os
import re
import numpy as np
import matplotlib.pyplot as plt

# Configuration
# NOTE: must match the designation used in MPCReport.txt
# / PhotReport.txt / Astrometrica.log. 
TARGET_ID = "AST0001"

FILE_PATH_MPC = "MPCReport.txt"
FILE_PATH_PHOT = "PhotReport.txt"
FILE_PATH_LOG = "Astrometrica.log"


def parse_all_logs(target_id):
    """
    Parses MPCReport.txt (astrometry), PhotReport.txt (magnitude/SNR), and
    Astrometrica.log (Flux, FWHM, Fit RMS, dRA/dDec residuals, plus per-image
    calibration quality: Zero Point and photometric dmag scatter).

    Astrometrica.log target-block matching is done by day_fraction rather than
    line order, because manually-added detections are logged in the order you
    added them (not necessarily chronological order).
    """
    observations = []

    # 1. Parse MPCReport.txt -> RA/Dec/day_fraction
    if os.path.exists(FILE_PATH_MPC):
        with open(FILE_PATH_MPC, 'r') as f:
            for line in f:
                if len(line) >= 60 and target_id in line:
                    ra_h = int(line[32:34])
                    ra_m = int(line[35:37])
                    ra_s = float(line[38:44])

                    dec_d_raw = line[44:47].strip()
                    dec_m = int(line[47:50])
                    dec_s = float(line[50:56])

                    ra_deg = (ra_h + ra_m / 60.0 + ra_s / 3600.0) * 15.0
                    sign = -1.0 if '-' in dec_d_raw else 1.0
                    dec_deg = sign * (abs(int(dec_d_raw)) + dec_m / 60.0 + dec_s / 3600.0)
                    day_fraction = float(line[15:32].split()[-1])

                    observations.append({
                        "day_fraction": day_fraction,
                        "ra_deg": ra_deg,
                        "dec_deg": dec_deg,
                        "snr": 0.0,
                        "mag": 0.0,
                        "fwhm": None,
                        "flux": None,
                        "fit_rms": None,
                        "dra_err": None,
                        "ddec_err": None,
                        "zero_point": None,
                        "dmag_scatter": None,
                        "ref_stars": None,
                    })

    observations.sort(key=lambda x: x["day_fraction"])

    # 2. Parse PhotReport.txt -> mag/SNR, matched by sorting on JD (each file's
    #    own time field), rather than assuming matching row order.
    if os.path.exists(FILE_PATH_PHOT):
        with open(FILE_PATH_PHOT, 'r') as f:
            phot_rows = [l.split() for l in f if target_id in l]
        phot_sorted = sorted(phot_rows, key=lambda r: float(r[0]))
        for obs, row in zip(observations, phot_sorted):
            try:
                obs["mag"] = float(row[1])
                obs["snr"] = float(row[4])
            except (ValueError, IndexError):
                pass

    if not os.path.exists(FILE_PATH_LOG):
        return observations

    with open(FILE_PATH_LOG, 'r', encoding='latin-1') as f:
        content = f.read()

    # 3. Parse per-target "Position added" blocks -> Flux, FWHM, Fit RMS,
    #    dRA/dDec residual errors, and the source image filename.
    block_re = re.compile(
        r'Position added [a-z]+ from file ([^\s:]+):\r?\n'   # filename
        r'([^\n]*)\n'                                        # data row
        r'([^\n]*)\n'                                        # residual/error row
        r'(\s*' + re.escape(target_id) + r'[^\n]*)\n'        # embedded designation line
    )

    matched, unmatched = 0, 0
    matched_filenames = []
    for fname, data_line, err_line, desig_line in block_re.findall(content):
        fields = data_line.split()
        if len(fields) < 4:
            continue
        try:
            # Trailing four columns are always Flux, FWHM, Peak(SNR), Fit RMS,
            # regardless of whether dRA/dDec/dG columns are present earlier.
            flux_val = float(fields[-4])
            fwhm_val = float(fields[-3])
            fit_rms_val = float(fields[-1])
        except ValueError:
            continue

        errs = re.findall(r'±([0-9.]+)', err_line)
        dra_err = float(errs[0]) if len(errs) >= 1 else None
        ddec_err = float(errs[1]) if len(errs) >= 2 else None

        day_fraction = float(desig_line[15:32].split()[-1])

        best_i = min(range(len(observations)),
                     key=lambda i: abs(observations[i]["day_fraction"] - day_fraction))
        if abs(observations[best_i]["day_fraction"] - day_fraction) < 1e-6:
            observations[best_i]["flux"] = flux_val
            observations[best_i]["fwhm"] = fwhm_val
            observations[best_i]["fit_rms"] = fit_rms_val
            observations[best_i]["dra_err"] = dra_err
            observations[best_i]["ddec_err"] = ddec_err
            matched_filenames.append((best_i, fname))
            matched += 1
        else:
            unmatched += 1

    if unmatched:
        print(f"[WARN] {unmatched} Astrometrica.log block(s) could not be matched "
              f"to an MPCReport observation by day_fraction.")
    print(f"[INFO] Matched Flux/FWHM/residuals for {matched} of {len(observations)} observations from Astrometrica.log.")

    # 4. Parse per-image photometric calibration quality: Zero Point and the
    #    dmag scatter of the reference-star solution. This describes how good
    #    the *frame's* overall calibration was, independent of the target's
    #    own detection.
    photo_block_re = re.compile(
        r'Photometry of Image \d+ \(([^)]+)\):\s*\r?\n'
        r'\s*(\d+) of \d+ Reference Stars used: dmag = ([0-9.]+)mag\s*\r?\n'
        r'\s*Zero Point: ([0-9.]+)mag'
    )
    frame_calib = {}
    for fname, n_stars, dmag, zp in photo_block_re.findall(content):
        frame_calib.setdefault(fname, {"dmag": float(dmag), "zp": float(zp), "n_stars": int(n_stars)})

    for obs_i, fname in matched_filenames:
        calib = frame_calib.get(fname)
        if calib:
            observations[obs_i]["zero_point"] = calib["zp"]
            observations[obs_i]["dmag_scatter"] = calib["dmag"]
            observations[obs_i]["ref_stars"] = calib["n_stars"]

    return observations


def compute_angular_velocity(data):
    """
    Pairwise angular speed (arcsec/hour) between each pair of consecutive
    frames, following the standard angular-distance convention (Richmond,
    2025). Returns (midpoint_time_hours, omega_arcsec_per_hr) arrays, one
    value per consecutive pair (n-1 values for n frames).
    """
    jd = np.array([d["day_fraction"] for d in data])
    ra = np.array([d["ra_deg"] for d in data])
    dec = np.array([d["dec_deg"] for d in data])
    time_hours = (jd - jd[0]) * 24.0

    mid_times, omegas = [], []
    for i in range(len(data) - 1):
        dalpha = ra[i + 1] - ra[i]
        ddelta = dec[i + 1] - dec[i]
        cosd = np.cos(np.radians((dec[i] + dec[i + 1]) / 2.0))
        theta_arcsec = np.hypot(dalpha * cosd, ddelta) * 3600.0
        dt_hr = (jd[i + 1] - jd[i]) * 24.0
        omegas.append(theta_arcsec / dt_hr)
        mid_times.append((time_hours[i] + time_hours[i + 1]) / 2.0)
    return np.array(mid_times), np.array(omegas)


def generate_plots(data):
    if not data:
        print("Data empty. Generation halted.")
        return

    jd = np.array([d["day_fraction"] for d in data])
    ra = np.array([d["ra_deg"] for d in data])
    dec = np.array([d["dec_deg"] for d in data])
    mag = np.array([d["mag"] for d in data])
    snr = np.array([d["snr"] for d in data])

    time_hours = (jd - jd[0]) * 24.0
    c_navy, c_slate, c_red = '#1B365D', '#5C768D', '#D62728'

    cos_dec = np.cos(np.radians(dec[0]))
    delta_ra_arcsec = (ra - ra[0]) * 3600.0 * cos_dec
    delta_dec_arcsec = (dec - dec[0]) * 3600.0

    # ----------------------------------------------------
    # FIGURE 1: Apparent Sky Motion and Linear Regression Model
    # ----------------------------------------------------
    fig1, ax1 = plt.subplots(figsize=(7, 6))
    slope, intercept = np.polyfit(delta_ra_arcsec, delta_dec_arcsec, 1)
    r_squared = np.corrcoef(delta_ra_arcsec, delta_dec_arcsec)[0, 1] ** 2
    ax1.plot(delta_ra_arcsec, delta_dec_arcsec, 'ko', markersize=6, label='Data points')
    ax1.plot(delta_ra_arcsec, slope * delta_ra_arcsec + intercept, 'r-', linewidth=1.5, label='Linear fit')
    equation_str = f"y = {slope:.4f}x + ({intercept:.4f})\n$R^2$ = {r_squared:.4f}"
    ax1.text(0.05, 0.90, equation_str, transform=ax1.transAxes, fontsize=10,
              bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="gray", alpha=0.8))
    ax1.set_title("Figure 1: Apparent Sky Motion and Linear Regression Model", fontsize=11, fontweight='bold', color=c_navy)
    ax1.set_xlabel("Right Ascension (arcsec)")
    ax1.set_ylabel("Declination (arcsec)")
    ax1.grid(True, linestyle=':', alpha=0.6)
    ax1.invert_xaxis()
    ax1.legend(loc='lower right')
    plt.tight_layout()
    fig1.savefig("Figure 1 Apparent Sky Motion and Linear Regression Model.png", dpi=300)
    plt.close(fig1)
    print(" -> 'Figure 1 Apparent Sky Motion and Linear Regression Model.png' generated successfully.")

    ra_rate = np.polyfit(time_hours, delta_ra_arcsec, 1)[0]
    dec_rate = np.polyfit(time_hours, delta_dec_arcsec, 1)[0]
    rate_arcsec_hr = float(np.hypot(ra_rate, dec_rate))
    pa_deg = float(np.degrees(np.arctan2(ra_rate, dec_rate)) % 360)
    print(f"[INFO] Derived rate of motion: ~{rate_arcsec_hr:.2f} arcsec/hr, position angle ~{pa_deg:.1f} deg")

    # ----------------------------------------------------
    # FIGURE 2: Photometric, Image Quality, Astrometric, and Kinematic Metrics
    #           3x2 grid grouped and color-coded by category:
    #           Row 1 (Photometric, blue):   Magnitude, Flux
    #           Row 2 (Image Quality, green): SNR, FWHM
    #           Row 3 (Astrometric/red, Kinematic/purple): Fit RMS, Angular Velocity
    # ----------------------------------------------------
    have_fwhm_flux = all(d["fwhm"] is not None and d["flux"] is not None for d in data)
    have_fit_rms = all(d["fit_rms"] is not None for d in data)

    c_blue, c_green, c_fitrms, c_purple = '#4C72B0', '#55A868', '#C44E52', '#8172B2'

    fig2, axs2 = plt.subplots(3, 2, figsize=(10, 14))
    fig2.suptitle("Figure 2: Photometric, Image Quality, Astrometric, and Kinematic Metrics",
                  fontsize=13, fontweight='bold', color=c_navy, y=0.98)

    axs2[0, 0].plot(time_hours, mag, 'o-', color=c_blue, linewidth=1.5)
    axs2[0, 0].set_xlabel("Elapsed Time (Hours)", fontweight='bold')
    axs2[0, 0].set_ylabel("Magnitude (mag)", fontweight='bold')
    axs2[0, 0].set_ylim(max(mag) + 0.2, min(mag) - 0.2)
    axs2[0, 0].grid(True, linestyle=':', alpha=0.6)
    axs2[0, 0].set_title("Magnitude vs Time")

    if have_fwhm_flux:
        flux = np.array([d["flux"] for d in data])
        axs2[0, 1].plot(time_hours, flux, 'd-', color=c_blue, linewidth=1.5)
        axs2[0, 1].set_xlabel("Elapsed Time (Hours)", fontweight='bold')
        axs2[0, 1].set_ylabel("Flux (ADU)", fontweight='bold')
        axs2[0, 1].grid(True, linestyle=':', alpha=0.6)
        axs2[0, 1].set_title("Flux vs Time")
    else:
        axs2[0, 1].text(0.5, 0.5, "No Flux data parsed\nfrom Astrometrica.log",
                         transform=axs2[0, 1].transAxes, ha='center', va='center', color='gray', style='italic')
        axs2[0, 1].set_title("Flux vs Time")
        axs2[0, 1].set_xlabel("Elapsed Time (Hours)", fontweight='bold')

    axs2[1, 0].plot(time_hours, snr, 's-', color=c_green, linewidth=1.5)
    axs2[1, 0].set_xlabel("Elapsed Time (Hours)", fontweight='bold')
    axs2[1, 0].set_ylabel("SNR", fontweight='bold')
    axs2[1, 0].grid(True, linestyle=':', alpha=0.6)
    axs2[1, 0].set_title("SNR vs Time")

    if have_fwhm_flux:
        fwhm = np.array([d["fwhm"] for d in data])
        axs2[1, 1].plot(time_hours, fwhm, '^-', color=c_green, linewidth=1.5)
        axs2[1, 1].set_xlabel("Elapsed Time (Hours)", fontweight='bold')
        axs2[1, 1].set_ylabel("FWHM (arcsec)", fontweight='bold')
        axs2[1, 1].grid(True, linestyle=':', alpha=0.6)
        axs2[1, 1].set_title("FWHM vs Time")
    else:
        axs2[1, 1].text(0.5, 0.5, "No FWHM data parsed\nfrom Astrometrica.log",
                         transform=axs2[1, 1].transAxes, ha='center', va='center', color='gray', style='italic')
        axs2[1, 1].set_title("FWHM vs Time")
        axs2[1, 1].set_xlabel("Elapsed Time (Hours)", fontweight='bold')

    if have_fit_rms:
        fit_rms = np.array([d["fit_rms"] for d in data])
        axs2[2, 0].plot(time_hours, fit_rms, 'v-', color=c_fitrms, linewidth=1.5)
        axs2[2, 0].set_xlabel("Elapsed Time (Hours)", fontweight='bold')
        axs2[2, 0].set_ylabel("Fit RMS (arcsec)", fontweight='bold')
        axs2[2, 0].grid(True, linestyle=':', alpha=0.6)
        axs2[2, 0].set_title("Fit RMS vs Time")
    else:
        axs2[2, 0].text(0.5, 0.5, "No Fit RMS data parsed\nfrom Astrometrica.log",
                         transform=axs2[2, 0].transAxes, ha='center', va='center', color='gray', style='italic')
        axs2[2, 0].set_title("Fit RMS vs Time")
        axs2[2, 0].set_xlabel("Elapsed Time (Hours)", fontweight='bold')

    omega_times, omega_vals = compute_angular_velocity(data)
    axs2[2, 1].plot(omega_times, omega_vals, 'h-', color=c_purple, linewidth=1.5)
    axs2[2, 1].set_xlabel("Elapsed Time (Hours)", fontweight='bold')
    axs2[2, 1].set_ylabel("Angular Velocity (arcsec/hr)", fontweight='bold')
    axs2[2, 1].grid(True, linestyle=':', alpha=0.6)
    axs2[2, 1].set_title("Angular Velocity vs Time")

    fig2.subplots_adjust(top=0.90, bottom=0.06, hspace=0.55, wspace=0.32, left=0.10, right=0.97)

    # Gray category labels sit close to the top of the row they describe
    # (not centered in the gap above them), so each label reads as
    # belonging to the row directly below it rather than the row above.
    label_kwargs = dict(ha='center', fontsize=11, color='gray', fontweight='bold')
    pos0 = (axs2[0, 0].get_position(), axs2[0, 1].get_position())
    pos1 = (axs2[1, 0].get_position(), axs2[1, 1].get_position())
    pos2 = (axs2[2, 0].get_position(), axs2[2, 1].get_position())

    def label_above(pos_top_of_row, y_bottom_of_row_above, x_center, text):
        # Weighted 80% of the way down from the row-above's bottom toward
        # this row's top, so the label sits just above this row's title.
        y = y_bottom_of_row_above + 0.70 * (pos_top_of_row - y_bottom_of_row_above)
        fig2.text(x_center, y, text, **label_kwargs)

    label_above(pos0[0].y1, 0.98, (pos0[0].x0 + pos0[1].x1) / 2, "Photometric")
    label_above(pos1[0].y1, pos0[0].y0, (pos1[0].x0 + pos1[1].x1) / 2, "Image Quality")
    label_above(pos2[0].y1, pos1[0].y0, (pos2[0].x0 + pos2[0].x1) / 2, "Astrometric")
    label_above(pos2[1].y1, pos1[0].y0, (pos2[1].x0 + pos2[1].x1) / 2, "Kinematic")

    fig2.savefig("Figure 2 Photometric Image Quality Astrometric and Kinematic Metrics.png", dpi=300)
    plt.close(fig2)
    print(" -> 'Figure 2 Photometric Image Quality Astrometric and Kinematic Metrics.png' generated successfully.")

    # ----------------------------------------------------
    # FIGURE 3: Astrometric Plate Reduction Residuals
    # ----------------------------------------------------
    have_residuals = all(d["dra_err"] is not None and d["ddec_err"] is not None for d in data)
    if have_residuals:
        ra_err = np.array([d["dra_err"] for d in data])
        dec_err = np.array([d["ddec_err"] for d in data])

        fig3, ax3 = plt.subplots(figsize=(7, 5))
        ax3.plot(time_hours, ra_err, 's--', color=c_red, label=r'dRA Residual ($\sigma$)')
        ax3.plot(time_hours, dec_err, '^--', color=c_slate, label=r'dDec Residual ($\sigma$)')
        ax3.axhline(0.40, color='black', linestyle='--', linewidth=1.0, alpha=0.7,
                    label='Typical Astrometric Acceptance Threshold (0.40\u2033)')
        ax3.set_title("Figure 3: Astrometric Plate Solution Residuals",
                      fontsize=12, fontweight='bold', color=c_navy)
        ax3.set_xlabel("Elapsed Time (Hours)")
        ax3.set_ylabel("Residual (arcsec)")
        ax3.set_ylim(0, 0.46)
        ax3.grid(True, linestyle=':', alpha=0.6)
        ax3.legend(fontsize=8, loc='center right')
        plt.tight_layout()
        fig3.savefig("Figure 3 Astrometric Plate Solution Residuals.png", dpi=300)
        plt.close(fig3)
        print(" -> 'Figure 3 Astrometric Plate Solution Residuals.png' generated successfully.")
    else:
        print("[WARN] Skipped Figure 3 - residual errors not fully parsed from Astrometrica.log.")

    # ----------------------------------------------------
    # FIGURE 4: Frame Calibration Quality (Zero Point + photometric dmag scatter)
    #           New: explains WHY SNR/mag behave the way they do, independent
    #           of the target itself. Real per-image data from the log.
    # ----------------------------------------------------
    have_calib = all(d["zero_point"] is not None and d["dmag_scatter"] is not None for d in data)
    if have_calib:
        zp = np.array([d["zero_point"] for d in data])
        dmag = np.array([d["dmag_scatter"] for d in data])

        fig4, axs4 = plt.subplots(1, 2, figsize=(11, 5))
        fig4.suptitle("Figure 4: Plate Calibration Stability During the Observing Session",
                      fontsize=12, fontweight='bold', color=c_navy)

        axs4[0].plot(time_hours, zp, 'o-', color=c_navy, linewidth=1.5)
        axs4[0].set_title("Photometric Zero Point vs Time", fontsize=11, fontweight='bold')
        axs4[0].set_xlabel("Elapsed Time (Hours)")
        axs4[0].set_ylabel("Zero Point (mag)")
        axs4[0].grid(True, linestyle=':', alpha=0.6)

        axs4[1].plot(time_hours, dmag, 's-', color=c_red, linewidth=1.5)
        axs4[1].set_title("Reference-Star Calibration Scatter vs Time", fontsize=11, fontweight='bold')
        axs4[1].set_xlabel("Elapsed Time (Hours)")
        axs4[1].set_ylabel("dmag (mag)")
        axs4[1].grid(True, linestyle=':', alpha=0.6)

        plt.tight_layout()
        fig4.savefig("Figure 4 Plate Calibration Stability.png", dpi=300)
        plt.close(fig4)
        print(" -> 'Figure 4 Plate Calibration Stability.png' generated successfully.")
    else:
        print("[WARN] Skipped Figure 4 - Zero Point / dmag not fully parsed from Astrometrica.log.")

    print("\n[SUCCESS] All figures generated.")


def compute_position_angle(data):
    """
    Pairwise position angle (degrees, east of north) between each pair of
    consecutive frames, matching the cardinality of compute_angular_velocity.
    """
    ra = np.array([d["ra_deg"] for d in data])
    dec = np.array([d["dec_deg"] for d in data])
    pas = []
    for i in range(len(data) - 1):
        dalpha = ra[i + 1] - ra[i]
        ddelta = dec[i + 1] - dec[i]
        cosd = np.cos(np.radians((dec[i] + dec[i + 1]) / 2.0))
        pas.append(np.degrees(np.arctan2(dalpha * cosd, ddelta)) % 360)
    return np.array(pas)


def compute_summary_statistics(data):
    """
    Mean / SD (sample, ddof=1) / min / max for every metric reported in
    Table 2 of the report, including the pairwise-derived Angular Velocity
    and Position Angle (n=3 for n=4 frames) and per-frame Reference Stars Used.
    """
    omega_times, omega_vals = compute_angular_velocity(data)
    pa_vals = compute_position_angle(data)

    fields = {
        "Magnitude (mag)": np.array([d["mag"] for d in data]),
        "Flux (ADU)": np.array([d["flux"] for d in data]),
        "SNR": np.array([d["snr"] for d in data]),
        "FWHM (arcsec)": np.array([d["fwhm"] for d in data]),
        "Fit RMS (arcsec)": np.array([d["fit_rms"] for d in data]),
        "dRA Residual (arcsec)": np.array([d["dra_err"] for d in data]),
        "dDec Residual (arcsec)": np.array([d["ddec_err"] for d in data]),
        "Zero Point (mag)": np.array([d["zero_point"] for d in data]),
        "Reference-Star Calibration Scatter (mag)": np.array([d["dmag_scatter"] for d in data]),
        "Angular Velocity (arcsec/hr)": omega_vals,
        "Position Angle (deg)": pa_vals,
        "Reference Stars Used": np.array([d["ref_stars"] for d in data], dtype=float),
    }
    stats = {}
    for name, arr in fields.items():
        stats[name] = {
            "mean": float(np.mean(arr)),
            "sd": float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0,
            "min": float(np.min(arr)),
            "max": float(np.max(arr)),
            "n": len(arr),
        }
    return stats


if __name__ == "__main__":
    print(f"Running report execution parsing engine for target profile: {TARGET_ID}")
    parsed_data = parse_all_logs(TARGET_ID)
    if len(parsed_data) >= 2:
        stats = compute_summary_statistics(parsed_data)
        print("\n[INFO] Summary statistics (mean / SD / min / max):")
        for name, s in stats.items():
            print(f"  {name:42s} n={s['n']}  mean={s['mean']:.4f}  sd={s['sd']:.4f}  min={s['min']:.4f}  max={s['max']:.4f}")
        print()
        generate_plots(parsed_data)
    else:
        print(f"[ERROR] Fewer than 2 observations found for TARGET_ID='{TARGET_ID}'. "
              f"Check that TARGET_ID matches the designation in MPCReport.txt.")
