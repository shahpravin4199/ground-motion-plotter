"""Ground motion plotter.

Reads an acceleration record from a CSV file, calculates basic intensity
measures (PGA, PGV, PGD, Arias intensity, significant duration) and saves
a plot of acceleration, velocity and displacement.

Usage:
    python ground_motion.py data/sample_near_fault.csv
    python ground_motion.py my_record.csv --units m/s2
    python ground_motion.py data/a.csv data/b.csv      # several records at once

CSV format: a header row, then two columns: time (s), acceleration.
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

G = 9.81  # m/s^2


def load_record(path, units="g"):
    """Return time (s) and acceleration (m/s^2) from a two-column CSV."""
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    t, acc = data[:, 0], data[:, 1]
    if units == "g":
        acc = acc * G
    elif units == "cm/s2":
        acc = acc / 100
    elif units != "m/s2":
        raise ValueError("units must be g, m/s2 or cm/s2")
    return t, acc


def integrate(y, dt):
    """Cumulative trapezoidal integral of y, starting from zero."""
    out = np.zeros_like(y)
    out[1:] = np.cumsum((y[1:] + y[:-1]) * dt / 2)
    return out


def baseline_correct(acc, dt, f_cut=0.1, order=4):
    """Simple baseline correction: remove the mean, then gently remove very
    slow drift (below f_cut Hz) with a smooth high-pass filter.

    Without this, tiny errors in the record grow into huge, fake
    displacements when we integrate twice. Real processing (for example
    PEER's) is more careful; this is a good first version.
    """
    acc = acc - np.mean(acc)
    spectrum = np.fft.rfft(acc)
    freqs = np.fft.rfftfreq(len(acc), dt)
    gain = 1 / np.sqrt(1 + (f_cut / np.maximum(freqs, 1e-9)) ** (2 * order))
    return np.fft.irfft(spectrum * gain, n=len(acc))


def analyse(t, acc):
    """Calculate velocity, displacement and the key intensity measures."""
    dt = t[1] - t[0]
    acc = baseline_correct(acc, dt)
    vel = integrate(acc, dt)
    disp = integrate(vel, dt)

    # Arias intensity and Husid curve (cumulative energy as a fraction)
    energy = integrate(acc**2, dt)
    arias = np.pi / (2 * G) * energy[-1]
    husid = energy / energy[-1]
    t5, t95 = np.interp([0.05, 0.95], husid, t)

    i_pga = np.argmax(np.abs(acc))
    return {
        "acc": acc, "vel": vel, "disp": disp,
        "pga_g": np.abs(acc[i_pga]) / G,
        "t_pga": t[i_pga],
        "pgv_cm": np.max(np.abs(vel)) * 100,
        "pgd_cm": np.max(np.abs(disp)) * 100,
        "arias": arias,
        "d595": t95 - t5,
        "t5": t5, "t95": t95,
        "dt": dt, "duration": t[-1] - t[0],
    }


def plot_record(t, r, title, out_path):
    """Save a three-panel plot: acceleration, velocity, displacement."""
    fig, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True)

    panels = [
        (r["acc"] / G, "Acceleration (g)", "#c2410c"),
        (r["vel"] * 100, "Velocity (cm/s)", "#23313f"),
        (r["disp"] * 100, "Displacement (cm)", "#23313f"),
    ]
    for ax, (y, label, colour) in zip(axes, panels):
        ax.plot(t, y, color=colour, linewidth=0.9)
        ax.axhline(0, color="grey", linewidth=0.5)
        ax.axvspan(r["t5"], r["t95"], color="#f26b21", alpha=0.08)  # 5-95 % energy window
        ax.set_ylabel(label)
        ax.grid(alpha=0.25)

    # mark the peak ground acceleration
    axes[0].plot(r["t_pga"], np.interp(r["t_pga"], t, r["acc"] / G), "o", color="black")
    axes[0].annotate(f"PGA = {r['pga_g']:.2f} g", (r["t_pga"], np.interp(r["t_pga"], t, r["acc"] / G)),
                     textcoords="offset points", xytext=(10, 8))

    axes[-1].set_xlabel("Time (s)")
    summary = (f"PGV {r['pgv_cm']:.1f} cm/s   PGD {r['pgd_cm']:.1f} cm   "
               f"Arias {r['arias']:.2f} m/s   D5-95 {r['d595']:.1f} s")
    fig.suptitle(f"{title}\n{summary}", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Plot a ground motion record and report intensity measures.")
    parser.add_argument("files", nargs="+", help="CSV file(s): time, acceleration")
    parser.add_argument("--units", default="g", choices=["g", "m/s2", "cm/s2"],
                        help="units of the acceleration column (default: g)")
    parser.add_argument("--outdir", default="outputs", help="folder for the saved plots")
    args = parser.parse_args()

    Path(args.outdir).mkdir(exist_ok=True)
    print(f"{'record':<22}{'PGA(g)':>8}{'PGV(cm/s)':>11}{'PGD(cm)':>9}{'Arias(m/s)':>12}{'D5-95(s)':>10}")

    for file in args.files:
        t, acc = load_record(file, args.units)
        r = analyse(t, acc)
        name = Path(file).stem
        plot_record(t, r, name, Path(args.outdir) / f"{name}.png")
        print(f"{name:<22}{r['pga_g']:>8.3f}{r['pgv_cm']:>11.1f}{r['pgd_cm']:>9.1f}"
              f"{r['arias']:>12.3f}{r['d595']:>10.1f}")

    print(f"\nPlots saved in the '{args.outdir}' folder.")


if __name__ == "__main__":
    main()
