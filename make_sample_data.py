"""Create two SYNTHETIC ground motion records so you can try the plotter
without downloading anything.

These are made-up signals for practice, not real earthquake recordings.
For real work, download records from PEER NGA-West2, COSMOS or similar.

Run:  python make_sample_data.py
"""
import numpy as np

DT = 0.01  # time step in seconds (100 samples per second)
G = 9.81   # m/s^2


def envelope(t, rise, peak_time, decay):
    """A smooth 'shaking builds up, then dies away' shape."""
    x = t / peak_time
    env = x ** rise * np.exp(rise * (1 - x))
    return env * np.exp(-np.maximum(t - peak_time, 0) / decay)


def smooth_noise(n, rng, width=5):
    """Random noise, lightly smoothed so it looks more like real shaking."""
    noise = rng.standard_normal(n)
    kernel = np.ones(width) / width
    return np.convolve(noise, kernel, mode="same")


def near_fault(rng):
    """Strong shaking with one big velocity pulse (period about 3 s)."""
    t = np.arange(0, 30, DT)
    tp, t0 = 3.0, 8.0  # pulse period and pulse centre time
    wp = 2 * np.pi / tp
    # velocity pulse (0.9 m/s): a cosine wave inside a bell-shaped (Gaussian) window
    vel = 0.9 * np.exp(-((t - t0) / (0.6 * tp)) ** 2) * np.cos(wp * (t - t0))
    pulse_acc = np.gradient(vel, DT)                 # acceleration = d(velocity)/dt
    background = smooth_noise(len(t), rng) * envelope(t, 2.5, 8, 6)
    background = background / np.max(np.abs(background)) * 0.5 * G   # about 0.5 g of shaking
    return t, (pulse_acc + background) / G


def far_fault(rng):
    """Longer, gentler shaking with no clear pulse."""
    t = np.arange(0, 60, DT)
    acc_g = smooth_noise(len(t), rng, width=9) * envelope(t, 3, 14, 14)
    return t, acc_g / np.max(np.abs(acc_g)) * 0.18   # scale so PGA = 0.18 g


def save(name, t, acc_g):
    data = np.column_stack([t, acc_g])
    path = f"data/{name}.csv"
    np.savetxt(path, data, delimiter=",", header="time_s,acc_g", comments="", fmt="%.5f")
    print("saved", path)


if __name__ == "__main__":
    rng = np.random.default_rng(42)  # fixed seed -> same files every time
    save("sample_near_fault", *near_fault(rng))
    save("sample_far_fault", *far_fault(rng))
