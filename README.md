# Ground Motion Plotter

A small Python tool that reads an earthquake acceleration record from a CSV file, plots acceleration, velocity and displacement, and reports the main intensity measures.

I built it as a beginner project while studying earthquake engineering. It is a first step towards comparing near-fault and far-fault ground motions.

![Near-fault example](outputs/sample_near_fault.png)

## What it calculates

| Measure | Meaning |
|---|---|
| **PGA** | Peak ground acceleration (g) |
| **PGV** | Peak ground velocity (cm/s) |
| **PGD** | Peak ground displacement (cm) |
| **Arias intensity** | Total shaking energy in the record (m/s) |
| **D5-95** | Significant duration: time between 5% and 95% of the energy |

## Quick start

```bash
git clone https://github.com/shahpravin4199/ground-motion-plotter.git
cd ground-motion-plotter
pip install -r requirements.txt

python make_sample_data.py                        # creates two practice records
python ground_motion.py data/sample_near_fault.csv data/sample_far_fault.csv
```

Plots are saved in the `outputs/` folder.

## Example output

| record | PGA (g) | PGV (cm/s) | PGD (cm) | Arias (m/s) | D5-95 (s) |
|---|---|---|---|---|---|
| sample_near_fault | 0.588 | 95.7 | 95.8 | 3.484 | 7.6 |
| sample_far_fault | 0.181 | 19.5 | 25.7 | 0.477 | 15.0 |

The near-fault record has a strong velocity pulse and short, intense shaking. The far-fault record is weaker and lasts longer.

## Using your own record

Make a CSV file with a header row and two columns, time and acceleration:

```
time_s,acc_g
0.00,0.0000
0.01,0.0012
...
```

Then run:

```bash
python ground_motion.py my_record.csv --units g
```

`--units` can be `g` (default), `m/s2` or `cm/s2`.

## How it works

1. **Load** the time and acceleration columns.
2. **Baseline correct**: remove the mean and filter out very slow drift (below 0.1 Hz).
3. **Integrate** acceleration once for velocity and again for displacement (trapezoidal rule).
4. **Measure** the peaks, the Arias intensity `Ia = π/(2g) ∫ a² dt`, and the 5-95% duration from the cumulative energy curve.
5. **Plot** the three time histories, with the 5-95% window shaded.

## Limitations

- The sample records are **synthetic** (made up for practice), not real earthquakes.
- The baseline correction is simple. Real records need more careful processing, and displacement is the measure most sensitive to it.
- Only one horizontal component is analysed at a time.

## Ideas for next steps

- [ ] Add a 5%-damped response spectrum (using the Newmark method)
- [ ] Compare a near-fault and a far-fault record on the same plot
- [ ] Detect the velocity pulse and estimate its period
- [ ] Read real PEER NGA-West2 files directly

## Where to get real records

- PEER NGA-West2 ground motion database
- COSMOS Virtual Data Center

## Built with

Python, NumPy, Matplotlib.
