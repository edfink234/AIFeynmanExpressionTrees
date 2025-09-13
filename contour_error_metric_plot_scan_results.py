# Read the CSV, build a 3D surface, and save high-res PNG and vector PDF

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the grid scan CSV
path = "/Users/edwardfinkelstein/Downloads/grid_scan_results.csv"
df = pd.read_csv(path, header=None)

# Detect header structure:
# Row 0: ["b \\ Omega", Om1, Om2, ...]
# Col 0 (from row 1 down): b values
# Numeric grid: rows 1.., cols 1..
# Coerce to numeric where possible
omega_vals = pd.to_numeric(df.iloc[0, 1:], errors="coerce").to_numpy(dtype=float)
b_vals = pd.to_numeric(df.iloc[1:, 0], errors="coerce").to_numpy(dtype=float)
Z = df.iloc[1:, 1:].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)

# --- Find the minimum D and its corresponding (b, Omega) ---------
min_idx = np.unravel_index(np.nanargmin(Z), Z.shape)
min_D = Z[min_idx]
min_b = b_vals[min_idx[0]]
min_omega = omega_vals[min_idx[1]]

print(f"Minimum 𝒟 = {min_D:.6g} at b = {min_b}, Ω = {min_omega}")
# ---------------------------------------------------------------


# Build meshgrid
Om_grid, b_grid = np.meshgrid(omega_vals, b_vals)

# Mask invalid entries for plotting
Z_masked = np.ma.array(Z, mask=~np.isfinite(Z))

# Generate a 2D contourf heatmap of log10(D) for publication-style figure

fig, ax = plt.subplots(figsize=(7, 5), dpi=200)

# Mask invalid values, ensure positive
Z_pos = np.where(Z <= 0, np.nan, Z)
logZ = np.log10(Z_pos)

# Contourf heatmap
c = ax.contourf(Om_grid, b_grid, logZ, levels=30, cmap="viridis")
cb = fig.colorbar(c, ax=ax)
cb.set_label(r'$\log_{10}\mathcal{D}$')

# Labels
ax.set_xlabel(r'$\Omega$')
ax.set_ylabel(r'$b$')

# Optional: contour lines
contours = ax.contour(Om_grid, b_grid, logZ, levels=10, colors='k', linewidths=0.5, alpha=0.5)
ax.clabel(contours, inline=True, fontsize=6, fmt="%.1f")

fig.tight_layout()

pdf_contour_path = "grid_scan_contour_log.pdf"
fig.savefig(pdf_contour_path, bbox_inches="tight")

from os import system
system(f"open {pdf_contour_path}")

