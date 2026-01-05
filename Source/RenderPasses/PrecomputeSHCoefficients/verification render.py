import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 1. Load Data
df = pd.read_csv("verification_data.csv")
resolution = int(np.sqrt(len(df)))

# 2. Reshape
X = df['WorldX'].values.reshape(resolution, resolution)
Z = df['WorldZ'].values.reshape(resolution, resolution)
grad_analytic = df['AnalyticGrad'].values.reshape(resolution, resolution)
grad_fd       = df['FDGrad'].values.reshape(resolution, resolution)
hess_analytic = df['AnalyticHess'].values.reshape(resolution, resolution)
hess_fd       = df['FDHess'].values.reshape(resolution, resolution)

# 3. Calculate 10x Scaled Error
grad_diff_10x = np.abs(grad_analytic - grad_fd) * 10
hess_diff_10x = np.abs(hess_analytic - hess_fd) * 10

# 4. Plotting
fig, axes = plt.subplots(2, 3, figsize=(12, 8)) 

def plot_smooth(ax, data, title, vmin=None, vmax=None, cmap='jet'):
    im = ax.imshow(data, extent=[X.min(), X.max(), Z.min(), Z.max()], 
                   origin='lower', cmap=cmap, vmin=vmin, vmax=vmax,
                   interpolation='bicubic')
    ax.set_title(title)
    ax.axis('off')
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

# --- Row 1: Gradients ---
g_min = min(grad_analytic.min(), grad_fd.min())
g_max = max(grad_analytic.max(), grad_fd.max())
g_range = g_max - g_min

plot_smooth(axes[0, 0], grad_analytic, "Analytic Gradient", vmin=g_min, vmax=g_max)
plot_smooth(axes[0, 1], grad_fd, "Finite Diff Gradient", vmin=g_min, vmax=g_max)
# Plot 10x Error against the Signal Range
plot_smooth(axes[0, 2], grad_diff_10x, "Abs Difference (10x)", cmap='gray', vmin=0, vmax=g_range)

# --- Row 2: Hessians ---
h_min = min(hess_analytic.min(), hess_fd.min())
h_max = max(hess_analytic.max(), hess_fd.max())
h_range = h_max - h_min

plot_smooth(axes[1, 0], hess_analytic, "Analytic Hessian", vmin=h_min, vmax=h_max)
plot_smooth(axes[1, 1], hess_fd, "Finite Diff Hessian", vmin=h_min, vmax=h_max)
# Plot 10x Error against the Signal Range
plot_smooth(axes[1, 2], hess_diff_10x, "Abs Difference (10x)", cmap='gray', vmin=0, vmax=h_range)

plt.tight_layout()
plt.savefig("verification_smooth_10x.png", dpi=300)
plt.show()