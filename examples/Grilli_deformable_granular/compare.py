import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
path='results/'


shift=0.3
data_WG1 =pd.read_csv(path+'probe_0001', delim_whitespace=True, header=None)
data_WG2 =pd.read_csv(path+'probe_0002', delim_whitespace=True, header=None)
data_WG3 =pd.read_csv(path+'probe_0003', delim_whitespace=True, header=None)
data_WG4 =pd.read_csv(path+'probe_0004', delim_whitespace=True, header=None)

t_wg1, h_wg1 = data_WG1[0]+shift, data_WG1[1]
t_wg2, h_wg2 = data_WG2[0]+shift, data_WG2[1]
t_wg3, h_wg3 = data_WG3[0]+shift, data_WG3[1]
t_wg4, h_wg4 = data_WG4[0]+shift, data_WG4[1]


##################################################
# Load the CSV file into a pandas DataFrame
data = pd.read_csv('Experiments/Grilli2017_WG1.csv', header=0)
data.columns = data.columns.str.strip()
time_WG1 = data['x']
h_WG1 = data['y']/1000

data = pd.read_csv('Experiments/Grilli2017_WG2.csv', header=0)
data.columns = data.columns.str.strip()
time_WG2 = data['x']
h_WG2 = data['y']/1000

data = pd.read_csv('Experiments/Grilli2017_WG3.csv', header=0)
data.columns = data.columns.str.strip()
time_WG3 = data['x']
h_WG3 = data['y']/1000

data = pd.read_csv('Experiments/Grilli2017_WG4.csv', header=0)
data.columns = data.columns.str.strip()
time_WG4 = data['x']
h_WG4 = data['y']/1000




######################### Figure #########################
fig, (ax0, ax1, ax2, ax3) = plt.subplots(nrows=4, ncols=1,
                                figsize=(12, 6))
ax0.plot(time_WG1, h_WG1, label='Grilli et al. (2017)', color='k',linestyle="None", marker="o",fillstyle="none")
ax0.plot(t_wg1, h_wg1, label='NHWAVE', color='royalblue')
ax0.grid(True)

ax1.plot(time_WG2, h_WG2, label='Grilli et al. (2017)', color='k',linestyle="None", marker="o",fillstyle="none")
ax1.plot(t_wg2, h_wg2, label='NHWAVE', color='royalblue')
ax1.grid(True)

ax2.plot(time_WG3, h_WG3, label='Grilli et al. (2017)', color='k',linestyle="None", marker="o",fillstyle="none")
ax2.plot(t_wg3, h_wg3, label='NHWAVE', color='royalblue')
ax2.grid(True)

ax3.plot(time_WG4, h_WG4, label='Grilli et al. (2017)', color='k',linestyle="None", marker="o",fillstyle="none")
ax3.plot(t_wg4, h_wg4, label='NHWAVE', color='royalblue')
ax3.set_xlabel('Time [s]')
ax3.grid(True)


for ax in [ax0, ax1, ax2, ax3]:
    ax.set_ylabel('$\\eta$ [m]')
    # Set limits
    ax.set_xlim([-0.2, 10.2])
    ax.set_ylim([-0.011, 0.011])
    if ax!=ax3:
        ax.set_xticklabels([])

ax3.legend(loc="lower left",framealpha=1 )

plt.tight_layout()
plt.savefig("Grilli_comparison_WaveElevation.jpg",dpi=300)





##################### profile ########


from pathlib import Path

def load_Ha_any(path, N=None, M=None):
    """
    Returns either:
      - ("matrix", Ha2D) if file is N x M (or can be reshaped to N x M)
      - ("xyz", (x,y,Ha)) if file has 3 columns
    """
    a = np.loadtxt(path)

    # xyz format: 3 columns
    if a.ndim == 2 and a.shape[1] == 3:
        x, y, Ha = a[:, 0], a[:, 1], a[:, 2]
        return "xyz", (x, y, Ha)

    # matrix format
    if a.ndim == 2:
        return "matrix", a

    # 1D: try reshape if N,M provided
    if a.ndim == 1 and (N is not None) and (M is not None) and a.size == N*M:
        return "matrix", a.reshape((N, M))

    raise ValueError(f"Unrecognized Ha format in {path}")

def profile_from_matrix(Ha2D, dx, y_index=None, y_mean=False):
    N, M = Ha2D.shape
    x = (np.arange(M) + 0.5) * dx  # cell-center x

    if y_mean:
        h = Ha2D.mean(axis=0)       # average over y
    else:
        if y_index is None:
            y_index = N // 2        # midline by default
        h = Ha2D[y_index, :]

    return x, h

def profile_from_xyz(x, y, Ha, y_target=None):
    """
    Extract a cross-section at y_target (nearest y-level).
    If y_target is None, uses the median y level.
    """
    # pick a y_target if not given
    if y_target is None:
        y_target = np.median(y)

    # infer a reasonable tolerance from unique y spacing
    yu = np.unique(np.round(y, 12))
    if len(yu) > 1:
        dy = np.median(np.diff(np.sort(yu)))
        tol = 0.25 * dy
    else:
        tol = 1e-12

    mask = np.abs(y - y_target) <= tol
    xs = x[mask]
    hs = Ha[mask]

    # sort by x
    idx = np.argsort(xs)
    return xs[idx], hs[idx], y_target





# -------------------
# user settings
# -------------------


fdir = Path("results")
files    = ["depth_00001", "depth_00044", "depth_00064"]
filesEta = ["eta_00001",   "eta_00044",   "eta_00064"]

dx = 0.01
N  = 10
M  = 600

plot_intv  = 0.01
use_y_mean = False
y_index    = None
y_target   = None

def get_profile(kind, data, dx, y_index=None, y_mean=False, y_target=None):
    if kind == "matrix":
        Ha2D = data
        return profile_from_matrix(Ha2D, dx, y_index=y_index, y_mean=y_mean)  # (xprof, hprof)
    else:  # xyz
        x, y, Ha = data
        xprof, hprof, _ = profile_from_xyz(x, y, Ha, y_target=y_target)
        return xprof, hprof

def time_label(fname: str, plot_intv: float) -> str:
    try:
        k = int(fname.split("_")[1])
        return f"t={k * plot_intv:.3f}s"
    except Exception:
        return fname

# -------------------- CHOOSE COLORS HERE --------------------
#   same color per time (eta solid, depth dashed)
colors_by_time = ["tab:blue", "tab:orange", "tab:green"]  # length can be >= len(files)


# ------------------------------------------------------------

plt.figure()

for i, (f_depth, f_eta) in enumerate(zip(files, filesEta)):
    # --- load ---
    kind_d, data_d = load_Ha_any(fdir / f_depth, N=N, M=M)
    kind_e, data_e = load_Ha_any(fdir / f_eta,   N=N, M=M)

    x_d, depth = get_profile(kind_d, data_d, dx, y_index=y_index, y_mean=use_y_mean, y_target=y_target)
    x_e, eta   = get_profile(kind_e, data_e, dx, y_index=y_index, y_mean=use_y_mean, y_target=y_target)

    if (len(x_d) != len(x_e)) or (not np.allclose(x_d, x_e)):
        x = x_d
        eta = np.interp(x_d, x_e, eta)
    else:
        x = x_d

    label = time_label(f_depth, plot_intv)

    c_eta = c_depth = colors_by_time[i % len(colors_by_time)]


    # --- plot ---
    plt.plot(x, eta,    color=c_eta,   lw=2, label=f"{label}  η")
    plt.plot(x, -depth, color=c_depth, lw=2, ls="--", label=f"{label}  -h")

plt.xlabel("x [m]")
plt.ylabel("[m]")
plt.xlim([0, 1])
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()






