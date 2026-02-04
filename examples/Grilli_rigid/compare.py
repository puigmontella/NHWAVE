import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
path='results/'


shift=0.3
data_WG1 =pd.read_csv(path+'probe_0001', delim_whitespace=True, header=None)
data_WG2 =pd.read_csv(path+'probe_0002', delim_whitespace=True, header=None)
data_WG3 =pd.read_csv(path+'probe_0003', delim_whitespace=True, header=None)
data_WG4 =pd.read_csv(path+'probe_0003', delim_whitespace=True, header=None)

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





# ##################### profile ########


# from pathlib import Path

# def load_Ha_any(path, N=None, M=None):
#     """
#     Returns either:
#       - ("matrix", Ha2D) if file is N x M (or can be reshaped to N x M)
#       - ("xyz", (x,y,Ha)) if file has 3 columns
#     """
#     a = np.loadtxt(path)

#     # xyz format: 3 columns
#     if a.ndim == 2 and a.shape[1] == 3:
#         x, y, Ha = a[:, 0], a[:, 1], a[:, 2]
#         return "xyz", (x, y, Ha)

#     # matrix format
#     if a.ndim == 2:
#         return "matrix", a

#     # 1D: try reshape if N,M provided
#     if a.ndim == 1 and (N is not None) and (M is not None) and a.size == N*M:
#         return "matrix", a.reshape((N, M))

#     raise ValueError(f"Unrecognized Ha format in {path}")

# def profile_from_matrix(Ha2D, dx, y_index=None, y_mean=False):
#     N, M = Ha2D.shape
#     x = (np.arange(M) + 0.5) * dx  # cell-center x

#     if y_mean:
#         h = Ha2D.mean(axis=0)       # average over y
#     else:
#         if y_index is None:
#             y_index = N // 2        # midline by default
#         h = Ha2D[y_index, :]

#     return x, h

# def profile_from_xyz(x, y, Ha, y_target=None):
#     """
#     Extract a cross-section at y_target (nearest y-level).
#     If y_target is None, uses the median y level.
#     """
#     # pick a y_target if not given
#     if y_target is None:
#         y_target = np.median(y)

#     # infer a reasonable tolerance from unique y spacing
#     yu = np.unique(np.round(y, 12))
#     if len(yu) > 1:
#         dy = np.median(np.diff(np.sort(yu)))
#         tol = 0.25 * dy
#     else:
#         tol = 1e-12

#     mask = np.abs(y - y_target) <= tol
#     xs = x[mask]
#     hs = Ha[mask]

#     # sort by x
#     idx = np.argsort(xs)
#     return xs[idx], hs[idx], y_target

# # -------------------
# # user settings
# # -------------------
# fdir = Path("results")
# files = ["Ha_00001", "Ha_00031", "Ha_00137"]   # pick your outputs
# dx = 0.0125                                     # match DX in input.txt
# N = 10                                         # match Nglob (only needed for reshape check)
# M = 500                                        # match Mglob

# plot_intv = 0.02                               # PLOT_INTV (optional, for labels)
# use_y_mean = False                             # True = average over y instead of one cross-section
# y_index = None                                 # e.g. 45 for N=90, or None for midline
# y_target = None                                # for xyz format (set a value if you want a specific y)

# plt.figure()

# for fname in files:
#     path = fdir / fname
#     kind, data = load_Ha_any(path, N=N, M=M)

#     if kind == "matrix":
#         Ha2D = data
#         xprof, hprof = profile_from_matrix(Ha2D, dx, y_index=y_index, y_mean=use_y_mean)

#     else:  # xyz
#         x, y, Ha = data
#         xprof, hprof, y_used = profile_from_xyz(x, y, Ha, y_target=y_target)

#     # label with time inferred from file number (works if Ha_00001 corresponds to 1*PLOT_INTV)
#     try:
#         k = int(fname.split("_")[1])
#         label = f"t={k*plot_intv:.3f}s ({fname})"
#     except Exception:
#         label = fname

#     plt.plot(xprof, hprof, label=label)

# plt.xlabel("x [m]")
# plt.ylabel("Ha [m]")
# plt.xlim([0, 1])
# plt.grid(True)
# plt.legend()
# plt.tight_layout()
plt.show()






















plt.show()




