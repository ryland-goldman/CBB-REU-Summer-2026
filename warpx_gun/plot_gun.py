"""
Figures for the WarpX RZ CESR-gun simulation (gun_sim.py).

Reads the applied gun field (warpx_gun/gun_field/gun_E.h5) and the openPMD beam
output under warpx_gun/diags/, and writes four figures to warpx_gun/results/:

  1. gun_field.png     — on-axis Ez(z) and implied potential of the scaled
                         CESR_gun.gdf map: the accelerating field the beam sees.
  2. beam_rz.png       — r–z distribution of the beam at three snapshots
                         (launch, mid-gun, exit): transport through the gun.
  3. energy_gain.png   — mean/max kinetic energy of the beam vs. ⟨z⟩, climbing
                         toward the ~150 keV gun voltage.
  4. exit_phase_space.png — longitudinal (z–KE) and the final energy spectrum.

Run with:
    conda run -n CBB python warpx_gun/plot_gun.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import openpmd_api as io
from openpmd_viewer import OpenPMDTimeSeries

MC2 = 0.51099895e3           # electron rest energy [keV]
GUN_FIELD = "warpx_gun/gun_field/gun_E.h5"
RESULTS = "warpx_gun/results"
os.makedirs(RESULTS, exist_ok=True)


def gamma_of(ux, uy, uz):
    """Lorentz γ from openPMD normalized momenta (γβ)."""
    return np.sqrt(1.0 + ux**2 + uy**2 + uz**2)


# ── Applied gun field (on axis) ───────────────────────────────────────────────
s = io.Series(GUN_FIELD, io.Access.read_only)
E = s.iterations[0].meshes["E"]
ez_map = E["z"].load_chunk()
s.flush()
ez_map = ez_map[0]                                   # (nr, nz), mode 0
dz_map = E.grid_spacing[1]
nz_map = ez_map.shape[1]
z_map = np.arange(nz_map) * dz_map
ez_axis = ez_map[0]                                  # r = 0 row
# Implied on-axis potential V(z) = -∫ Ez dz (referenced to the exit).
V_axis = -np.cumsum(ez_axis[::-1]) [::-1] * dz_map

fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.3), constrained_layout=True)
a1.plot(z_map * 1e3, ez_axis / 1e6, color="C3")
a1.axhline(0, color="k", lw=0.6)
a1.set_xlabel("z  [mm]"); a1.set_ylabel(r"$E_z$ on axis  [MV/m]")
a1.set_title("Applied gun field (scaled CESR_gun.gdf)")
a2.plot(z_map * 1e3, V_axis / 1e3, color="C0")
a2.set_xlabel("z  [mm]"); a2.set_ylabel("implied potential  [kV]")
a2.set_title("On-axis potential (cathode → exit)")
fig.savefig(f"{RESULTS}/gun_field.png", dpi=140)
print(f"wrote {RESULTS}/gun_field.png")

# ── Beam time series ──────────────────────────────────────────────────────────
ts = OpenPMDTimeSeries("warpx_gun/diags/particles")
iters = ts.iterations

t_ns, zmean, ke_mean, ke_max, n_live = [], [], [], [], []
for i, it in enumerate(iters):
    z, ux, uy, uz = ts.get_particle(
        ["z", "ux", "uy", "uz"], species="electrons", iteration=it)
    n_live.append(len(z))
    t_ns.append(ts.t[i] * 1e9)
    if len(z) == 0:
        zmean.append(np.nan); ke_mean.append(np.nan); ke_max.append(np.nan)
        continue
    ke = (gamma_of(ux, uy, uz) - 1.0) * MC2
    zmean.append(z.mean() * 1e3)
    ke_mean.append(ke.mean()); ke_max.append(ke.max())

t_ns = np.array(t_ns); zmean = np.array(zmean)
ke_mean = np.array(ke_mean); ke_max = np.array(ke_max)
print(f"beam: {n_live[0]} launched, {n_live[-1]} at last dump; "
      f"peak ⟨KE⟩ {np.nanmax(ke_mean):.1f} keV, max KE {np.nanmax(ke_max):.1f} keV")

# ── Fig 2: r–z at launch / mid / exit ─────────────────────────────────────────
live = [it for it, n in zip(iters, n_live) if n > 0]
picks = [live[0], live[len(live)//2], live[-1]] if len(live) >= 3 else live
fig, axs = plt.subplots(1, len(picks), figsize=(4.2*len(picks), 4.0),
                        constrained_layout=True, squeeze=False)
for ax, it in zip(axs[0], picks):
    z, x, y = ts.get_particle(["z", "x", "y"], species="electrons", iteration=it)
    r = np.sqrt(x**2 + y**2)
    ax.hist2d(z*1e3, r*1e3, bins=[120, 60], cmap="viridis",
              norm=LogNorm(), cmin=1)
    ti = ts.t[list(iters).index(it)] * 1e9
    ax.set_title(f"t = {ti:.2f} ns  (N={len(z)})")
    ax.set_xlabel("z  [mm]"); ax.set_ylabel("r  [mm]")
fig.suptitle("Beam transport through the gun (r–z)", fontsize=12)
fig.savefig(f"{RESULTS}/beam_rz.png", dpi=140)
print(f"wrote {RESULTS}/beam_rz.png")

# ── Fig 3: energy gain vs ⟨z⟩ ─────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
ok = np.isfinite(zmean)
ax.plot(zmean[ok], ke_mean[ok], "o-", color="C2", ms=3, label="mean KE")
ax.plot(zmean[ok], ke_max[ok], "^--", color="C1", ms=3, label="max KE")
ax.axhline(150.0, color="k", ls=":", label="150 keV (gun voltage)")
ax.set_xlabel("mean beam position  ⟨z⟩  [mm]")
ax.set_ylabel("kinetic energy  [keV]")
ax.set_title("Beam energy gain along the gun")
ax.legend()
fig.savefig(f"{RESULTS}/energy_gain.png", dpi=140)
print(f"wrote {RESULTS}/energy_gain.png")

# ── Fig 4: exit longitudinal phase space + energy spectrum ────────────────────
it_exit = live[-1]
z, ux, uy, uz = ts.get_particle(
    ["z", "ux", "uy", "uz"], species="electrons", iteration=it_exit)
ke = (gamma_of(ux, uy, uz) - 1.0) * MC2
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.3), constrained_layout=True)
a1.scatter(z*1e3, ke, s=2, alpha=0.2, color="C0")
a1.set_xlabel("z  [mm]"); a1.set_ylabel("KE  [keV]")
a1.set_title(f"Longitudinal phase space  (t = {ts.t[list(iters).index(it_exit)]*1e9:.2f} ns)")
a2.hist(ke, bins=60, color="C3", alpha=0.85)
a2.axvline(ke.mean(), color="k", ls="--", label=f"⟨KE⟩ = {ke.mean():.1f} keV")
a2.set_xlabel("KE  [keV]"); a2.set_ylabel("count")
a2.set_title("Energy spectrum at last dump"); a2.legend()
fig.savefig(f"{RESULTS}/exit_phase_space.png", dpi=140)
print(f"wrote {RESULTS}/exit_phase_space.png")
