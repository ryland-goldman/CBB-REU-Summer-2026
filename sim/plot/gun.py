"""
Figures for the WarpX RZ CESR-gun simulation over logs/diags/gun/. Writes PNGs to
logs/plots/gun/.

Figures: phase_space_z_KE, energy_spectrum, transverse_r_pr, evolution_vs_z (mean KE / eps_n,x /
sigma_x along the gun via fixed-z virtual screens), plus the stage-specific rich figures (on-axis
applied field gun_field, r–z transport beam_rz, and the beam self-field space_charge). See
docs/gun.md for the physics each figure shows.

main() runs ONLY the plotting; sim/gun.py runs the simulation.
"""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import openpmd_api as io
from openpmd_viewer import OpenPMDTimeSeries

from sim.helpers.tools import prepare_env
from sim.plot import common as px

CONFIG = "config/gun.yaml"
RESULTS = "logs/plots/gun"
GUN_FIELD = "fieldmaps/h5/gun_E.h5"
PARTICLES = "logs/diags/gun/particles"
FIELDS = "logs/diags/gun/fields"


def _last_populated(diag, species="electrons"):
    ts = OpenPMDTimeSeries(diag)
    for it in reversed([int(i) for i in ts.iterations]):
        try:
            x, = ts.get_particle(["x"], species=species, iteration=it)
        except Exception:
            continue
        if len(x):
            return it
    return int(ts.iterations[-1])


def _save(fig, name):
    fig.savefig(f"{RESULTS}/{name}.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {RESULTS}/{name}.png")


def applied_field_figure():
    """On-axis applied E_z and the potential it implies (cathode → exit)."""
    s = io.Series(GUN_FIELD, io.Access.read_only)
    mesh = s.iterations[0].meshes["E"]
    ez = mesh["z"].load_chunk()
    s.flush()
    ez_axis = ez[0][0]                                  # mode 0, r = 0 row
    dz = mesh.grid_spacing[1]
    z_mm = np.arange(ez_axis.size) * dz * 1e3
    # Exit-referenced on-axis potential V(z) = ∫_z^exit E_z dz' (V(exit) = 0).
    v_kv = np.cumsum(ez_axis[::-1])[::-1] * dz / 1e3

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.3), constrained_layout=True)
    a1.plot(z_mm, ez_axis / 1e6, color="C3")
    a1.axhline(0, color="k", lw=0.6)
    a1.set_xlabel("z  [mm]"); a1.set_ylabel(r"$E_z$ on axis  [MV/m]")
    a1.set_title("Applied gun field (scaled CESR_gun.gdf)")
    a2.plot(z_mm, v_kv, color="C0")
    a2.set_xlabel("z  [mm]"); a2.set_ylabel("implied potential  [kV]")
    a2.set_title("On-axis potential (cathode → exit)")
    _save(fig, "gun_field")


def transport_figure(ts, populated_iters):
    """Beam shape in r–z at launch / mid-gun / exit (dumps past the one-particle seed)."""
    picks = ([populated_iters[0], populated_iters[len(populated_iters) // 2], populated_iters[-1]]
             if len(populated_iters) >= 3 else populated_iters)
    fig, axs = plt.subplots(1, len(picks), figsize=(4.2 * len(picks), 4.0),
                            constrained_layout=True, squeeze=False)
    times = list(ts.iterations)
    for ax, it in zip(axs[0], picks):
        z, x, y = ts.get_particle(["z", "x", "y"], species="electrons", iteration=it)
        r = np.hypot(x, y)
        ax.hist2d(z * 1e3, r * 1e3, bins=[120, 60], cmap="viridis", norm=LogNorm(), cmin=1)
        ax.set_title(f"t = {ts.t[times.index(it)] * 1e9:.2f} ns  (N={len(z)})")
        ax.set_xlabel("z  [mm]"); ax.set_ylabel("r  [mm]")
    fig.suptitle("Beam transport through the gun (r–z)", fontsize=12)
    _save(fig, "beam_rz")


def self_field_figure(ts, zmean_by_it):
    """Beam self-field ρ and φ (the dumped self-consistent fields) near launch."""
    fs = io.Series(f"{FIELDS}/openpmd_%06T.h5", io.Access.read_only)
    z_target = 0.4e-3                                    # near-cathode but off the wall [m]
    cand = [it for it in (int(k) for k in fs.iterations)
            if zmean_by_it.get(it, np.nan) > 0.0]
    if not cand:
        print("skipping space_charge.png: no field snapshot with a positive-⟨z⟩ beam")
        return
    it = min(cand, key=lambda it: abs(zmean_by_it[it] - z_target))
    zmean = zmean_by_it[it]

    itr = fs.iterations[it]
    t_ns = itr.time * itr.time_unit_SI * 1e9
    rho_m, phi_m = itr.meshes["rho"], itr.meshes["phi"]
    rho = rho_m[io.Mesh_Record_Component.SCALAR].load_chunk()
    phi = phi_m[io.Mesh_Record_Component.SCALAR].load_chunk()
    fs.flush()
    rho, phi = rho[0].T, phi[0].T                       # (nz, nr) → (nr, nz) for imshow
    dz, dr = rho_m.grid_spacing                         # axis order [z, r]
    extent = [0.0, rho.shape[1] * dz * 1e3, 0.0, rho.shape[0] * dr * 1e3]

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.6, 6.4), constrained_layout=True, sharex=True)
    rho_uc = rho * 1e6                                  # C/m³ → µC/m³
    rmax = np.nanmax(np.abs(rho_uc)) or 1.0
    im1 = a1.imshow(rho_uc, origin="lower", extent=extent, aspect="auto",
                    cmap="RdBu_r", vmin=-rmax, vmax=rmax)
    a1.set_ylabel("r  [mm]")
    a1.set_title(f"Beam self charge density  ρ(r, z)   (t = {t_ns:.3f} ns,  ⟨z⟩ = {zmean*1e3:.3f} mm)")
    fig.colorbar(im1, ax=a1, label=r"ρ  [µC/m$^3$]")
    im2 = a2.imshow(phi, origin="lower", extent=extent, aspect="auto", cmap="viridis")
    a2.set_xlabel("z  [mm]"); a2.set_ylabel("r  [mm]")
    a2.set_title("Space-charge potential well  φ(r, z)   (beam self-field only)")
    fig.colorbar(im2, ax=a2, label="φ  [V]")
    a1.set_xlim(0.0, min(max(5.0, 6.0 * zmean * 1e3), extent[1]))
    fig.suptitle("Beam self-field near launch (separate from the applied gun field)", fontsize=12)
    _save(fig, "space_charge")


def main():
    prepare_env()
    from warpx import WarpX
    os.makedirs(RESULTS, exist_ok=True)

    w = WarpX(input_file=CONFIG, path="logs/diags/gun")
    it = _last_populated(PARTICLES)

    # Generic phase-space / spectrum figures (lume-warpx helpers + shared plot common).
    w.load_output(diag_dir=PARTICLES)
    pg = w._particle_group(iteration=it)
    for name, fig in [
        ("phase_space_z_KE", w.plot2D("z", "kinetic_energy", iteration=it)),
        ("energy_spectrum",  px.energy_spectrum(pg)),
    ]:
        _save(fig, name)

    # Stage-specific rich figures (raw openPMD: applied field, transport, evolution, self-field).
    ts = OpenPMDTimeSeries(PARTICLES)
    live_iters, populated_iters, zmean_by_it = [], [], {}
    for itr in ts.iterations:
        z, = ts.get_particle(["z"], species="electrons", iteration=itr)
        if len(z):
            live_iters.append(itr)
            zmean_by_it[int(itr)] = z.mean()
        if len(z) > 50:                                  # past the one-particle release seed
            populated_iters.append(itr)

    x, y, ux, uy, wgt = ts.get_particle(["x", "y", "ux", "uy", "w"],
                                        species="electrons", iteration=it)
    _save(px.transverse_rpr(x, y, ux, uy, wgt,
                            title="Gun exit transverse phase space  (r, p_r)"),
          "transverse_r_pr")

    applied_field_figure()
    if live_iters:
        transport_figure(ts, populated_iters or live_iters)
        z_m, ke, emit, sigma, q_pc = px.evolution_screens(px.pool_trajectories(ts, live_iters))
        _save(px.evolution_vs_z(z_m, ke, emit, sigma, charge_pc=q_pc,
                                title="Beam evolution along the gun  (fixed-z virtual screens)"),
              "evolution_vs_z")
        self_field_figure(ts, zmean_by_it)


if __name__ == "__main__":
    main()
