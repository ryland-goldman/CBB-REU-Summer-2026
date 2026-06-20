"""
Figures for the WarpX RZ CESR-gun simulation: reads gun/gun_field/gun_E.h5,
gun/diags/fields/ (phi, rho) and gun/diags/particles/, and writes six PNGs to
gun/results/ (field, r-z transport, energy gain, exit phase space, envelope,
self-field).

See gun/README.md for physics, parameters, and gotchas.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import openpmd_api as io
from openpmd_viewer import OpenPMDTimeSeries

from pipeline.beam_metrics import screen_profile
from pipeline.constants import MC2_EV
from gun.build_gun_field import GUN_VOLTAGE   # single-source the display voltage

MC2 = MC2_EV / 1e3           # electron rest energy [keV]
GUN_FIELD = "gun/gun_field/gun_E.h5"
RESULTS = "gun/results"


def gamma_of(ux, uy, uz):
    """Lorentz γ from openPMD normalized momenta (u = γβ)."""
    return np.sqrt(1.0 + ux**2 + uy**2 + uz**2)


def main():
    os.makedirs(RESULTS, exist_ok=True)

    # ── Applied gun field (on axis) ───────────────────────────────────────────
    s = io.Series(GUN_FIELD, io.Access.read_only)
    E = s.iterations[0].meshes["E"]
    ez_map = E["z"].load_chunk()
    s.flush()
    ez_map = ez_map[0]                                   # (nr, nz), mode 0
    dz_map = E.grid_spacing[1]
    nz_map = ez_map.shape[1]
    z_map = np.arange(nz_map) * dz_map
    ez_axis = ez_map[0]                                  # r = 0 row
    # Implied on-axis potential, exit-referenced (V(exit)=0): V(z) = +∫_z^exit Ez dz'.
    V_axis = np.cumsum(ez_axis[::-1])[::-1] * dz_map

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

    # ── Beam time series ──────────────────────────────────────────────────────
    ts = OpenPMDTimeSeries("gun/diags/particles")
    iters = ts.iterations

    t_ns, zmean, ke_mean, ke_max, n_live = [], [], [], [], []
    # Pooled (all-dump) columns for the per-z profile, kept WITH the particle id so the profile
    # is reconstructed as fixed-z virtual SCREENS (id-track each particle's crossing of every
    # z-plane), not a z-histogram: a per-dump ⟨z⟩ aggregate is meaningless for the timed quasi-DC
    # stream (see README -> "Beam representation"). See screen_profile().
    id_pool, z_pool, ke_pool, x_pool, ux_pool, w_pool = [], [], [], [], [], []
    for i, it in enumerate(iters):
        idp, z, x, ux, uy, uz, w = ts.get_particle(
            ["id", "z", "x", "ux", "uy", "uz", "w"], species="electrons", iteration=it)
        n_live.append(len(z))
        t_ns.append(ts.t[i] * 1e9)
        if len(z) == 0:
            zmean.append(np.nan); ke_mean.append(np.nan); ke_max.append(np.nan)
            continue
        ke = (gamma_of(ux, uy, uz) - 1.0) * MC2
        zmean.append(z.mean() * 1e3)
        ke_mean.append(ke.mean()); ke_max.append(ke.max())
        id_pool.append(idp.astype(np.int64))
        z_pool.append(z); ke_pool.append(ke); x_pool.append(x); ux_pool.append(ux); w_pool.append(w)

    t_ns = np.array(t_ns); zmean = np.array(zmean)
    ke_mean = np.array(ke_mean); ke_max = np.array(ke_max)
    print(f"beam: {n_live[0]} launched, {n_live[-1]} at last dump; "
          f"peak ⟨KE⟩ {np.nanmax(ke_mean):.1f} keV, max KE {np.nanmax(ke_max):.1f} keV")

    # ── Per-z PROFILE via fixed-z virtual SCREENS (id-tracked plane crossings) ──
    zc = kez = kezmax = sxz = enz = np.array([])
    if z_pool:
        idP = np.concatenate(id_pool); zP = np.concatenate(z_pool)
        xP = np.concatenate(x_pool); uxP = np.concatenate(ux_pool)
        keP = np.concatenate(ke_pool); wP = np.concatenate(w_pool)
        screens, prof = screen_profile(
            idP, zP, wP, {"x": xP, "ux": uxP, "ke": keP}, emit_pairs=[("x", "ux")])
        zc = screens * 1e3                                  # [mm]
        sxz = prof["rms"]["x"] * 1e3                        # σ_x [mm]
        enz = prof["emit"][("x", "ux")] * 1e6              # εn,x [mm·mrad]
        kez = prof["mean"]["ke"]; kezmax = prof["max"]["ke"]

    # ── Fig 2: r–z at launch / mid / exit ─────────────────────────────────────
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

    # ── Fig 3: energy gain vs z (fixed-z virtual screens) ─────────────────────
    fig, ax = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
    okz = np.isfinite(kez)
    ax.plot(zc[okz], kez[okz], "o-", color="C2", ms=3, label="mean KE")
    ax.plot(zc[okz], kezmax[okz], "^--", color="C1", ms=3, label="max KE")
    ax.axhline(GUN_VOLTAGE / 1e3, color="k", ls=":", label=f"{GUN_VOLTAGE/1e3:.0f} keV (gun voltage)")
    ax.set_xlabel("beam position  z  [mm]")
    ax.set_ylabel("kinetic energy  [keV]")
    ax.set_title("Beam energy gain along the gun  (fixed-z virtual screens)")
    ax.legend()
    fig.savefig(f"{RESULTS}/energy_gain.png", dpi=140)
    print(f"wrote {RESULTS}/energy_gain.png")

    # ── Fig 4: exit longitudinal phase space + energy spectrum ────────────────
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

    # ══════════════════════════════════════════════════════════════════════════
    # Fig 5: beam_envelope.png — per-plane RMS size σ_x and emittance εn,x vs z
    # ──────────────────────────────────────────────────────────────────────────
    okz = np.isfinite(sxz) & np.isfinite(enz)
    fig, ax = plt.subplots(figsize=(7.6, 4.6), constrained_layout=True)
    l1, = ax.plot(zc[okz], sxz[okz], "o-", color="C0", ms=3,
                  label=r"RMS size  $\sigma_x=\sqrt{\langle x^2\rangle}$")
    ax.set_xlabel("beam position  z  [mm]")
    ax.set_ylabel(r"per-plane RMS size  $\sigma_x$  [mm]", color="C0")
    ax.tick_params(axis="y", labelcolor="C0")
    ax.set_title("Transverse envelope and emittance along the gun  (fixed-z virtual screens)")
    ax2 = ax.twinx()
    l2, = ax2.plot(zc[okz], enz[okz], "s--", color="C3", ms=3,
                   label=r"norm. emittance  $\varepsilon_{n,x}$")
    ax2.set_ylabel(r"$\varepsilon_{n,x}$  [mm·mrad]", color="C3")
    ax2.tick_params(axis="y", labelcolor="C3")
    ax.legend(handles=[l1, l2], loc="best")
    fig.savefig(f"{RESULTS}/beam_envelope.png", dpi=140)
    print(f"wrote {RESULTS}/beam_envelope.png")

    # ══════════════════════════════════════════════════════════════════════════
    # Fig 6: space_charge.png — the beam SELF-FIELD (dumped ρ and φ) at near launch
    # ══════════════════════════════════════════════════════════════════════════
    fs = io.Series("gun/diags/fields/openpmd_%06T.h5", io.Access.read_only)
    field_iters = [int(k) for k in fs.iterations]
    particle_iters = list(iters)
    zmean_by_it = {int(it): zm for it, zm in zip(particle_iters, zmean)}
    Z_TARGET = 0.4                                  # mm — near-cathode but off the wall
    cand = [it for it in field_iters
            if it in zmean_by_it and np.isfinite(zmean_by_it[it])
            and zmean_by_it[it] > 0.0]
    if not cand:
        print("skipping space_charge.png: no field snapshot with a positive-⟨z⟩ beam")
        return
    it_sc = min(cand, key=lambda it: abs(zmean_by_it[it] - Z_TARGET))
    zmean_sc = zmean_by_it[it_sc]

    itr = fs.iterations[it_sc]
    t_sc = itr.time * itr.time_unit_SI * 1e9            # snapshot time [ns]
    rho_m = itr.meshes["rho"]; phi_m = itr.meshes["phi"]
    rho = rho_m[io.Mesh_Record_Component.SCALAR].load_chunk()
    phi = phi_m[io.Mesh_Record_Component.SCALAR].load_chunk()
    fs.flush()
    rho = rho[0]; phi = phi[0]                           # (nz, nr) mode-0 half-plane
    dz_sc, dr_sc = rho_m.grid_spacing                    # axis order ['z','r']
    nz_sc, nr_sc = rho.shape
    extent = [0.0, nz_sc * dz_sc * 1e3, 0.0, nr_sc * dr_sc * 1e3]   # [z0,z1,r0,r1] mm
    rho_img = rho.T * 1e6                                 # C/m³ → µC/m³ for readability
    phi_img = phi.T                                       # V

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.6, 6.4), constrained_layout=True,
                                 sharex=True)
    rmax = np.nanmax(np.abs(rho_img)) or 1.0
    im1 = a1.imshow(rho_img, origin="lower", extent=extent, aspect="auto",
                    cmap="RdBu_r", vmin=-rmax, vmax=rmax)
    a1.set_ylabel("r  [mm]")
    a1.set_title(f"Beam self charge density  ρ(r, z)   "
                 f"(t = {t_sc:.3f} ns,  ⟨z⟩ = {zmean_sc:.3f} mm)")
    cb1 = fig.colorbar(im1, ax=a1); cb1.set_label(r"ρ  [µC/m$^3$]")
    im2 = a2.imshow(phi_img, origin="lower", extent=extent, aspect="auto",
                    cmap="viridis")
    a2.set_xlabel("z  [mm]"); a2.set_ylabel("r  [mm]")
    a2.set_title("Space-charge potential well  φ(r, z)   (beam self-field only)")
    cb2 = fig.colorbar(im2, ax=a2); cb2.set_label("φ  [V]")
    z_zoom = max(5.0, 6.0 * zmean_sc)                    # mm; show a few mm past ⟨z⟩
    a1.set_xlim(0.0, min(z_zoom, extent[1]))
    fig.suptitle("Beam self-field near launch (separate from the applied gun field)",
                 fontsize=12)
    fig.savefig(f"{RESULTS}/space_charge.png", dpi=140)
    print(f"wrote {RESULTS}/space_charge.png")


if __name__ == "__main__":
    main()
