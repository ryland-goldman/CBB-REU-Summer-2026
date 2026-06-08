"""
Figures for the WarpX RZ CESR-gun simulation (gun_sim.py).

Reads the applied gun field (gun/gun_field/gun_E.h5), the dumped beam
self-field (gun/diags/fields/: phi, rho), and the openPMD beam output
(gun/diags/particles/), and writes six figures to gun/results/:

  1. gun_field.png     — on-axis Ez(z) and implied potential of the scaled
                         CESR_gun.gdf map: the accelerating field the beam sees.
  2. beam_rz.png       — r–z distribution of the beam at three snapshots
                         (launch, mid-gun, exit): transport through the gun.
  3. energy_gain.png   — mean/max kinetic energy of the beam vs. ⟨z⟩, climbing
                         toward the ~150 keV gun voltage.
  4. exit_phase_space.png — longitudinal (z–KE) and the final energy spectrum.
  5. beam_envelope.png — per-plane RMS size σ_x and normalized transverse emittance
                         εn,x vs. ⟨z⟩: the near-cathode focusing, quantified.
  6. space_charge.png  — r–z maps of the beam self charge density ρ and its
                         space-charge potential well φ at a near-launch snapshot.

Run with:
    conda run -n CBB python gun/plot_gun.py
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
GUN_FIELD = "gun/gun_field/gun_E.h5"
GUN_VOLTAGE = 150.0e3        # [V]; module-level so gun.config(GUN_VOLTAGE=...) is reflected
                             # in the energy-gain reference line (mirrors gun_sim.py)
RESULTS = "gun/results"


def gamma_of(ux, uy, uz):
    """Lorentz γ from openPMD normalized momenta (γβ)."""
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
    # Implied on-axis potential, exit-referenced (V(exit) = 0):
    #   V(z) = -∫_exit^z Ez dz' = +∫_z^exit Ez dz'.
    # With the accelerating field (Ez < 0 on axis) this gives V ≈ -GUN_VOLTAGE at the
    # cathode rising to 0 at the exit — the physical negative-cathode potential.
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
    sig_r, emit_nx = [], []           # for beam_envelope.png (filled in the same loop)
    for i, it in enumerate(iters):
        z, x, ux, uy, uz = ts.get_particle(
            ["z", "x", "ux", "uy", "uz"], species="electrons", iteration=it)
        n_live.append(len(z))
        t_ns.append(ts.t[i] * 1e9)
        if len(z) == 0:
            zmean.append(np.nan); ke_mean.append(np.nan); ke_max.append(np.nan)
            sig_r.append(np.nan); emit_nx.append(np.nan)
            continue
        ke = (gamma_of(ux, uy, uz) - 1.0) * MC2
        zmean.append(z.mean() * 1e3)
        ke_mean.append(ke.mean()); ke_max.append(ke.max())
        # ── Envelope diagnostics (one transverse plane; RZ ⇒ x and y are equivalent) ──
        # Per-plane RMS size σ_x = sqrt(⟨x²⟩), reported in mm. This is the SINGLE-PLANE RMS
        # (the convention paired with εn,x below and used by plot_chain / linac), NOT the
        # radial RMS: for an axisymmetric beam √⟨r²⟩ = √2·σ_x, so do not label this σ_r.
        sig_r.append(np.sqrt(np.mean(x**2)) * 1e3)
        # Normalized transverse emittance εn,x = sqrt(⟨x²⟩⟨ux²⟩ − ⟨x·ux⟩²): a phase-space
        # area that is invariant under linear (and acceleration) forces; growth ⇒ nonlinear
        # space-charge / field aberrations. ux is the normalized momentum γβx, so x[m]·ux
        # gives [m·rad]; ×1e6 → mm·mrad.
        e2 = np.mean(x**2) * np.mean(ux**2) - np.mean(x * ux)**2
        emit_nx.append(np.sqrt(max(e2, 0.0)) * 1e6)

    t_ns = np.array(t_ns); zmean = np.array(zmean)
    ke_mean = np.array(ke_mean); ke_max = np.array(ke_max)
    sig_r = np.array(sig_r); emit_nx = np.array(emit_nx)
    print(f"beam: {n_live[0]} launched, {n_live[-1]} at last dump; "
          f"peak ⟨KE⟩ {np.nanmax(ke_mean):.1f} keV, max KE {np.nanmax(ke_max):.1f} keV")

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

    # ── Fig 3: energy gain vs ⟨z⟩ ─────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
    ok = np.isfinite(zmean)
    ax.plot(zmean[ok], ke_mean[ok], "o-", color="C2", ms=3, label="mean KE")
    ax.plot(zmean[ok], ke_max[ok], "^--", color="C1", ms=3, label="max KE")
    ax.axhline(GUN_VOLTAGE / 1e3, color="k", ls=":", label=f"{GUN_VOLTAGE/1e3:.0f} keV (gun voltage)")
    ax.set_xlabel("mean beam position  ⟨z⟩  [mm]")
    ax.set_ylabel("kinetic energy  [keV]")
    ax.set_title("Beam energy gain along the gun")
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
    # Fig 5: beam_envelope.png — per-plane RMS size σ_x and emittance εn,x vs ⟨z⟩
    # ──────────────────────────────────────────────────────────────────────────
    ok = np.isfinite(zmean) & np.isfinite(sig_r)
    fig, ax = plt.subplots(figsize=(7.6, 4.6), constrained_layout=True)
    l1, = ax.plot(zmean[ok], sig_r[ok], "o-", color="C0", ms=3,
                  label=r"RMS size  $\sigma_x=\sqrt{\langle x^2\rangle}$")
    ax.set_xlabel("mean beam position  ⟨z⟩  [mm]")
    ax.set_ylabel(r"per-plane RMS size  $\sigma_x$  [mm]", color="C0")
    ax.tick_params(axis="y", labelcolor="C0")
    ax.set_title("Transverse envelope and emittance along the gun")
    ax2 = ax.twinx()
    l2, = ax2.plot(zmean[ok], emit_nx[ok], "s--", color="C3", ms=3,
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
