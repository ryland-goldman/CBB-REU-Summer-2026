"""
SPEBS gun in IMPACT-T: bunched cathode emission across the 50 kV Pierce-gun gap,
driven from Python via lume-impact.

This is the IMPACT-T companion to ``warpx_cathode/spebs_gun.py``. WarpX solved the
SAME gun (V = 50 kV, gap d = 23.1 mm, R = 4 mm cathode) as a *steady-state*
space-charge-limited diode and validated the Child-Langmuir current limit. IMPACT-T
is a photoinjector / bunched-beam tracking code, so here we model the complementary
picture: a finite electron *bunch* is emitted from the cathode over a short pulse and
tracked across the gap with 3D space charge **and the cathode image charge**
(``Flagimg``). The deliverables are the bunch dynamics — energy gain, transverse
beam size, thermal-emittance growth — and how the cathode image charge changes them.

Field model: a uniform DC accelerating field E0 = V/d = 2.165 MV/m is built as an
on-axis ``FieldMesh`` (flat across the gap, smooth cosine ramps that sit *inside* the
cathode at z<0 and *past* the anode, so the tracked beam only ever sees the flat
region) and attached as a standard ``solrf`` (type 105) gun element with RF
frequency 0. The field map is written in IMPACT-T's discrete derivative format, which
requires the solrf file ID to be > 1000 (otherwise IMPACT-T misreads it as Fourier
coefficients and the field collapses to ~0).

Usage (run once per image-charge setting):
    python impactt_cathode/spebs_gun_impactt.py img    # Flagimg = 1  (image charge ON)
    python impactt_cathode/spebs_gun_impactt.py noimg  # Flagimg = 0  (image charge OFF)

Each run drops an ImpactT.in + rfdata1001 (the actual IMPACT-T input it executed) and
a stats .json into impactt_cathode/results/. Run plot_impactt.py afterwards to render
the comparison figures.
"""

import json
import os
import sys

import numpy as np
from impact import Impact
from pmd_beamphysics import FieldMesh

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

# ── CLI: image-charge tag ────────────────────────────────────────────────────────
tag = sys.argv[1] if len(sys.argv) > 1 else "img"
flag_img = 0 if tag == "noimg" else 1

# ── SPEBS gun parameters (match warpx_cathode/spebs_gun.py) ──────────────────────
V_anode = 50.0e3        # 50 kV bias
gap_d = 23.1e-3         # cathode -> anode gap [m]
R_cathode = 4.0e-3      # 8 mm-diameter cathode -> 4 mm radius
T_cathode = 1500.0      # thermionic cathode temperature [K]
E0 = V_anode / gap_d    # uniform DC field magnitude = 2.165 MV/m

# physical constants (eV-based, matching IMPACT-T conventions)
m_e_eV = 510998.95      # electron rest mass [eV]
kB = 1.380649e-23
c = 299792458.0
q_e = 1.602176634e-19

# thermal momentum spread gamma*beta = sqrt(kB T / m_e c^2)
v_th_over_c = np.sqrt(kB * T_cathode / (9.1093837015e-31 * c**2))

# bunch charge: SPEBS runs ~1 A continuously; we emit a 0.1 nC slice
# (Bcurr / Bfreq = 1 A / 1e10 Hz = 0.1 nC) over a 120 ps pulse.
Bcurr, Bfreq = 1.0, 1.0e10
Q_bunch = Bcurr / Bfreq

# Child-Langmuir limit, for reference / comparison with the WarpX run
J_CL = (4.0 / 9.0) * 8.8541878128e-12 * np.sqrt(2.0 * q_e / 9.1093837015e-31) \
    * V_anode**1.5 / gap_d**2

print(f"[{tag}] V={V_anode/1e3:.0f} kV  gap={gap_d*1e3:.1f} mm  R={R_cathode*1e3:.0f} mm  "
      f"E0={E0/1e6:.3f} MV/m  Flagimg={flag_img}")
print(f"[{tag}] bunch charge Q={Q_bunch*1e9:.3f} nC,  J_CL={J_CL:.3e} A/m^2,  "
      f"thermal gamma*beta={v_th_over_c:.2e}")

# ── Build the uniform DC gun field as an on-axis FieldMesh ───────────────────────
# Ez = -E0 (negative so the force -eE accelerates electrons toward +z). The field is
# flat across the gap; the ramp-up sits at z<0 (inside the cathode) so that at the
# emission plane z=0 the field is already flat (dEz/dz=0 -> no spurious radial kick).
#
# The flat region must extend well past the anode: the bunch is long (early-emitted
# electrons accelerate longer, so the head runs several mm ahead of the centroid). If
# the head reaches the field ramp-down before the run stops, the ramp's dEz/dz gives
# off-axis radial kicks that artificially inflate the normalized emittance. We keep
# the field uniform out to 40 mm (the head reaches ~28 mm when the centroid hits the
# 23.1 mm anode) and ramp down only beyond there, outside the tracked region.
z = np.linspace(-0.004, 0.044, 961)
Ez = np.full_like(z, -E0)
up = z < -0.001
Ez[up] = -E0 * 0.5 * (1 - np.cos(np.pi * (z[up] + 0.004) / 0.003))
dn = z > 0.040
Ez[dn] = -E0 * 0.5 * (1 + np.cos(np.pi * (z[dn] - 0.040) / 0.004))

fmesh = FieldMesh.from_onaxis(z=z, Ez=Ez, frequency=0)
# discrete-derivative solrf map; file_id>1000 is REQUIRED for IMPACT-T to read it.
pkg = fmesh.to_impact_solrf(style="derivatives", zedge=-0.004, scale=E0, file_id=1001)
gun = pkg["ele"]
gun["L"] = 0.048
gun["name"] = "dc_gun"

# ── Assemble the IMPACT-T run via lume-impact ────────────────────────────────────
# Start from the bundled tesla deck only for a valid header skeleton, then overwrite
# every field that matters for this problem.
BASE = os.path.join(os.path.dirname(__import__("impact").__file__),
                    "tests", "input", "tesla_9cell_cavity", "ImpactT.in")
I = Impact(BASE)
I.input["lattice"] = [gun, {"type": "stop", "s": gap_d, "name": "stop_1"}]
I.input["fieldmaps"] = {gun["filename"]: pkg["fmap"]}

h = I.header
h.update(dict(
    Npcol=1, Nprow=1,
    Dt=1.0e-12, Ntstep=50000, Nbunch=1,
    Np=10000, Flagmap=1, Flagerr=0, Flagdiag=1,
    Flagimg=flag_img, Zimage=0.005,
    Nx=32, Ny=32, Nz=128, Flagbc=1, Xrad=0.010, Yrad=0.010, Perdlen=0.060,
    Flagdist=2, Rstartflg=0, Nemission=100, Temission=120e-12,
    Bcurr=Bcurr, Bkenergy=1.0, Bmass=m_e_eV, Bcharge=-1.0, Bfreq=Bfreq, Tini=0.0,
))
# transverse cathode spot + thermal momentum; longitudinal pulse slice
h["sigx(m)"], h["sigy(m)"] = 2.0e-3, 2.0e-3
h["sigpx"], h["sigpy"] = v_th_over_c, v_th_over_c
h["sigz(m)"], h["sigpz"] = 1.0e-4, 1.0e-3
# IMPORTANT: zero the longitudinal centroid (the tesla base ships a 10 MeV pz here).
h["zmu1(m)"], h["zmu2"] = 0.0, 0.0

# ── Run ──────────────────────────────────────────────────────────────────────────
I.verbose = False
print(f"[{tag}] running IMPACT-T ...")
I.run()
info = I.output["run_info"]
print(f"[{tag}] done in {info['run_time']:.1f} s  (error={info['error']})")

# Persist the exact deck + field map that ran, for inspection.
deck_dir = os.path.join(RESULTS, f"deck_{tag}")
os.makedirs(deck_dir, exist_ok=True)
I.write_input("ImpactT.in", path=deck_dir)

# ── Save stats + final particles for plotting ────────────────────────────────────
s = I.output["stats"]
final = I.particles["final_particles"]
out = {
    "tag": tag, "flag_img": flag_img,
    "E0": E0, "V": V_anode, "gap": gap_d, "R": R_cathode, "Q": Q_bunch,
    "z": s["mean_z"].tolist(),
    "KE": s["mean_kinetic_energy"].tolist(),
    "sigma_x": s["sigma_x"].tolist(),
    "sigma_y": s["sigma_y"].tolist(),
    "norm_emit_x": s["norm_emit_x"].tolist(),
    "norm_emit_y": s["norm_emit_y"].tolist(),
    "n_particle": s["n_particle"].tolist(),
    "final_KE_eV": float(s["mean_kinetic_energy"][-1]),
    "final_z": float(s["mean_z"][-1]),
    "transmission": float(s["n_particle"][-1] / s["n_particle"][0]),
    # final phase space (subsample for compact storage)
    "fp_x": final["x"][::5].tolist(),
    "fp_px": final["px"][::5].tolist(),
    "fp_z": final["z"][::5].tolist(),
    "fp_pz": final["pz"][::5].tolist(),
    "final_norm_emit_x": float(final["norm_emit_x"]),
}
with open(os.path.join(RESULTS, f"stats_{tag}.json"), "w") as f:
    json.dump(out, f)

print(f"[{tag}] final KE={out['final_KE_eV']/1e3:.2f} keV  (anode V={V_anode/1e3:.0f} kV)  "
      f"transmission={out['transmission']*100:.0f}%  "
      f"norm_emit_x={out['final_norm_emit_x']*1e6:.3f} um")
print(f"[{tag}] wrote results/stats_{tag}.json")
