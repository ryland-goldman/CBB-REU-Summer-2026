"""One-off driver: re-run the injector + linac_sec1 at Tan condition (i) for Deliverable B.

Tan cond (i): Preb1 50 kV @ phase-null, Preb2 150 kV @ crest, default gun/Sec1.
kV->kW (injector_sim scale = sqrt(1e3*Q*P/(2*pi*f)), V_gap = scale*V1J_KEV; inverted):
    Preb1 50 kV  (Q=3000) -> PREB1_KW = 5.8298 kW
    Preb2 150 kV (Q=4300) -> PREB2_KW = 36.6054 kW
Phase-null = crest +/- 90 deg => PREB1_PHI_OFF=-90 (Preb1 not accelerating, ~0 net kick);
Preb2 @ crest => PREB2_PHI_OFF=0. Writes injector/diags/cond_i and linac_sec1/diags/cond_i.

Two-pass note (plan §1.4): at PREB1_PHI_OFF=-90 the analytic Preb-1 mean kick
    kick1 = -cos(base + phi_off)*scale1*V1J = -cos(pi - pi/2)*... = 0,
so the Preb-2 two-segment timing (v_after_preb1 = v_beam) is already self-consistent — no
beta walk to correct. Single pass suffices; the printed Preb-2 Δφ diagnostic is checked.

Run (CBB env, repo root):  python pipeline/_run_cond_i.py
"""

import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
os.chdir(_ROOT)

import injector
import linac_sec1

PREB1_KW_CONDI = 5.8298     # 50 kV @ Q=3000
PREB2_KW_CONDI = 36.6054    # 150 kV @ Q=4300


def main():
    # The dense cond_i bunch (150 kV crest Preb-2) diverged the long-thin-box MLMG solve
    # ~20 cm short of the 2.03 m handoff in the first attempt (aborted at ⟨z⟩≈1.83 m, no
    # summary). Per CLAUDE.md: relax the solve (more iters, looser tol) and give the transit
    # more headroom so a dump lands at the handoff before the bunch drains into the absorber.
    injector.config(
        PREB1_KW=PREB1_KW_CONDI, PREB1_PHI_OFF=-90.0,   # 50 kV @ phase-null
        PREB2_KW=PREB2_KW_CONDI, PREB2_PHI_OFF=0.0,     # 150 kV @ crest
        OUTDIR="injector/diags/cond_i",
        MAX_ITERS=1000,                                 # more MLMG headroom for the dense bunch
        REQUIRED_PRECISION=3e-4,                        # looser tol (long-thin box, dense self-field)
        CFL=0.6,                                        # finer dt → smaller per-step field change
    )
    injector.run(plots=False)

    linac_sec1.config(
        INJECTOR_DIAG="injector/diags/cond_i/particles",
        OUTDIR="linac_sec1/diags/cond_i",
    )
    linac_sec1.run(plots=False)


if __name__ == "__main__":
    main()
