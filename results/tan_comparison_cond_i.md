# Chapter 10 comparison — cond_i

**Operating point:** Tan condition (i) — Preb1 50 kV @ phase-null (PREB1_PHI_OFF=−90, 5.83 kW), Preb2 150 kV @ crest (PREB2_PHI_OFF=0, 36.61 kW); achieved V_gap 50.0 / 150.0 kV. **This IS Tan's operating point — the repo numbers below may be compared NUMERICALLY to Table 10.2** (within the plan §6 tolerances; σ_E and capture diverge in a documented direction).


### A at Gun
| location | zbar_m | Ebar_MeV | sigE_MeV | sigz_deg@214 | sigz_deg@2856 | capture_in_bucket_pct | capture_all_buckets_pct |
| --- | --- | --- | --- | --- | --- | --- | --- |
| A at Gun | 0.0467 | 0.146 | 0.005 | 1.07 | 14.25 |  |  |
| Tan (published) |  | 0.150 | 0.000 | 31.70 | 428.80 | 100.0 | 100.0 |
| Δ (repo−Tan) |  | -0.004 | 0.005 | -30.631 | -414.545 |  |  |

### B before Preb2
| location | zbar_m | Ebar_MeV | sigE_MeV | sigz_deg@214 | sigz_deg@2856 | capture_in_bucket_pct | capture_all_buckets_pct |
| --- | --- | --- | --- | --- | --- | --- | --- |
| B before Preb2 | 1.1634 | 0.141 | 0.015 | 16.20 | 216.00 |  |  |
| Tan (published) |  | 0.139 | 0.023 | 14.40 | 194.80 | 100.0 | 100.0 |
| Δ (repo−Tan) |  | 0.002 | -0.008 | 1.799 | 21.200 |  |  |

### C before Sec1
| location | zbar_m | Ebar_MeV | sigE_MeV | sigz_deg@214 | sigz_deg@2856 | capture_in_bucket_pct | capture_all_buckets_pct |
| --- | --- | --- | --- | --- | --- | --- | --- |
| C before Sec1 | 2.0132 | 0.281 | 0.024 | 20.84 | 277.87 | n/a (pre-iris) | n/a (pre-iris) |
| Tan (published) |  | 0.253 | 0.043 | 5.57 | 75.30 | 100.0 | 100.0 |
| Δ (repo−Tan) |  | 0.028 | -0.019 | 15.268 | 202.565 |  |  |

### D after Sec1
| location | zbar_m | Ebar_MeV | sigE_MeV | sigz_deg@214 | sigz_deg@2856 | capture_in_bucket_pct | capture_all_buckets_pct |
| --- | --- | --- | --- | --- | --- | --- | --- |
| D after Sec1 | 3.3993 | 26.274 | 5.343 |  | 163.87 | 2.34 | 2.90 |
| Tan (published) |  | 27.200 | 3.500 |  | 11.10 | 89.4 | 96.8 |
| Δ (repo−Tan) |  | -0.926 | 1.843 |  | 152.766 |  |  |

## Capture bookkeeping (denominators)

- `q_injected_C` (honest, pre-iris) = 6.7576e-10 C
- `q_in_bore_C` (post-iris)         = 9.6200e-11 C
- `q_in_domain_C`                   = 1.5282e-10 C
- after-Sec1 in-bucket / q_in_bore  = 16.47% (separates iris loss from capture loss)

## Footnotes

- **σ_z internal ratio:** the repo @2856/@214 columns derive from ONE σ_t, so their ratio is BY CONSTRUCTION F_2856/F_214 = 13.335 — an internal-consistency check only. Tan's own ratio is 13.52 (a different spatial-σ_z definition with a per-location β); the 13.335 check is NOT agreement with Tan.
- **Denominator mismatch:** Tan's 100% (upstream) and 89.4/96.8 (postSec1) reference the gun-emitted bunch with no loss yet. The repo's `q_injected_C` already sits past the converging halo and includes the real 9.547 mm iris scrape. Repo capture = 2.3% in-bucket / 2.9% all-buckets (both vs `q_injected_C`); iris transmission = 14.2% (`q_in_bore`/`q_injected`) is a SEPARATE quantity (do NOT conflate it with capture). The repo numbers are far below Tan's BY CONSTRUCTION (γ² self-field + real iris + operating point). Annotate, do NOT 'fix'.
- **Location C capture = n/a (pre-iris):** C is the pre-scrape injector population; its charge differs from D's iris survivors, so the C→D ratio is not a clean repo capture fraction.
- **after-Sec1 σ_z@214 blank:** Tan leaves it blank; the captured-core σ_z is reported @2856 only.
