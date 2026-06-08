# Chapter 10 comparison — repo_default

**Operating point:** repo default 8 kW / 10 kW two-cavity (Preb1 −70°, Preb2 −45°). **QUALITATIVE comparison only — NOT Tan cond (i)** (50 kV @ phase-null / 150 kV @ crest). Tan's published numbers are a reference column, NOT a target.


### A at Gun
| location | zbar_m | Ebar_MeV | sigE_MeV | sigz_deg@214 | sigz_deg@2856 | capture_in_bucket_pct | capture_all_buckets_pct |
| --- | --- | --- | --- | --- | --- | --- | --- |
| A at Gun | 0.0467 | 0.146 | 0.005 | 1.07 | 14.25 |  |  |
| Tan (published) |  | 0.150 | 0.000 | 31.70 | 428.80 | 100.0 | 100.0 |
| Δ (repo−Tan) |  | -0.004 | 0.005 | -30.631 | -414.545 |  |  |

### B before Preb2
| location | zbar_m | Ebar_MeV | sigE_MeV | sigz_deg@214 | sigz_deg@2856 | capture_in_bucket_pct | capture_all_buckets_pct |
| --- | --- | --- | --- | --- | --- | --- | --- |
| B before Preb2 | 1.1682 | 0.160 | 0.014 | 14.32 | 191.00 |  |  |
| Tan (published) |  | 0.139 | 0.023 | 14.40 | 194.80 | 100.0 | 100.0 |
| Δ (repo−Tan) |  | 0.021 | -0.009 | -0.076 | -3.800 |  |  |

### C before Sec1
| location | zbar_m | Ebar_MeV | sigE_MeV | sigz_deg@214 | sigz_deg@2856 | capture_in_bucket_pct | capture_all_buckets_pct |
| --- | --- | --- | --- | --- | --- | --- | --- |
| C before Sec1 | 2.0308 | 0.220 | 0.006 | 16.07 | 214.30 | n/a (pre-iris) | n/a (pre-iris) |
| Tan (published) |  | 0.253 | 0.043 | 5.57 | 75.30 | 100.0 | 100.0 |
| Δ (repo−Tan) |  | -0.033 | -0.037 | 10.501 | 138.999 |  |  |

### D after Sec1
| location | zbar_m | Ebar_MeV | sigE_MeV | sigz_deg@214 | sigz_deg@2856 | capture_in_bucket_pct | capture_all_buckets_pct |
| --- | --- | --- | --- | --- | --- | --- | --- |
| D after Sec1 | 3.3134 | 27.130 | 5.424 |  | 121.42 | 5.40 | 6.09 |
| Tan (published) |  | 27.200 | 3.500 |  | 11.10 | 89.4 | 96.8 |
| Δ (repo−Tan) |  | -0.070 | 1.924 |  | 110.319 |  |  |

## Capture bookkeeping (denominators)

- `q_injected_C` (honest, pre-iris) = 7.4584e-10 C
- `q_in_bore_C` (post-iris)         = 2.3810e-10 C
- `q_in_domain_C`                   = 2.4150e-10 C
- after-Sec1 in-bucket / q_in_bore  = 16.90% (separates iris loss from capture loss)

## Footnotes

- **σ_z internal ratio:** the repo @2856/@214 columns derive from ONE σ_t, so their ratio is BY CONSTRUCTION F_2856/F_214 = 13.335 — an internal-consistency check only. Tan's own ratio is 13.52 (a different spatial-σ_z definition with a per-location β); the 13.335 check is NOT agreement with Tan.
- **Denominator mismatch:** Tan's 100% (upstream) and 89.4/96.8 (postSec1) reference the gun-emitted bunch with no loss yet. The repo's `q_injected_C` already sits past the converging halo and includes the real 9.547 mm iris scrape. Repo capture = 5.4% in-bucket / 6.1% all-buckets (both vs `q_injected_C`); iris transmission = 31.9% (`q_in_bore`/`q_injected`) is a SEPARATE quantity (do NOT conflate it with capture). The repo numbers are far below Tan's BY CONSTRUCTION (γ² self-field + real iris + operating point). Annotate, do NOT 'fix'.
- **Location C capture = n/a (pre-iris):** C is the pre-scrape injector population; its charge differs from D's iris survivors, so the C→D ratio is not a clean repo capture fraction.
- **after-Sec1 σ_z@214 blank:** Tan leaves it blank; the captured-core σ_z is reported @2856 only.
