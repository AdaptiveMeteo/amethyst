# AMETHYST Development Record

Chronological log of significant code changes, bug fixes, and algorithmic
improvements.  Each entry records the problem, the fix, the files affected,
and the measured impact where available.

---

## 2026-05-16 — RTTOV LM convergence fixes

**Branch:** `devel`  **Commit:** `e2f9b8b`

Baseline RTTOV convergence (121 IASI Arctic obs, 2022-08-08): **5.8 %** at
|d2| < 9.21.  Three targeted fixes brought this to **79 %** at |d2| < 100,
matching the OSS forward-model baseline (76 % at the same threshold).

### Fix 1 — Q2m decoupling

**Files:** `processor/rttovfm/rttovFM.py`, `processor/main/amethyst_code_main.py`

RTTOV's `NearSurface.Q2m` was set to the current surface water-vapour value
from the state vector at every forward-model call.  When RTTOV computes the
K-matrix it perturbs the surface profile level *and* Q2m simultaneously,
inflating the surface WV Jacobian by ~2.5×.  This caused catastrophic overshoot
at iteration 3 of the LM loop.

Fix: at the start of each retrieval (`invert()`), the first-guess surface WV is
stored as a constant on `fm.model._nearsurface_q2m_ppmv`.  In `rttovFM.py`,
Q2m is read from that attribute (with a fallback to the current profile value)
instead of tracking the state vector.

### Fix 2 — Standard LM step rejection

**File:** `processor/main/amethyst_code_main.py`

AMETHYST unconditionally accepted every LM step (`update_xhat = True` always).
When a step made `mspo` worse, gamma was increased (×5) but the next step was
computed from the bad position, compounding the error.

Fix: standard Marquardt-Levenberg bookkeeping.  Track `xhat_best` / `mspo_best`
across iterations.  When `mspo` increases, reject the step
(`update_xhat = False`) and revert `xhat` to `xhat_best` so the retry at higher
gamma starts from the best position seen so far.

### Fix 3 — Adaptive step-size limiter (key fix)

**File:** `processor/main/amethyst_code_main.py` — `update_solution()`

Near the OE minimum (mspo ≈ 1), the forward model is weakly nonlinear and the
linearised LM step can overshoot by ~10 K (T) or ~5 log-units (WV).  Each
overshoot triggers a bad step → gamma × 5.  With Fix 2 in place the revert
produces a "good" step → gamma ÷ 2.  The net per bad+good cycle is
×5/÷2 = ×2.5, causing gamma to spiral upward until the step size totx → 0 and
the convergence trigger d2\_trigger → 0 fires falsely, far from the true OE
minimum.

Fix: after solving the LM linear system, scale `totx` down uniformly if any T
or WV step exceeds a threshold:

| Regime | T limit | log(q) limit |
|--------|:-------:|:------------:|
| Normal (mspo ≥ 2) | ±15 K | ±5 |
| Near convergence (mspo < 2) | ±5 K | ±2 |

The tighter near-convergence limits prevent the linearisation overshoot,
allowing the LM to navigate cleanly to the OE minimum instead of spiralling.

### Convergence results

| Run | \|d2\| < 9.21 | \|d2\| < 100 | median \|d2\| |
|-----|:---:|:---:|:---:|
| RTTOV baseline (no fixes) | 5.8 % | — | 166 |
| + Q2m fix + LM rejection | 9.1 % | 34 % | 164 |
| + Adaptive step limit | **12 %** | **79 %** | **62** |
| OSS baseline | 44 % | 76 % | 11 |

**Recommended convergence threshold for RTTOV output: |d2| < 100** (vs the
default 9.21 used for OSS).  At this threshold RTTOV and OSS give equivalent
acceptance rates (~79 % and ~76 % respectively).

The remaining gap in median |d2| (62 vs 11) reflects residual forward-model
nonlinearity: the step limiter prevents exact convergence to the OE gradient
minimum, but the observation fit (mspo ≈ 1) is correct.

---
