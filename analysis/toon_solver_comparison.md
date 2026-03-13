# Toon-McKay 1989 Longwave Solver — Comparison: pyharp vs PICASO

## Overview

This document records the systematic comparison between the longwave (thermal
emission) two-stream Toon-McKay 1989 solver as implemented in
**pyharp** (`src/rtsolver/toon_mckay89_longwave_impl.h`, repo
[chengcli/pyharp](https://github.com/chengcli/pyharp)) and **PICASO**
(`picaso/fluxes.py` → `get_thermal_1d`, repo
[natashabatalha/picaso](https://github.com/natashabatalha/picaso)).
PICASO is a well-tested radiative-transfer code that correctly implements the
Toon et al. (1989) source-function technique.

---

## Primary Bug: Delta-Eddington Scaling Applied to Longwave

### What is delta-Eddington scaling?

Delta-Eddington (or "delta" truncation) is a pre-processing step that
replaces the actual optical properties with scaled equivalents that better
represent a strongly forward-peaked scattering phase function:

```
w0_new  = (1 - g²) * w0  / (1 - w0 * g²)
dtau_new = (1 - w0 * g²) * dtau
g_new   = g / (1 + g)
```

This scaling is designed for the **shortwave** (solar) regime, where
radiation encounters highly anisotropic cloud/aerosol particles with a
dominant forward-scattering peak.

### What pyharp does (incorrect for longwave)

In `toon_mckay89_longwave_impl.h`, the very first operation in the solver is:

```cpp
// Delta-Eddington Scaling  ← APPLIED TO LONGWAVE (WRONG)
for (int i = 0; i < nlay; i++) {
    T gsq = G_IN(i) * G_IN(i);
    w0[i]   = (1.0 - gsq) * W_IN(i) / (1.0 - W_IN(i) * gsq);
    dtau[i] = (1.0 - W_IN(i) * gsq) * DTAU_IN(i);
    hg[i]   = G_IN(i) / (1.0 + G_IN(i));   // reduced asymmetry
}
```

### What PICASO does (correct)

In `get_thermal_1d`, PICASO uses the **original** optical properties directly:

```python
g1 = 2.0 - w0 * (1 + cosb)   # hemispheric mean, no delta-Eddington
g2 = w0 * (1 - cosb)
alpha = sqrt((1. - w0) / (1. - w0 * cosb))
lamda = sqrt(g1**2 - g2**2)   # Eq 21 Toon 1989
```

No delta-Eddington transformation is performed.

### Why this matters

| Effect | Delta-Eddington applied | No delta-Eddington (correct) |
|--------|------------------------|------------------------------|
| Optical depth | `dtau_new < dtau` | `dtau` unchanged |
| Effective w0  | Modified | Unchanged |
| Atmosphere apparent opacity | **Too transparent** | Correct |
| Outgoing thermal flux | **Too large** | Correct |
| Inferred temperatures | **Too cold** | Correct |

For a typical infrared cloud layer with `w0 = 0.9`, `g = 0.85`:
- `dtau_new = (1 − 0.9 × 0.7225) × dtau = 0.350 × dtau`
  (optical depth reduced to **35 %** of its true value!)
- The atmosphere behaves as if it is nearly transparent, letting far too much
  thermal radiation escape.

---

## Secondary Difference: Planck Source Function Terms

Once the delta-Eddington bug is fixed (i.e., original `g` is used directly),
the formulae for `alpha1`, `sigma1`, `term_val` in pyharp become
**mathematically equivalent** to the corresponding expressions in PICASO.

**pyharp** (corrected, using `g` directly):
```cpp
term_val = ubari * w0 * g / (1.0 - w0 * g)
alpha1   = 2π * (B0 + B1 * term_val)
sigma1   = 2π * (B0 - B1 * term_val)
```

**PICASO**:
```python
g1_plus_g2 = 1.0 / (g1 + g2)     # = 0.5 / (1 - w0*g)
alpha1 = 2*pi * (b0 + b1 * (g1_plus_g2 - mu1))
                                   # = 2π*(B0 + B1 * 0.5*w0*g/(1-w0*g))
```

Since `0.5 * w0 * g / (1 − w0 * g) = term_val` when using `mu1 = ubari = 0.5`,
both expressions are identical ✓.

---

## Secondary Difference: Scattering Source Coefficients G, H, J, K

Pyharp computes:
```cpp
G = 2π * w0 * xk1 * (1 + g*alp) / (1 + alp)   // "g" in pyharp code
H = 2π * w0 * xk2 * (1 - g*alp) / (1 + alp)
J = 2π * w0 * xk1 * (1 - g*alp) / (1 + alp)
K = 2π * w0 * xk2 * (1 + g*alp) / (1 + alp)
```

PICASO follows Toon Table 3 directly:
```python
G = (1/mu1 - lamda) * positive      # = (2 - λ) * positive
H = gama * (lamda + 1/mu1) * negative
J = gama * (lamda + 1/mu1) * positive
K = (1/mu1 - lamda) * negative
```

These look different, but they are **numerically equal** once the different
matrix-RHS normalizations are accounted for.  PICASO's `C±` terms include a
`2π·μ₁ = π` prefactor, so `positive ≈ π · xk1` from pyharp's solve.
Substituting:

```
pyharp G = 2π · w0 · xk1 · (1 + g·α)/(1 + α)
         = π · (2 - λ) · xk1          [since w0*(1+g*α)/(1+α) = (2-λ)/2]
         = (2 - λ) · positive = PICASO G  ✓
```

Numerical verification (w0 = 0.5, g = 0.5):

| Quantity | Value |
|----------|-------|
| α = √((1−w0)/(1−w0g)) | 0.8165 |
| λ = 2α(1−w0g) | 1.2247 |
| PICASO: (2−λ)/2 | 0.3877 |
| pyharp: w0(1+g·α)/(1+α) | 0.3876 ✓ |

---

## Minor Difference: TOA Downward Boundary Condition

**PICASO** corrects for the finite optical depth above the model top:
```python
tau_top = dtau[0] * plevel[0] / (plevel[1] - plevel[0])
b_top   = (1 - exp(-tau_top / mu1)) * B_TOA * pi
```
This gives a near-zero downward flux at the TOA for small `tau_top`.

**pyharp** uses:
```cpp
lw_down_g[0] = 2π * B(T_TOA)
```
This is non-zero and approximately equals PICASO when `tau_top → ∞`, but
exceeds PICASO for thin upper atmospheres.  For typical model configurations
with a cold, optically thin TOA, the difference is small.

---

## Minor Difference: Surface Upwelling Boundary

**PICASO** (non-hard surface) adds a first-order Planck gradient correction:
```python
b_surface = (B_surf + b1[-1] * mu1) * pi
```

**pyharp** uses only the surface Planck function:
```cpp
bsurf_flux = B(T_surface)
```

For most applications with coarse vertical resolution near the surface, the
`b1 * mu1` term is negligible.  For high-resolution near-surface layers it
could matter.

---

## Required Fix in pyharp

In `src/rtsolver/toon_mckay89_longwave_impl.h`, replace the
delta-Eddington scaling block:

```cpp
// Delta-Eddington Scaling  ← REMOVE THIS FOR LONGWAVE
for (int i = 0; i < nlay; i++) {
    T gsq = G_IN(i) * G_IN(i);
    w0[i]   = (1.0 - gsq) * W_IN(i) / (1.0 - W_IN(i) * gsq);
    dtau[i] = (1.0 - W_IN(i) * gsq) * DTAU_IN(i);
    hg[i]   = G_IN(i) / (1.0 + G_IN(i));
}
```

with direct assignment of the original optical properties:

```cpp
// No delta-Eddington scaling for longwave thermal emission
for (int i = 0; i < nlay; i++) {
    w0[i]   = W_IN(i);
    dtau[i] = DTAU_IN(i);
    hg[i]   = G_IN(i);   // use g directly, not g/(1+g)
}
```

All downstream formulae that use `hg` remain correct because they compute
the same dimensionless ratios — they just now operate on the original
asymmetry parameter `g` rather than the delta-Eddington-reduced `g/(1+g)`.

A reference Python implementation of the corrected solver is provided in
`src/longwave_toon_corrected.py`.

---

## Summary Table

| Issue | pyharp (original) | PICASO | Fix |
|-------|-------------------|--------|-----|
| Delta-Eddington scaling | ✗ Applied to longwave | ✓ Not applied | Remove scaling in LW init loop |
| Two-stream λ, γ formulas | Correct form, wrong inputs | Correct | Use original g, w0 |
| G/H/J/K scattering terms | Numerically equivalent | Reference | No change needed |
| α₁/σ₁ Planck terms | Correct form, wrong inputs | Correct | Use original g, w0 |
| TOA downward BC | 2π·B(T_TOA) | tau_top corrected, ≈0 | Secondary; pyharp is conservative |
| Surface BC | B(T_surf) only | (B_surf + B1·μ₁)·π | Minor; affects thin near-sfc layers |

---

## Shortwave Solver Check (pyharp vs PICASO)

This section extends the same methodology to the shortwave (solar/reflection)
Toon solver in:

- **pyharp**: `src/rtsolver/toon_mckay89_shortwave_impl.h`
- **PICASO**: `picaso/fluxes.py` (`get_reflected_1d`) + `picaso/optics.py`

### What is consistent

For the diffuse two-stream core (quadrature closure), pyharp is largely
consistent with PICASO:

- delta-Eddington optical-property scaling is present
- two-stream coefficients use the same quadrature-form Toon Table 1 equations
- tridiagonal matrix assembly follows the same rotated-layer form

### Inconsistency 1 (important): zero-scattering branch drops surface reflection

In pyharp, when all `w0 == 0`, the code takes a direct-beam-only shortcut:

```cpp
FLX_DN(k) = F0_in * mu * exp(-tau/mu);
FLX_DN(nlev - 1) *= (1.0 - w_surf_in);
for (int k = 0; k < nlev; k++) FLX_UP(k) = 0.0;
```

This removes reflected energy from downward flux at the bottom level but does
not put it into `FLX_UP`. So energy reflected by the surface disappears,
creating an energy-conservation error that scales with surface albedo and
incident direct-beam flux.

PICASO keeps the surface reflection as an explicit lower boundary source:

```python
b_surface = surf_reflect*u0*F0PI*exp(-tau[-1,:]/u0)
```

which then generates upward diffuse flux through the matrix solve.

**Suggested fix in pyharp**:
- in the zero-scattering shortcut, keep `FLX_DN` as the incoming direct beam
- add reflected component to `FLX_UP` at the surface (or route through the same
  boundary-condition machinery as the general branch)
- do not absorb reflected energy by scaling only `FLX_DN(nlev-1)`

### Inconsistency 2 (modeling choice): corrected optical path used exclusively for this solver path

PICASO carries **both** delta-Eddington-corrected and original optical
properties (`dtau/tau/w0/cosb` and `dtau_og/tau_og/w0_og/cosb_og`). It uses the
uncorrected set for single-scattering source terms in reflected-light intensity
calculation.

pyharp shortwave currently applies delta-Eddington up front and then uses the
scaled properties everywhere in this solver path.

This is not always wrong for diffuse fluxes, but it is less flexible than
PICASO’s split treatment and can bias direct/single-scattering contributions in
forward-peaked cases.

**Suggested improvement in pyharp**:
- keep current scaled path for diffuse multiple-scattering solve
- optionally carry original (`*_og`) optical properties for direct-beam /
  single-scattering terms, following PICASO’s separation

### Inconsistency 3 (minor): no closure toggle

PICASO exposes both quadrature and Eddington Toon coefficients
(`toon_coefficients`), while pyharp hardcodes the quadrature form.

This is not a correctness bug (quadrature is PICASO default), but adding an
optional closure mode would improve cross-model reproducibility for sensitivity
tests.

### Shortwave comparison summary

| Item | pyharp shortwave | PICASO shortwave | Assessment |
|------|------------------|------------------|------------|
| Diffuse two-stream core (quadrature) | Implemented | Implemented | Consistent |
| Delta-Eddington support | Implemented | Implemented | Consistent |
| Zero-scattering + reflective surface handling | Reflected energy dropped in shortcut branch | Reflection enters lower BC source | **Inconsistent (bug)** |
| Uncorrected optical path for direct/single terms | Not separated in this solver path | Explicit `*_og` path available | Potential inconsistency |
| Closure options (quadrature/Eddington) | Quadrature only | Toggle available | Minor capability gap |

---

## Minimal Necessary pyharp Changes for Full PICASO Parity (Longwave + Shortwave)

This is the **minimal** safe change set if the target is numerical consistency
with PICASO Toon solvers, including boundary treatment.

### A. What was verified outside PICASO solver kernels (important)

To avoid copying behavior that actually comes from upstream preprocessing, the
following PICASO pipeline pieces were checked:

- `picaso/optics.py::compute_opacity`
  - computes both delta-Eddington-corrected (`dtau, tau, w0, cosb`) and
    original (`dtau_og, tau_og, w0_og, cosb_og`) optical properties.
- `picaso/justdoit.py`
  - passes corrected properties to reflected-light two-stream core;
  - passes original/no-raman properties to thermal solver;
  - passes both corrected and original sets to reflected solver because
    single-scattering source terms use original properties.

**Implication for pyharp changes:**  
For solver-level flux parity, you do **not** need to replicate all PICASO
upstream chemistry/opacity machinery. You only need the same optical-property
choice at each solver call and the same boundary conditions.

### B. Minimal pyharp code changes

#### 1) Longwave (`toon_mckay89_longwave_impl.h`)

1. Keep/remove delta-Eddington exactly as already concluded:
   - no delta-Eddington in longwave path;
   - use original `w0`, `dtau`, `g`.
2. Add TOA finite-optical-depth boundary option:
   - replace hard `Btop` usage with
     `Btop_eff = (1 - exp(-tau_top / ubari)) * Btop`,
     where `tau_top = dtau[0] * p_top / (p_1 - p_top)`.
   - if pressure metadata is unavailable, keep `tau_top=0` default path
     (compatibility mode), but expose input for parity mode.
3. Add PICASO-consistent lower thermal boundary:
   - gas-giant/non-hard-surface mode:
     `Bsurf_eff = Bsurf + B1_last * ubari`.
   - hard-surface mode:
     `Bsurf_eff = (1 - a_surf) * Bsurf`.
   - use `Bsurf_eff` in the last-row RHS.

These three items are sufficient for longwave parity at solver level.

#### 2) Shortwave (`toon_mckay89_shortwave_impl.h`)

1. Fix `all_zero_w` shortcut branch to conserve reflected energy:
   - do **not** attenuate `FLX_DN(nlev-1)` by `(1 - w_surf_in)` as a sink;
   - set `FLX_UP(nlev-1)` to reflected direct beam:
     `w_surf_in * FLX_DN(nlev-1)` (using incoming direct beam before any sink);
   - keep `FLX_UP(k)=0` for upper levels in this no-scattering branch.
2. Keep existing general two-stream branch unchanged (already aligned in form).

This is the only shortwave bug fix required for boundary-consistent parity in
the current pyharp flux-style implementation.

### C. Exact patch intent (pyharp-side) in compact form

#### Longwave BC intent

```cpp
// New optional inputs: p_top, p_next, hard_surface (or tau_top directly)
T tau_top = /* from pressure geometry, else 0 */;
T Btop_eff = (1.0 - exp(-tau_top / ubari)) * Btop;
Df[0] = Btop_eff - Cmm1[0];

T Bsurf_eff;
if (hard_surface) {
  Bsurf_eff = (1.0 - a_surf) * Bsurf;
} else {
  Bsurf_eff = Bsurf + B1_last * ubari;
}
Df[l-1] = Bsurf_eff - Cp[nlay-1] + a_surf * Cm[nlay-1];
```

#### Shortwave zero-scattering branch intent

```cpp
// Keep direct beam as downwelling
for (int k = 0; k < nlev; ++k) FLX_DN(k) = Fdir(k);
for (int k = 0; k < nlev; ++k) FLX_UP(k) = 0.0;
FLX_UP(nlev - 1) = w_surf_in * FLX_DN(nlev - 1);
```

### D. Careful validation plan (what to test before/after pyharp edits)

Run these as regression cases in pyharp after patching:

1. **Longwave, pure absorption, thin TOA**
   - expect top downwelling flux to follow finite-`tau_top` boundary
     (near zero when `tau_top << 1`).
2. **Longwave, non-hard-surface**
   - verify bottom upwelling differs by expected `B1_last * ubari` term.
3. **Longwave, hard-surface**
   - verify emissivity scaling `(1-a_surf)` at lower boundary.
4. **Shortwave, all-zero-SSA with reflective surface**
   - confirm energy conservation:
     reflected power appears in `FLX_UP(surface)` and is no longer lost.
5. **Shortwave, scattering atmosphere**
   - ensure no regressions in general branch outputs.

Recommended numerical tolerances against PICASO for matched inputs:

- layer/interface fluxes: relative error `<= 1e-5` (or absolute `<= 1e-8`
  where fluxes are tiny),
- net flux residual (energy closure): `<= 1e-6` of incident flux scale.

### E. What is intentionally *not* required for minimal parity

- Reproducing PICASO chemistry/opacity table generation internals.
- Reproducing Raman/correlated-k plumbing outside solver calls.
- Adding Eddington-closure toggle (useful, but not required for default-parity
  target because PICASO default reflected Toon path is quadrature).

---

## Reference

Toon, O. B., McKay, C. P., Ackerman, T. P., & Santhanam, K. (1989).
Rapid calculation of radiative heating rates and photodissociation rates
in inhomogeneous multiple scattering atmospheres.
*Journal of Geophysical Research*, 94(D13), 16287-16301.
https://doi.org/10.1029/JD094iD13p16287
