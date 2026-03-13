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

## Reference

Toon, O. B., McKay, C. P., Ackerman, T. P., & Santhanam, K. (1989).
Rapid calculation of radiative heating rates and photodissociation rates
in inhomogeneous multiple scattering atmospheres.
*Journal of Geophysical Research*, 94(D13), 16287-16301.
https://doi.org/10.1029/JD094iD13p16287
