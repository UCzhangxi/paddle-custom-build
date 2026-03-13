"""Corrected Toon-McKay 1989 Longwave (Thermal Emission) Solver

This module provides a corrected Python reference implementation of the
two-stream longwave (thermal emission) solver based on Toon et al. 1989,
intended to fix an error in the pyharp ``toon_mckay89_longwave_impl.h``
implementation.

Key fix
-------
The original pyharp longwave implementation applied a **delta-Eddington
scaling** to the optical properties (optical depth ``dtau``, single-scattering
albedo ``w0``, and asymmetry parameter ``g``) before the two-stream
calculation.  Delta-Eddington scaling is a correction for the **strongly
forward-peaked phase function** encountered in the *shortwave* (solar)
regime.  It must **not** be applied to longwave (thermal emission)
calculations because:

1. There is no direct solar beam whose forward scattering needs to be
   re-mapped onto a reduced optical depth.
2. Delta-Eddington artificially reduces ``dtau`` → ``(1 - w*g²) dtau``, making
   the atmosphere appear more transparent than it really is and producing
   systematically too-large outgoing fluxes / too-low temperatures.
3. PICASO (Batalha et al. 2019), the reference implementation also based on
   Toon 1989, does **not** apply delta-Eddington for the longwave path.

Reference
---------
Toon, O. B., McKay, C. P., Ackerman, T. P., & Santhanam, K. (1989).
Rapid calculation of radiative heating rates and photodissociation rates
in inhomogeneous multiple scattering atmospheres.
*Journal of Geophysical Research*, 94(D13), 16287–16301.
"""

import numpy as np


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def toon_mckay89_longwave(nlay, be, dtau_in, w0_in, g_in, a_surf, nmu=5):
    """Toon-McKay 1989 longwave two-stream solver (corrected).

    Parameters
    ----------
    nlay : int
        Number of atmospheric layers.
    be : array_like, shape (nlay+1,)
        Planck function at each *level* in consistent flux units.
        ``be[0]`` is the top-of-atmosphere (TOA) level;
        ``be[nlay]`` is the surface level.
    dtau_in : array_like, shape (nlay,)
        Layer optical depth.  **No delta-Eddington pre-scaling should be
        applied before calling this function.**
    w0_in : array_like, shape (nlay,)
        Single-scattering albedo per layer.
    g_in : array_like, shape (nlay,)
        Asymmetry parameter per layer.
    a_surf : float
        Surface albedo (diffuse reflectivity).  Surface emissivity = 1 − a_surf.
    nmu : int, optional
        Number of Gaussian quadrature angles (default 5; maximum 5).

    Returns
    -------
    flx_up : ndarray, shape (nlay+1,)
        Upward flux at each level (TOA at index 0, surface at index nlay).
    flx_dn : ndarray, shape (nlay+1,)
        Downward flux at each level.

    Notes
    -----
    All array quantities are indexed **top-to-bottom**:
    ``be[0]`` = TOA, ``be[nlay]`` = surface.
    """
    nmu = min(nmu, 5)
    nlev = nlay + 1
    twopi = 2.0 * np.pi
    ubari = 0.5  # hemispheric mean cosine angle (Table 1, Toon 1989)

    # Gaussian quadrature angles and weights for half-sphere integration
    # (Toon 1989, Table 1)
    _uarr  = np.array([0.0985350858, 0.3045357266, 0.5620251898,
                       0.8019865821, 0.9601901429])
    _wuarr = np.array([0.0157479145, 0.0739088701, 0.1463869871,
                       0.1671746381, 0.0967815902])
    uarr  = _uarr[:nmu]
    wuarr = _wuarr[:nmu]

    # ------------------------------------------------------------------
    # CORRECTED: use original optical properties — no delta-Eddington.
    #
    # The pyharp C++ implementation previously did:
    #   gsq   = G_IN(i) ** 2
    #   w0[i] = (1 - gsq) * W_IN(i) / (1 - W_IN(i) * gsq)  # delta-Edd
    #   dtau[i] = (1 - W_IN(i) * gsq) * DTAU_IN(i)          # delta-Edd
    #   hg[i]  = G_IN(i) / (1 + G_IN(i))                    # reduced g
    #
    # That is appropriate only for shortwave calculations.  For longwave
    # we use the raw inputs directly.
    # ------------------------------------------------------------------
    be   = np.asarray(be,    dtype=float)
    dtau = np.asarray(dtau_in, dtype=float).copy()
    w0   = np.asarray(w0_in,   dtype=float).copy()
    g    = np.asarray(g_in,    dtype=float).copy()   # use g directly (not g/(1+g))

    # ------------------------------------------------------------------
    # Two-stream coefficients — hemispheric mean (Toon 1989 Table 1)
    # For the hemispheric-mean approximation:
    #   g1 = 2 - w0*(1+g),   g2 = w0*(1-g)
    #   lambda = sqrt(g1^2 - g2^2) = 2*sqrt((1-w0)*(1-w0*g))
    #   gamma  = (g1 - lambda) / g2   [Eq 22]
    # We compute the equivalent quantities via the "alpha" intermediate:
    #   alpha = sqrt((1-w0)/(1-w0*g))
    #   lambda = alpha * (1-w0*g) / ubari  [= 2*sqrt((1-w0)*(1-w0*g))]
    #   gamma  = (1-alpha)/(1+alpha)
    # ------------------------------------------------------------------
    alp  = np.sqrt(np.clip((1.0 - w0) / np.where(1.0 - w0 * g != 0,
                                                   1.0 - w0 * g, 1e-30), 0, None))
    lam  = alp * (1.0 - w0 * g) / ubari        # Eq 21 (Toon 1989)
    gam  = (1.0 - alp) / (1.0 + alp)           # Eq 22 (Toon 1989)
    term = ubari / np.where(1.0 - w0 * g != 0, 1.0 - w0 * g, 1e-30)
    # term = 1/(g1+g2) = 0.5/(1-w0*g)  [used in C± particular solution, Eq 27]

    # ------------------------------------------------------------------
    # Planck source function: B(tau) = B0 + B1*tau  (Toon Eq 26)
    # B0 = Planck at top of layer,  B1 = dB/dtau  (positive going down
    # for a typical atmosphere where T increases with depth)
    # ------------------------------------------------------------------
    B0 = be[:-1].copy()
    B1 = np.where(dtau > 1.0e-6,
                  (be[1:] - be[:-1]) / dtau,
                  0.0)
    # For optically thin layers use the layer-mean Planck value
    B0 = np.where(dtau <= 1.0e-6, 0.5 * (be[1:] + be[:-1]), B0)

    # Particular-solution source functions at layer top (Cpm1/Cmm1) and
    # bottom (Cp/Cm)  — Toon Eq 27
    Cpm1 = B0 + B1 * term                    # C+ at layer top
    Cmm1 = B0 - B1 * term                    # C- at layer top
    Cp   = B0 + B1 * dtau + B1 * term        # C+ at layer bottom
    Cm   = B0 + B1 * dtau - B1 * term        # C- at layer bottom

    # Exponential terms  (Toon Eq 44)
    exptrm = np.minimum(lam * dtau, 35.0)
    Ep = np.exp(exptrm)          # exp(+lambda * dtau)
    Em = 1.0 / Ep                # exp(-lambda * dtau)
    E1 = Ep + gam * Em
    E2 = Ep - gam * Em
    E3 = gam * Ep + Em
    E4 = gam * Ep - Em

    # ------------------------------------------------------------------
    # Boundary values
    # ------------------------------------------------------------------
    Btop  = be[0]    # Planck at TOA level
    Bsurf = be[-1]   # Planck at surface level

    # ------------------------------------------------------------------
    # Build tridiagonal matrix  (Toon Eq 44 / setup_tri_diag in PICASO)
    # Unknowns: xkk[2n], xkk[2n+1]  for n = 0..nlay-1
    # The system is (2*nlay) × (2*nlay) with sub/main/super diagonals.
    # ------------------------------------------------------------------
    l  = 2 * nlay
    Af = np.zeros(l)
    Bf = np.zeros(l)
    Cf = np.zeros(l)
    Df = np.zeros(l)

    # Row 0 — TOA boundary: downward diffuse flux = Btop
    Af[0] = 0.0
    Bf[0] = gam[0] + 1.0
    Cf[0] = gam[0] - 1.0
    Df[0] = Btop - Cmm1[0]

    # Rows 1, 3, 5, … (odd) — layer interface continuity (flux-minus)
    for n in range(nlay - 1):
        i = 2 * n + 1
        Af[i] = (E1[n] + E3[n]) * (gam[n + 1] - 1.0)
        Bf[i] = (E2[n] + E4[n]) * (gam[n + 1] - 1.0)
        Cf[i] = 2.0 * (1.0 - gam[n + 1] ** 2)
        Df[i] = ((gam[n + 1] - 1.0) * (Cpm1[n + 1] - Cp[n]) +
                 (1.0 - gam[n + 1]) * (Cm[n] - Cmm1[n + 1]))

    # Rows 2, 4, 6, … (even, excluding 0) — layer interface continuity (flux-plus)
    for n in range(nlay - 1):
        i = 2 * n + 2
        Af[i] = 2.0 * (1.0 - gam[n] ** 2)
        Bf[i] = (E1[n] - E3[n]) * (1.0 + gam[n + 1])
        Cf[i] = (E1[n] + E3[n]) * (gam[n + 1] - 1.0)
        Df[i] = (E3[n] * (Cpm1[n + 1] - Cp[n]) +
                 E1[n] * (Cm[n] - Cmm1[n + 1]))

    # Row l-1 — surface boundary: upward flux = surface emission + reflection
    Af[l - 1] = E1[nlay - 1] - a_surf * E3[nlay - 1]
    Bf[l - 1] = E2[nlay - 1] - a_surf * E4[nlay - 1]
    Cf[l - 1] = 0.0
    Df[l - 1] = Bsurf - Cp[nlay - 1] + a_surf * Cm[nlay - 1]

    # Solve the tridiagonal system
    xkk = _dtridgl(l, Af, Bf, Cf, Df)

    # Extract Toon Table 3 coefficients  (Y1+Y2 and Y1-Y2)
    xk1 = xkk[0::2] + xkk[1::2]   # positive combination
    xk2 = xkk[0::2] - xkk[1::2]   # negative combination
    xk2[np.abs(xk2) < 1e-30 * np.abs(xkk[1::2])] = 0.0

    em1 = np.exp(-np.minimum(lam * dtau, 35.0))   # exp(-lambda*dtau)

    # ------------------------------------------------------------------
    # Source-function technique coefficients  (Toon Table 3)
    #
    # For the upwelling stream  (Eq 55):
    #   G / (lambda*mu - 1) * (Ep*em2 - 1)  +  H / (lambda*mu + 1) * (1 - Em*em2)
    #   + alpha1*(1-em2)  +  alpha2*(mu - (dtau+mu)*em2)
    #
    # For the downwelling stream (Eq 55):
    #   J / (lambda*mu + 1) * (Ep - em2)  +  K / (lambda*mu - 1) * (em2 - Em)
    #   + sigma1*(1-em2)  +  sigma2*(mu*em2 + dtau - mu)
    #
    # G, H, J, K are the scattering source contributions.
    # alpha1/alpha2 and sigma1/sigma2 are the thermal emission contributions.
    #
    # Mathematical note: pyharp uses the equivalent formulation
    #   G = 2*pi * w0 * xk1 * (1 + g*alp) / (1 + alp)
    # which is numerically identical to PICASO's
    #   G = (1/mu1 - lambda) * positive
    # because the two matrix systems have different RHS scalings (pyharp
    # omits the 2*pi*mu1 factor from C± while PICASO includes it, so the
    # solved coefficients differ by that factor, which is compensated by
    # the explicit 2*pi in pyharp's G/H/J/K expressions).
    # ------------------------------------------------------------------
    G      = np.zeros(nlay)
    H      = np.zeros(nlay)
    J      = np.zeros(nlay)
    K      = np.zeros(nlay)
    alpha1 = np.zeros(nlay)
    alpha2 = np.zeros(nlay)
    sigma1 = np.zeros(nlay)
    sigma2 = np.zeros(nlay)

    for n in range(nlay):
        if w0[n] <= 1.0e-4:
            # Pure-absorption limit — no scattering source terms
            alpha1[n] = twopi * B0[n]
            alpha2[n] = twopi * B1[n]
            sigma1[n] = alpha1[n]
            sigma2[n] = alpha2[n]
        else:
            den = 1.0 + alp[n]
            # Scattering source-function amplitudes
            G[n] = twopi * w0[n] * xk1[n] * (1.0 + g[n] * alp[n]) / den
            H[n] = twopi * w0[n] * xk2[n] * (1.0 - g[n] * alp[n]) / den
            J[n] = twopi * w0[n] * xk1[n] * (1.0 - g[n] * alp[n]) / den
            K[n] = twopi * w0[n] * xk2[n] * (1.0 + g[n] * alp[n]) / den
            # Planck source terms (corrected: use g directly, not g/(1+g))
            # term_val = ubari * w0 * g / (1 - w0*g)  [= 1/(g1+g2) - mu1]
            term_val   = ubari * w0[n] * g[n] / (1.0 - w0[n] * g[n])
            alpha1[n]  = twopi * (B0[n] + B1[n] * term_val)
            alpha2[n]  = twopi * B1[n]
            sigma1[n]  = twopi * (B0[n] - B1[n] * term_val)
            sigma2[n]  = alpha2[n]

    # ------------------------------------------------------------------
    # Gaussian quadrature — accumulate upward and downward fluxes
    # ------------------------------------------------------------------
    flx_up = np.zeros(nlev)
    flx_dn = np.zeros(nlev)

    for m in range(nmu):
        u = uarr[m]
        lw_dn = np.zeros(nlev)
        lw_up = np.zeros(nlev)

        # ---- Downward pass (TOA → surface) ----
        lw_dn[0] = twopi * Btop   # downward intensity entering at TOA
        for k in range(nlay):
            em2    = np.exp(-dtau[k] / u)
            l_u_p1 = lam[k] * u + 1.0
            l_u_m1 = lam[k] * u - 1.0

            lw_dn[k + 1] = (lw_dn[k] * em2
                            + (J[k] / l_u_p1) * (Ep[k] - em2)
                            + (K[k] / l_u_m1) * (em2 - em1[k])
                            + sigma1[k] * (1.0 - em2)
                            + sigma2[k] * (u * em2 + dtau[k] - u))

        # ---- Upward pass (surface → TOA) ----
        lw_up[nlev - 1] = twopi * (Bsurf + B1[nlay - 1] * u)  # surface emission
        for k in range(nlay - 1, -1, -1):
            em2    = np.exp(-dtau[k] / u)
            em3    = em1[k] * em2          # exp(-lambda*dtau) * exp(-dtau/u)
            l_u_m1 = lam[k] * u - 1.0
            l_u_p1 = lam[k] * u + 1.0

            lw_up[k] = (lw_up[k + 1] * em2
                        + (G[k] / l_u_m1) * (Ep[k] * em2 - 1.0)
                        + (H[k] / l_u_p1) * (1.0 - em3)
                        + alpha1[k] * (1.0 - em2)
                        + alpha2[k] * (u - (dtau[k] + u) * em2))

        flx_dn += lw_dn * wuarr[m]
        flx_up += lw_up * wuarr[m]

    return flx_up, flx_dn


# ---------------------------------------------------------------------------
# Internal helper
# ---------------------------------------------------------------------------

def _dtridgl(l, a, b, c, d):
    """Thomas algorithm for a tridiagonal system.

    Solves  a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1] = d[i].
    """
    x       = np.zeros(l)
    c_prime = np.zeros(l)
    d_prime = np.zeros(l)

    c_prime[0] = c[0] / b[0]
    d_prime[0] = d[0] / b[0]

    for i in range(1, l):
        denom = b[i] - a[i] * c_prime[i - 1]
        if i < l - 1:
            c_prime[i] = c[i] / denom
        d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) / denom

    x[l - 1] = d_prime[l - 1]
    for i in range(l - 2, -1, -1):
        x[i] = d_prime[i] - c_prime[i] * x[i + 1]

    return x
