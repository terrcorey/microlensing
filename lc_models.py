"""Functional forms used to model a microlensing light curve.

Point-source point-lens (PSPL / Paczynski 1986), plus a point-source
binary-lens (2L1S) extension. Further extensions (finite source, annual
parallax) can be added here later as their own functions once we actually
need them.
"""

import numpy as np
import torch
from astropy.coordinates import get_body_barycentric_posvel, SkyCoord
from astropy.time import Time
import astropy.units as u

ZERO_POINT_MAG = 18.0  # arbitrary; cancels out of every fitted flux ratio


def mag_to_flux(mag, mag_err):
    """Calibrated magnitude + error -> flux + error on the ZERO_POINT_MAG scale."""
    flux = 10 ** (-0.4 * (mag - ZERO_POINT_MAG))
    flux_err = flux * mag_err * (np.log(10) / 2.5)
    return flux, flux_err


def trajectory(t, t0, u0, tE):
    """Source-lens separation u(t), in units of the Einstein radius."""
    tau = (t - t0) / tE
    return np.sqrt(u0**2 + tau**2)


def magnification(u):
    """Paczynski point-lens magnification A(u)."""
    u2 = u**2
    return (u2 + 2) / (u * np.sqrt(u2 + 4))


def flux(t, t0, u0, tE, f_source, f_blend):
    """Observed flux: magnified source plus a constant blend flux."""
    A = magnification(trajectory(t, t0, u0, tE))
    return f_source * A + f_blend


def magnitude(t, t0, u0, tE, f_source, f_blend):
    """Observed magnitude, for overplotting against real photometry."""
    return ZERO_POINT_MAG - 2.5 * np.log10(flux(t, t0, u0, tE, f_source, f_blend))


def binary_trajectory(t, t0, u0, tE, alpha):
    """Source-lens separation for binary lenses on the lens plane, given in
    complex form ζ(t), in units of the Einstein radius."""
    tau = (t - t0) / tE
    return (tau + 1j * u0) * np.exp(1j * alpha)


def lens_position(s, q):
    """Returns the lens masses and positions from separation s and mass ratio q"""
    m1 = 1 / (1 + q)
    m2 = q / (1 + q)
    z1 = -q * s / (1 + q)
    z2 = s / (1 + q)
    return m1, m2, z1, z2


def caustic_curve(s, q, n_phi=600):
    """Caustic curve(s) in the source plane: the image, under the lens
    equation, of the critical curve where the lens map's Jacobian vanishes.

    The critical curve is m1/(z-z1)^2 + m2/(z-z2)^2 = e^(i*phi) for
    phi in [0, 2pi) -- a quartic in z for each phi, with no conjugate to
    clear (unlike the quintic lens equation), so its coefficients are built
    directly via polynomial multiplication (np.convolve) rather than a
    symbolic derivation.
    """
    m1, m2, z1, z2 = lens_position(s, q)
    c1_sq = np.convolve([1, -z1], [1, -z1])  # (z-z1)^2, degree 2
    c2_sq = np.convolve([1, -z2], [1, -z2])  # (z-z2)^2, degree 2
    quad_coeffs = np.convolve(c1_sq, c2_sq)  # (z-z1)^2 (z-z2)^2, degree 4

    base = np.zeros(5, dtype=complex)
    base[-3:] = m1 * c2_sq + m2 * c1_sq  # degree 2, right-aligned into the degree-4 slot

    phis = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)
    points = []
    for phi in phis:
        z_crit = np.roots(base - complex(np.exp(1j * phi)) * quad_coeffs)
        points.append(z_crit - m1 / np.conj(z_crit - z1) - m2 / np.conj(z_crit - z2))
    return np.concatenate(points)


def _quintic_coefficients(zeta, zetab, m1, z1, z2):
    """Coefficients (highest degree first) of the degree-5 polynomial in z
    equivalent to the binary lens equation
    zeta = z - m1/conj(z-z1) - m2/conj(z-z2), for fixed zeta/s/q.

    Mechanically derived (not hand-transcribed) via sympy from that lens
    equation -- see derive_binary_quintic.py to regenerate/audit.
    """
    x0 = zetab**2
    x1 = z1*z2
    x2 = z1*zetab
    x3 = z2*zetab
    x4 = m1*z1
    x5 = m1*z2
    x6 = x1*zeta
    x7 = x0*zeta
    x8 = z2**2
    x9 = x8*z1
    x10 = 2*x9
    x11 = 2*x0
    x12 = z1**2
    x13 = 2*x12*z2
    x14 = m1*x12
    x15 = z1*zeta
    x16 = z2**3
    x17 = z2*zeta
    x18 = z1**3
    x19 = x18*z2
    x20 = x0*x12
    x21 = x0*x8
    x22 = x18*zetab
    x23 = x16*zetab
    x24 = 2*x2
    x25 = 2*zetab
    x26 = x25*zeta
    x27 = x12*x8
    x28 = 2*x3
    x29 = x8*zeta
    x30 = 2*x7
    x31 = 2*x12
    x32 = 4*x1
    x33 = m1*x8
    x34 = x12*zetab
    x35 = 2*x15
    x36 = -x33 - 2*x34*zeta + x35*x8 - 4*x6*zetab
    x37 = x12*zeta
    x38 = 2*x37
    x39 = x12*x29
    x40 = 2*x2*zeta
    x41 = -x12*z2 + x17*x18
    x42 = x18*zeta
    x43 = m1**2
    x44 = x12*x43
    x45 = 2*x1
    x46 = x16*zeta
    x47 = 2*zeta
    return [
        x0 + x1 - x2 - x3,
        -x10 - x11*z1 - x11*z2 + 2*x12*zetab - x13 - x4 + x5 - x6 - x7 + 2*x8*zetab + 4*z1*z2*zetab + z1*zeta*zetab + z2*zeta*zetab - z2 + zetab,
        m1*x24 - m1*x28 + x0*x32 + x1 - 5*x12*x3 + x14 + x15 + x16*z1 + x17*x31 + x17 + x19 - 5*x2*x8 + x20 + x21 - x22 - x23 - x24 - x25*x29 - x26 + 4*x27 + x30*z1 + x30*z2 + x36 + x8,
        m1*x12*zeta - m1*x40 + 2*m1*x8*z1 + 2*m1*x8*zetab + 2*m1*z2*zeta*zetab - x12*x7 + 4*x12*x8*zetab + 5*x12*z2*zeta*zetab + x12*zetab - x14*x25 - x15*x16 - x16*x31 + 2*x16*z1*zetab + x16*zeta*zetab - 2*x18*x8 + 2*x18*z2*zetab + x18*zeta*zetab - 2*x20*z2 - 2*x21*z1 - x29 - x31*x5 - x32*x7 - x33*zeta - x38 - 4*x39 + x4 - x41 - x5 - 3*x6 - x7*x8 + 5*x8*z1*zeta*zetab - x8*zetab - x9 + 4*z1*zeta*zetab + 2*z2*zeta*zetab + z2 - zeta,
        -m1*x42 + m1*x45 + m1*x46 + x10*x7 + 3*x12*x17 - x12*x23 + x13*x7 + x14*x26 + x14*x28 - x14 + x16*x18 + x16*x38 - x16*x4 - x16*x40 + 2*x18*x29 + x18*x5 - x19 + x20*x8 - x22*x8 - x24*x33 + x24*x8 - x26*x33 - x27 - x28*x42 - 4*x29*x34 + x29*x4 + x35 + x36 - x37*x5 - x4*x47 + x42 - x43*x45 + x43*x8 + x44 - x45 + x47*x5,
        m1*x12*x16 + 2*m1*x12*zeta + m1*x18*z2*zeta - 2*m1*x6 + 2*m1*x8*z1*zeta*zetab + m1*x8*z1 + x12*x16*zeta*zetab - x12*x5 + 2*x12*z2*zeta*zetab - 2*x14*x3*zeta - x16*x42 - x18*x33 + x18*x8*zeta*zetab + x18*x8 - x27*x7 - x29*x43 - x34*x8 - x37 - x39 - x4*x46 - x41 + 2*x43*z1*z2*zeta - x44*zeta,
    ]


_TORCH_DEVICE = "cuda" if torch.cuda.is_available() else "cpu"


def _companion_eigvals(coeffs):
    """Roots of a batch of degree-5 polynomials via a batched companion-matrix
    eigensolve on `_TORCH_DEVICE` (GPU if available) -- same construction
    numpy.roots uses per-polynomial, done here for the whole batch at once.
    `coeffs`: (N, 6) complex128 tensor, highest degree first.
    """
    n = coeffs.shape[-1] - 1
    normalized = coeffs / coeffs[:, :1]
    companion = torch.zeros(coeffs.shape[0], n, n, dtype=coeffs.dtype, device=coeffs.device)
    companion[:, 0, :] = -normalized[:, 1:]
    companion[:, torch.arange(1, n), torch.arange(n - 1)] = 1.0
    return torch.linalg.eigvals(companion)


def _binary_images_batch(zeta, s, q, tol=1e-6):
    """Shared batched core for binary_images/binary_magnification.

    Solves the binary lens equation as a polynomial (see
    _quintic_coefficients) for every point in `zeta` at once via a single
    batched eigensolve, instead of looping a per-point np.roots call.
    Clearing the complex conjugate also introduces spurious extra roots --
    every candidate is checked against the original lens equation; `valid`
    marks the survivors (3 or 5 per point, depending on whether that point
    sits inside a caustic).

    Returns (roots, valid, (m1, m2, z1, z2)): roots and valid are (N, 5).
    """
    m1, m2, z1, z2 = lens_position(s, q)
    zetab = np.conj(zeta)
    coeffs = np.stack(_quintic_coefficients(zeta, zetab, m1, z1, z2), axis=-1)

    coeffs_t = torch.as_tensor(coeffs, dtype=torch.complex128, device=_TORCH_DEVICE)
    roots = _companion_eigvals(coeffs_t).cpu().numpy()

    residual = roots - m1 / (np.conj(roots) - z1) - m2 / (np.conj(roots) - z2) - zeta[:, None]
    valid = np.abs(residual) < tol
    return roots, valid, (m1, m2, z1, z2)


def binary_images(zeta, s, q, tol=1e-6):
    """True image positions for a binary lens at source position(s) `zeta`."""
    zeta_arr = np.atleast_1d(zeta).astype(complex)
    roots, valid, _ = _binary_images_batch(zeta_arr, s, q, tol)
    images = [roots[i][valid[i]] for i in range(len(zeta_arr))]
    return images if np.ndim(zeta) else images[0]


def binary_magnification(zeta, s, q):
    """Total binary-lens magnification at source position(s) `zeta`.

    Sums 1/|J| over every true image -- each image contributes to the
    total observed flux regardless of its parity, so this is a sum of the
    per-image magnifications' absolute values, not a signed sum.
    """
    zeta_arr = np.atleast_1d(zeta).astype(complex)
    roots, valid, (m1, m2, z1, z2) = _binary_images_batch(zeta_arr, s, q)
    jacobian = 1 - np.abs(m1 / (roots - z1) ** 2 + m2 / (roots - z2) ** 2) ** 2
    A = np.where(valid, 1 / np.abs(jacobian), 0.0).sum(axis=1)
    return A if np.ndim(zeta) else A[0]

def _project(vec, n_hat, e_hat):
    """Project a Cartesian position/velocity vector onto the sky's
    (N, E) tangent-plane basis at the target."""
    return (vec.x * n_hat[0] + vec.y * n_hat[1] + vec.z * n_hat[2],
            vec.x * e_hat[0] + vec.y * e_hat[1] + vec.z * e_hat[2])

def sun_earth_projection(time, ra_str, dec_str, t0_par):
    """Returns (delta_s_n, delta_s_e): Earth's sky-projected position
    relative to the Sun, in AU, with the constant-velocity part at
    t0_par subtracted out -- i.e. only the non-uniform (curvature)
    signal that's NOT already degenerate with (t0, u0, tE) remains.
    """
    time_hjd = Time(time + 2450000, format="jd")
    t0 = Time(t0_par + 2450000, format="jd")
        
    earth_pos, _ = get_body_barycentric_posvel("earth", time_hjd)
    sun_pos, _ = get_body_barycentric_posvel("sun", time_hjd)
    earth_t0pos, earth_t0vel = get_body_barycentric_posvel("earth", t0)
    sun_t0pos, sun_t0vel = get_body_barycentric_posvel("sun", t0)

    s = earth_pos - sun_pos
    s0 = earth_t0pos - sun_t0pos
    v0 = earth_t0vel - sun_t0vel

    coords = SkyCoord(ra_str + " " + dec_str)
    ra = coords.ra.radian
    dec = coords.dec.radian
    n_hat = (-np.sin(dec) * np.cos(ra), -np.sin(dec) * np.sin(ra), np.cos(dec))
    e_hat = (-np.sin(ra), np.cos(ra), 0)

    s_N, s_E = _project(s, n_hat, e_hat)
    s0_N, s0_E = _project(s0, n_hat, e_hat)
    v0_N, v0_E = _project(v0, n_hat, e_hat)

    dt_days = (time_hjd - t0).to(u.day)
    delta_sN = s_N - s0_N - dt_days * v0_N
    delta_sE = s_E - s0_E - dt_days * v0_E

    return delta_sN, delta_sE