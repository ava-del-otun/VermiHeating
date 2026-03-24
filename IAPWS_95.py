# IAPWS_95.py
"""
IAPWS-95 ancillary routines used in this project.

Implemented here (SI outputs):
- Dynamic viscosity mu(rho, T) -> Pa·s  (IAPWS 2008 viscosity formulation)
- Thermal conductivity k(rho, T) -> W/m·K (IAPWS 2011 thermal conductivity formulation)
- Helmholtz core derivatives needed for cp/cv and critical enhancement terms
- Saturation-pressure ancillary p_sat(T) -> Pa
- Saturated liquid/vapor density helpers from the Wagner-Pruss auxiliary saturation fits
- Saturated liquid/vapor enthalpy helpers evaluated at those auxiliary saturation densities
- Latent heat along saturation from h_v,sat(T) - h_l,sat(T)

Notes on units:
- rho: kg/m^3
- T: K
- Pc: Pa
- R: J/kg/K
- viscosity intermediate mu0 is in microPa·s and converted to Pa·s at the end.
- thermal conductivity lambda0 and lambda2 intermediates are in mW/m·K and converted to W/m·K.

This file intentionally keeps the algebra close to the MATLAB script you provided to
minimize transcription risk.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import numpy as np

# -----------------------
# Fluid constants (water)
# -----------------------
M_W = 18.015268e-3          # kg/mol
TcW = 647.096               # K
rho_cW = 322.0              # kg/m^3
PcW = 22.064e6              # Pa
R_95 = 461.51805            # J/kg/K

# -----------------------
# Coefficients (as given)
# -----------------------
n0w = np.array([-8.32044648201, 6.6832105268, 3.00632, 0.012436, 0.97315, 1.27950, 0.96956, 0.24873], dtype=float)
ga0w = np.array([1, 1, 1, 1.28728967, 3.53734222, 7.74073708, 9.24437796, 27.5075105], dtype=float)

nw = np.array([
    0.12533547935523e-1, 0.78957634722828e1, -0.87803203303561e1, 0.31802509345418, -0.26145533859358,
    -0.78199751687981e-2, 0.88089493102134e-2,
    -0.66856572307965, 0.20433810950965, -0.66212605039687e-4, -0.19232721156002, -0.25709043003438,
    0.16074868486251, -0.40092828925807e-1, 0.39343422603254e-6, -0.75941377088144e-5,
    0.56250979351888e-3, -0.15608652257135e-4, 0.11537996422951e-8, 0.36582165144204e-6,
    -0.13251180074668e-11, -0.62639586912454e-9,
    -0.10793600908932, 0.17611491008752e-1, 0.22132295167546, -0.40247669763528,
    0.58083399985759, 0.49969146990806e-2, -0.31358700712549e-1, -0.74315929710341,
    0.47807329915480, 0.20527940895948e-1, -0.13636435110343, 0.14180634400617e-1,
    0.83326504880713e-2, -0.29052336009585e-1, 0.38615085574206e-1, -0.20393486513704e-1,
    -0.16554050063734e-2, 0.19955571979541e-2, 0.15870308324157e-3, -0.16388568342530e-4,
    0.43613615723811e-1, 0.34994005463765e-1, -0.76788197844621e-1, 0.22446277332006e-1,
    -0.62689710414685e-4, -0.55711118565645e-9,
    -0.19905718354408, 0.31777497330738, -0.11841182425981,
    -0.31306260323435e2, 0.31546140237781e2, -0.25213154341695e4,
    -0.14874640856724, 0.31806110878444
], dtype=float)

tw = np.array([
    -0.5, 0.875, 1, 0.5, 0.75, 0.375, 1, 4, 6, 12, 1, 5, 4, 2, 13, 9, 3, 4, 11, 4,
    13, 1, 7, 1, 9, 10, 10, 3, 7, 10, 10, 6, 10, 10, 1, 2, 3, 4, 8, 6, 9, 8,
    16, 22, 23, 23, 10, 50, 44, 46, 50, 0, 1, 4
], dtype=float)

dw = np.array([
    1, 1, 1, 2, 2, 3, 4, 1, 1, 1, 2, 2, 3, 4, 4, 5, 7, 9, 10, 11, 13, 13, 1, 2,
    2, 2, 3, 4, 4, 4, 5, 6, 6, 7, 9, 9, 9, 9, 9, 10, 10, 12, 3, 4, 4, 5, 14,
    3, 6, 6, 6, 3, 3, 3
], dtype=float)

cw = np.array([
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    3, 3, 3, 3,
    4, 6, 6, 6, 6
], dtype=float)

bw = np.array([*([1.0]*54), 0.85, 0.95], dtype=float)
aw = np.array([*([1.0]*54), 3.5, 3.5], dtype=float)
Aw = np.array([*([1.0]*54), 0.32, 0.32], dtype=float)
Bw = np.array([*([1.0]*54), 0.2, 0.2], dtype=float)
Cw = np.array([*([1.0]*54), 28.0, 32.0], dtype=float)
Dw = np.array([*([1.0]*54), 700.0, 800.0], dtype=float)

alphw = np.array([*([1.0]*51), 20.0, 20.0, 20.0], dtype=float)       # length 54
betw  = np.array([*([1.0]*51), 150.0, 150.0, 250.0, 0.3, 0.3], dtype=float)  # length 56
gaw   = np.array([*([1.0]*51), 1.21, 1.21, 1.25], dtype=float)       # length 54
epsw  = np.array([*([1.0]*54)], dtype=float)                         # length 54


# -----------------------
# Core helpers
# -----------------------
def helmholtz_all(tau: float, delta: float):
    """
    Returns:
        phi0, phi0_t, phi0_tt,
        phiR, phiR_d, phiR_t, phiR_dd, phiR_tt, phiR_dt
    All derivatives are with respect to tau and/or delta in the IAPWS-95 dimensionless Helmholtz form.
    """
    tau = float(tau)
    delta = float(delta)

    # ---- Ideal part ----
    # MATLAB: sumPhiE0 over i=4:8 (1-based) -> python indices 3..7
    sumPhiE0 = 0.0
    for i in range(3, 8):
        sumPhiE0 += n0w[i] * np.log(1.0 - np.exp(-ga0w[i] * tau))
    phi0 = np.log(delta) + n0w[0] + n0w[1] * tau + n0w[2] * np.log(tau) + sumPhiE0

    sumPhi0_t = 0.0
    for i in range(3, 8):
        sumPhi0_t += n0w[i] * ga0w[i] * ((1.0 - np.exp(-ga0w[i] * tau)) ** (-1) - 1.0)
    phi0_t = sumPhi0_t + n0w[1] + n0w[2] / tau

    sumPhi0_tt = 0.0
    for i in range(3, 8):
        sumPhi0_tt += n0w[i] * (ga0w[i] ** 2) * np.exp(-ga0w[i] * tau) * (1.0 - np.exp(-ga0w[i] * tau)) ** (-2)
    phi0_tt = -sumPhi0_tt - n0w[2] / (tau ** 2)

    # ---- Residual part (split into 4 groups) ----
    sum1 = sum2 = sum3 = sum4 = 0.0
    sum1_d = sum1_t = sum1_dd = sum1_tt = sum1_dt = 0.0
    sum2_d = sum2_t = sum2_dd = sum2_tt = sum2_dt = 0.0
    sum3_d = sum3_t = sum3_dd = sum3_tt = sum3_dt = 0.0
    sum4_d = sum4_t = sum4_dd = sum4_tt = sum4_dt = 0.0

    # Group 1: i=1:7 (MATLAB) -> indices 0..6
    for idx in range(0, 7):
        n = nw[idx]
        t = tw[idx]
        d = dw[idx]
        sum1 += n * (delta ** d) * (tau ** t)
        sum1_d += n * d * (delta ** (d - 1.0)) * (tau ** t)
        sum1_t += n * t * (delta ** d) * (tau ** (t - 1.0))
        sum1_dd += n * d * (d - 1.0) * (delta ** (d - 2.0)) * (tau ** t)
        sum1_tt += n * t * (t - 1.0) * (delta ** d) * (tau ** (t - 2.0))
        sum1_dt += n * d * t * (delta ** (d - 1.0)) * (tau ** (t - 1.0))

    # Group 2: i=8:51 -> indices 7..50 (uses cw)
    for idx in range(7, 51):
        n = nw[idx]
        t = tw[idx]
        d = dw[idx]
        c = cw[idx]
        e = np.exp(-(delta ** c))
        sum2 += n * (delta ** d) * (tau ** t) * e
        A = (d - c * (delta ** c))
        sum2_d += n * e * (delta ** (d - 1.0)) * (tau ** t) * A
        sum2_t += n * t * (delta ** d) * (tau ** (t - 1.0)) * e

        sum2_dd += n * e * (tau ** t) * (delta ** (d - 2.0)) * (
            A * (d - 1.0 - c * (delta ** c)) - (c ** 2) * (delta ** c)
        )
        sum2_tt += n * t * (t - 1.0) * (delta ** d) * (tau ** (t - 2.0)) * e
        sum2_dt += n * t * (delta ** (d - 1.0)) * (tau ** (t - 1.0)) * A * e

    # Group 3: i=52:54 -> indices 51..53
    for idx in range(51, 54):
        n = nw[idx]
        t = tw[idx]
        d = dw[idx]
        e = np.exp(-alphw[idx] * (delta - epsw[idx]) ** 2 - betw[idx] * (tau - gaw[idx]) ** 2)
        sum3 += n * (delta ** d) * (tau ** t) * e
        sum3_d += n * (delta ** d) * (tau ** t) * e * ((d / delta) - 2.0 * alphw[idx] * (delta - epsw[idx]))
        sum3_t += n * (delta ** d) * (tau ** t) * e * ((t / tau) - 2.0 * betw[idx] * (tau - gaw[idx]))
        sum3_dd += n * (delta ** d) * (tau ** t) * e * (
            ((d / delta) - 2.0 * alphw[idx] * (delta - epsw[idx])) ** 2 - d / (delta ** 2) - 2.0 * alphw[idx]
        )
        sum3_tt += n * (delta ** d) * (tau ** t) * e * (
            ((t / tau) - 2.0 * betw[idx] * (tau - gaw[idx])) ** 2 - t / (tau ** 2) - 2.0 * betw[idx]
        )
        sum3_dt += n * (delta ** d) * (tau ** t) * e * (
            ((d / delta) - 2.0 * alphw[idx] * (delta - epsw[idx])) * ((t / tau) - 2.0 * betw[idx] * (tau - gaw[idx]))
        )

    # Group 4: i=55:56 -> indices 54..55
    for idx in range(54, 56):
        n = nw[idx]
        A = Aw[idx]
        B = Bw[idx]
        C = Cw[idx]
        D = Dw[idx]
        a = aw[idx]
        b = bw[idx]
        beta = betw[idx]

        dm1 = delta - 1.0
        dm12 = dm1 ** 2

        Theta = (1.0 - tau) + A * (dm12) ** (1.0 / (2.0 * beta))
        Delta = Theta ** 2 + B * (dm12) ** a

        if abs(dm1) < 1e-14:
            Delta_d = 0.0
            Delta_dd = 0.0
        else:
            Delta_d = dm1 * (
                A * Theta * (2.0 / beta) * (dm12) ** (1.0 / (2.0 * beta) - 1.0)
                + 2.0 * B * a * (dm12) ** (a - 1.0)
            )
            Delta_dd = (Delta_d / dm1) + dm12 * (
                4.0 * B * a * (a - 1.0) * (dm12) ** (a - 2.0)
                + 2.0 * (A ** 2) * (1.0 / beta) ** 2 * ((dm12) ** (1.0 / (2.0 * beta) - 1.0)) ** 2
                + A * Theta * (4.0 / beta) * ((1.0 / (2.0 * beta)) - 1.0) * (dm12) ** (1.0 / (2.0 * beta) - 2.0)
            )

        Deltab_d = b * (Delta ** (b - 1.0)) * Delta_d
        Deltab_dd = b * ((Delta ** (b - 1.0)) * Delta_dd + (b - 1.0) * (Delta ** (b - 2.0)) * (Delta_d ** 2))
        Deltab_t = -2.0 * Theta * b * (Delta ** (b - 1.0))
        Deltab_tt = 2.0 * b * (Delta ** (b - 1.0)) + 4.0 * (Theta ** 2) * b * (b - 1.0) * (Delta ** (b - 2.0))
        Deltab_dt = (
            -A * b * (2.0 / beta) * (Delta ** (b - 1.0)) * dm1 * (dm12) ** (1.0 / (2.0 * beta) - 1.0)
            - 2.0 * Theta * b * (b - 1.0) * (Delta ** (b - 2.0)) * Delta_d
        )

        Psi = np.exp(-C * dm12 - D * ((tau - 1.0) ** 2))
        Psi_d = -2.0 * C * dm1 * Psi
        Psi_dd = (2.0 * C * (dm1 ** 2) - 1.0) * 2.0 * C * Psi
        Psi_t = -2.0 * D * (tau - 1.0) * Psi
        Psi_tt = (2.0 * D * ((tau - 1.0) ** 2) - 1.0) * 2.0 * D * Psi
        Psi_dt = 4.0 * C * D * dm1 * (tau - 1.0) * Psi

        sum4 += n * (Delta ** b) * delta * Psi
        sum4_d += n * ((Delta ** b) * (Psi + delta * Psi_d) + Deltab_d * delta * Psi)
        sum4_t += n * delta * (Deltab_t * Psi + (Delta ** b) * Psi_t)
        sum4_dd += n * (
            (Delta ** b) * (2.0 * Psi_d + delta * Psi_dd)
            + 2.0 * Deltab_d * (Psi + delta * Psi_d)
            + Deltab_dd * delta * Psi
        )
        sum4_tt += n * delta * (Deltab_tt * Psi + 2.0 * Deltab_t * Psi_t + (Delta ** b) * Psi_tt)
        sum4_dt += n * (
            (Delta ** b) * (Psi_t + delta * Psi_dt)
            + delta * Deltab_d * Psi_t
            + Deltab_t * (Psi + delta * Psi_d)
            + Deltab_dt * delta * Psi
        )

    phiR = sum1 + sum2 + sum3 + sum4
    phiR_d = sum1_d + sum2_d + sum3_d + sum4_d
    phiR_t = sum1_t + sum2_t + sum3_t + sum4_t
    phiR_dd = sum1_dd + sum2_dd + sum3_dd + sum4_dd
    phiR_tt = sum1_tt + sum2_tt + sum3_tt + sum4_tt
    phiR_dt = sum1_dt + sum2_dt + sum3_dt + sum4_dt

    return phi0, phi0_t, phi0_tt, phiR, phiR_d, phiR_t, phiR_dd, phiR_tt, phiR_dt


def phase_props(rho: float, T: float, Tc: float = TcW, rhoc: float = rho_cW, Rg: float = R_95):
    """
    Returns:
      cp [J/kg/K], cv [J/kg/K], cp/cv [-], (drho/dP)_T [kg/m^3/Pa]
    """
    rho = float(rho)
    T = float(T)
    tau = Tc / T
    delta = rho / rhoc

    phi0, _, phi0_tt, _, phiR_d, _, phiR_dd, phiR_tt, phiR_dt = helmholtz_all(tau, delta)

    cv = -Rg * (tau ** 2) * (phi0_tt + phiR_tt)

    num = (1.0 + delta * phiR_d - delta * tau * phiR_dt) ** 2
    den = (1.0 + 2.0 * delta * phiR_d + (delta ** 2) * phiR_dd)
    cp = Rg * (-(tau ** 2) * (phi0_tt + phiR_tt) + num / den)

    cp_cv = cp / cv

    dpdrho = Rg * T * (1.0 + 2.0 * delta * phiR_d + (delta ** 2) * phiR_dd)
    drhodP_T = 1.0 / dpdrho

    return cp, cv, cp_cv, drhodP_T


def drho_industrial(d: float, rhoc: float = rho_cW, Pc: float = PcW) -> float:
    """
    Industrial-region reference derivative used in the critical enhancement terms.

    Returns drho_ref = (drho/dP)_T,ref in [kg/m^3/Pa]
    """
    d = float(d)
    if d <= 0.310559006:
        ai = [6.53786807199516, -5.61149954923348, 3.39624167361325, -2.27492629730878, 10.2631854662709, 1.97815050331519]
    elif d <= 0.776397516:
        ai = [6.52717759281799, -6.30816983387575, 8.08379285492595, -9.82240510197603, 12.1358413791395, -5.54349664571295]
    elif d <= 1.242236025:
        ai = [5.35500529896124, -3.96415689925446, 8.91990208918795, -12.0338729505790, 9.19494865194302, -2.16866274479712]
    elif d <= 1.863354037:
        ai = [1.55225959906681, 0.464621290821181, 8.93237374861479, -11.0321960061126, 6.16780999933360, -0.965458722086812]
    else:
        ai = [1.11999926419994, 0.595748562571649, 9.88952565078920, -10.3255051147040, 4.66861294457414, -0.503243546373828]

    den = 0.0
    for i in range(0, 6):
        den += ai[i] * (d ** i)
    return (1.0 / den) * rhoc / Pc


def viscosity(rho: float, T: float, Tc: float = TcW, rhoc: float = rho_cW, Pc: float = PcW, Rg: float = R_95) -> float:
    """
    Dynamic viscosity mu [Pa·s] for water.
    """
    rho = float(rho)
    T = float(T)

    Tr = T / Tc
    Dr = rho / rhoc

    # mu0 in microPa*s
    H = np.array([1.67752, 2.20462, 0.6366564, -0.241605], dtype=float)
    mu0den = 0.0
    for i in range(0, 4):
        mu0den += H[i] / (Tr ** i)
    mu0 = 100.0 * np.sqrt(Tr) / mu0den

    # mu1 (finite-density)
    li = np.array([0, 1, 2, 3, 0, 1, 2, 3, 5, 0, 1, 2, 3, 4, 0, 1, 0, 3, 4, 3, 5], dtype=int)
    lj = np.array([0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 4, 4, 5, 6, 6], dtype=int)
    Hij = np.array([
        0.520094, 0.0850895, -1.08374, -0.289555, 0.222531, 0.999115,
        1.88797, 1.26613, 0.120573, -0.281378, -0.906851, -0.772479,
        -0.489837, -0.257040, 0.161913, 0.257399, -0.0325372, 0.0698452,
        0.00872102, -0.00435673, -0.000593264
    ], dtype=float)

    s = 0.0
    for k in range(Hij.size):
        i = li[k]
        j = lj[k]
        h = Hij[k]
        s += ((1.0 / Tr) - 1.0) ** i * h * (Dr - 1.0) ** j
    mu1 = np.exp(Dr * s)

    # mu2 (critical enhancement)
    mu2 = 1.0
    if rho > 0.0:
        tau = Tc / T
        delta = rho / rhoc
        _, _, _, _, phiR_d, _, phiR_dd, _, _ = helmholtz_all(tau, delta)

        dpdrho = Rg * T * (1.0 + 2.0 * delta * phiR_d + (delta ** 2) * phiR_dd)
        drhodP_T = 1.0 / dpdrho

        drho_ref = drho_industrial(delta, rhoc=rhoc, Pc=Pc)

        TR_over_T = 1.5 / Tr
        DeltaChi = delta * ((Pc / rhoc) * drhodP_T - (Pc / rhoc) * drho_ref * TR_over_T)
        if DeltaChi < 0.0:
            DeltaChi = 0.0

        xi_nm = 0.13 * (DeltaChi / 0.06) ** (0.63 / 1.239)

        qc = 1.0 / 1.9
        qd = 1.0 / 1.1

        if xi_nm <= 0.3817016416:
            Y = (qc / 5.0) * xi_nm * (qd * xi_nm) ** 5 * (
                1.0 - qc * xi_nm + (qc * xi_nm) ** 2 - (765.0 / 504.0) * (qd * xi_nm) ** 2
            )
        else:
            # Fid = acos((1+qd^2*xi_nm^2)^(-0.5))
            arg = (1.0 + (qd ** 2) * (xi_nm ** 2)) ** (-0.5)
            arg = float(np.clip(arg, -1.0, 1.0))
            Fid = np.arccos(arg)
            w = np.sqrt(abs((qc * xi_nm - 1.0) / (qc * xi_nm + 1.0))) * np.tan(Fid / 2.0)
            if qc * xi_nm > 1.0:
                Lw = np.log((1.0 + w) / (1.0 - w))
            else:
                Lw = 2.0 * np.arctan(abs(w))
            Y = (
                np.sin(3.0 * Fid) / 12.0
                - np.sin(2.0 * Fid) / (4.0 * qc * xi_nm)
                + (1.0 - 5.0 / 4.0 * (qc * xi_nm) ** 2) / (qc * xi_nm) ** 2 * np.sin(Fid)
                - (
                    (1.0 - 3.0 / 2.0 * (qc * xi_nm) ** 2) * Fid
                    - abs((qc * xi_nm) ** 2 - 1.0) ** (3.0 / 2.0) * Lw
                )
                / (qc * xi_nm) ** 3
            )

        mu2 = np.exp(0.068 * Y)

    # microPa*s -> Pa*s
    return float(mu0 * mu1 * mu2 * 1e-6)


def thermal_conductivity(rho: float, T: float, Tc: float = TcW, rhoc: float = rho_cW, Pc: float = PcW,
                         Rg: float = R_95, return_parts: bool = False):
    """
    Thermal conductivity k [W/m/K] for water.

    If return_parts=True returns (k, lambda0, lambda1, lambda2) where:
      lambda0 and lambda2 are in W/m/K, lambda1 is dimensionless.
    """
    rho = float(rho)
    T = float(T)

    d = rho / rhoc
    Tr = T / Tc

    # lambda0 in mW/mK
    no = np.array([2.443221e-3, 1.323095e-2, 6.770357e-3, -3.454586e-3, 4.096266e-4], dtype=float)
    k0den = 0.0
    for i in range(0, 5):
        k0den += no[i] / (Tr ** i)
    lambda0_mWmK = np.sqrt(Tr) / k0den

    # lambda1 dimensionless
    li = np.array([0, 0, 0, 0, 0, 0,
                   1, 1, 1, 1, 1, 1,
                   2, 2, 2, 2, 2, 2,
                   3, 3, 3, 3,
                   4, 4, 4, 4, 4, 4], dtype=int)
    lj = np.array([0, 1, 2, 3, 4, 5,
                   0, 1, 2, 3, 4, 5,
                   0, 1, 2, 3, 4, 5,
                   0, 1, 2, 3,
                   0, 1, 2, 3, 4, 5], dtype=int)
    nij = np.array([
        1.60397357, -0.646013523, 0.111443906, 0.102997357, -0.0504123634, 0.00609859258,
        2.33771842, -2.78843778, 1.53616167, -0.463045512, 0.0832827019, -0.00719201245,
        2.19650529, -4.54580785, 3.55777244, -1.40944978, 0.275418278, -0.0205938816,
        -1.21051378, 1.60812989, -0.621178141, 0.0716373224,
        -2.7203370, 4.57586331, -3.18369245, 1.1168348, -0.19268305, 0.012913842
    ], dtype=float)

    s = 0.0
    for m in range(nij.size):
        i = li[m]
        j = lj[m]
        n = nij[m]
        s += ((1.0 / Tr) - 1.0) ** i * n * (d - 1.0) ** j
    lambda1 = np.exp(d * s)

    # lambda2 (critical enhancement) in mW/mK
    lambda2_mWmK = 0.0
    if rho > 0.0:
        cp, cv, cp_cv, drhodP_T = phase_props(rho, T, Tc=Tc, rhoc=rhoc, Rg=Rg)
        mu = viscosity(rho, T, Tc=Tc, rhoc=rhoc, Pc=Pc, Rg=Rg)
        drho_ref = drho_industrial(d, rhoc=rhoc, Pc=Pc)

        TR_over_T = 1.5 / Tr
        DeltaChi = d * ((Pc / rhoc) * drhodP_T - (Pc / rhoc) * drho_ref * TR_over_T)
        if DeltaChi < 0.0:
            DeltaChi = 0.0

        xi_nm = 0.13 * (DeltaChi / 0.06) ** (0.63 / 1.239)
        y = xi_nm / 0.4

        if y < 1.2e-7:
            Z = 0.0
        else:
            # guard d near zero
            d_safe = max(d, 1e-30)
            Z = (2.0 / np.pi / y) * (
                ((1.0 - 1.0 / cp_cv) * np.arctan(y) + y / cp_cv)
                - (1.0 - np.exp(-1.0 / (1.0 / y + y ** 2 / (3.0 * d_safe ** 2))))
            )

        # 177.8514 factor matches the MATLAB transcription; keep as-is
        lambda2_mWmK = 177.8514 * d * (cp / Rg) * Tr / mu * 1e-6 * Z

    lambda_mWmK = lambda0_mWmK * lambda1 + lambda2_mWmK

    # convert to SI
    k = float(1e-3 * lambda_mWmK)
    lambda0 = float(1e-3 * lambda0_mWmK)
    lambda2 = float(1e-3 * lambda2_mWmK)

    if return_parts:
        return k, lambda0, float(lambda1), lambda2
    return k


def cp_mass(rho: float, T: float) -> float:
    """Convenience wrapper for cp in J/kg/K."""
    cp, _, _, _ = phase_props(rho, T)
    return float(cp)


def cp_molar(rho: float, T: float) -> float:
    """Convenience wrapper for cp in J/mol/K."""
    return cp_mass(rho, T) * M_W


# -----------------------
# Optional self-test when run as a script
# -----------------------

# -----------------------
# State functions (p, h) and rho(T,p) inversion
# -----------------------
def pressure(rho: float, T: float) -> float:
    """Pressure p [Pa] from IAPWS-95 Helmholtz formulation at (rho [kg/m^3], T [K])."""
    if rho <= 0.0 or T <= 0.0:
        raise ValueError("pressure: rho and T must be positive.")
    tau = TcW / T
    delta = rho / rho_cW
    # helmholtz_all returns (phi0, phi0_t, phi0_tt, phiR, phiR_d, phiR_t, phiR_dd, phiR_tt, phiR_dt)
    _, _, _, _, phiR_d, _, _, _, _ = helmholtz_all(tau, delta)
    return rho * R_95 * T * (1.0 + delta * phiR_d)

def dpdrho_T(rho: float, T: float) -> float:
    """(∂p/∂rho)_T [Pa·m^3/kg] from IAPWS-95 at (rho,T)."""
    if rho <= 0.0 or T <= 0.0:
        raise ValueError("dpdrho_T: rho and T must be positive.")
    tau = TcW / T
    delta = rho / rho_cW
    _, _, _, _, phiR_d, _, phiR_dd, _, _ = helmholtz_all(tau, delta)
    return R_95 * T * (1.0 + 2.0 * delta * phiR_d + (delta ** 2) * phiR_dd)

def enthalpy_mass(rho: float, T: float) -> float:
    """Specific enthalpy h [J/kg] from IAPWS-95 at (rho,T)."""
    if rho <= 0.0 or T <= 0.0:
        raise ValueError("enthalpy_mass: rho and T must be positive.")
    tau = TcW / T
    delta = rho / rho_cW
    _, phi0_t, _, _, phiR_d, phiR_t, _, _, _ = helmholtz_all(tau, delta)
    return R_95 * T * (1.0 + tau * (phi0_t + phiR_t) + delta * phiR_d)

def entropy_mass(rho: float, T: float) -> float:
    """Specific entropy s [J/kg/K] from IAPWS-95 at (rho,T)."""
    if rho <= 0.0 or T <= 0.0:
        raise ValueError("entropy_mass: rho and T must be positive.")
    tau = TcW / T
    delta = rho / rho_cW
    phi0, phi0_t, _, phiR, _, phiR_t, _, _, _ = helmholtz_all(tau, delta)
    return R_95 * (tau * (phi0_t + phiR_t) - phi0 - phiR)

def enthalpy_molar(rho: float, T: float) -> float:
    """Molar enthalpy h [J/mol] from IAPWS-95 at (rho,T)."""
    return enthalpy_mass(rho, T) * M_W

def entropy_molar(rho: float, T: float) -> float:
    """Molar entropy s [J/mol/K] from IAPWS-95 at (rho,T)."""
    return entropy_mass(rho, T) * M_W

def rho_from_Tp_newton(T: float,
                       p_target: float,
                       rho0: float | None = None,
                       *,
                       max_iter: int = 50,
                       rtol_p: float = 1e-10,
                       atol_p: float = 1e-2) -> float:
    """
    Solve rho [kg/m^3] for given (T [K], p_target [Pa]) using a damped Newton method
    on the IAPWS-95 pressure equation. Hard-fails on nonconvergence.

    Notes:
    - For subcritical liquid water, providing rho0 ~ 800-1200 greatly improves robustness.
    - This routine does not attempt multiphase root selection; callers should seed rho0
      on the desired branch (liquid vs vapor).
    """
    if T <= 0.0:
        raise ValueError("rho_from_Tp_newton: T must be positive.")
    if p_target <= 0.0:
        raise ValueError("rho_from_Tp_newton: p_target must be positive.")
    # default initial guess: liquid-like if below critical, else ideal-gas
    if rho0 is None:
        if T < TcW:
            rho = 1000.0
        else:
            rho = max(p_target / (R_95 * T), 0.1)
    else:
        rho = float(rho0)

    # Keep rho in a sane positive range
    rho = max(rho, 1e-6)

    for it in range(max_iter):
        p_calc = pressure(rho, T)
        f = p_calc - p_target
        # convergence checks
        if abs(f) <= atol_p or abs(f) <= rtol_p * max(p_target, 1.0):
            return float(rho)

        dfd = dpdrho_T(rho, T)
        if (not np.isfinite(dfd)) or dfd <= 0.0:
            raise RuntimeError(f"rho_from_Tp_newton: nonpositive/invalid dpdrho_T at iter {it}: {dfd}")

        step = -f / dfd

        # Damping / step limiting for robustness
        # Limit relative step to 50% of current rho magnitude.
        max_step = 0.5 * rho
        if abs(step) > max_step:
            step = np.sign(step) * max_step

        rho_new = rho + step
        # Ensure rho stays positive; backtrack if needed
        if rho_new <= 0.0 or (not np.isfinite(rho_new)):
            # halve step until positive
            step_bt = step
            for _ in range(20):
                step_bt *= 0.5
                rho_new = rho + step_bt
                if rho_new > 0.0 and np.isfinite(rho_new):
                    break
            else:
                raise RuntimeError("rho_from_Tp_newton: failed to keep rho positive during damping.")

        rho = rho_new

    raise RuntimeError(f"rho_from_Tp_newton: did not converge in {max_iter} iterations at T={T}, p={p_target}.")


@lru_cache(maxsize=4096)
def saturation_pressure_ancillary_Pa(T: float) -> float:
    """
    Saturation pressure p_sat [Pa] from the IAPWS saturation-pressure ancillary.

    This is the standard ancillary relation used around the IAPWS-95 Helmholtz core,
    and is the project entry point for water-vapor psychrometric calculations.
    """
    T = float(T)
    if T <= 0.0:
        raise ValueError("saturation_pressure_ancillary_Pa: T must be positive.")
    tau = 1.0 - T / TcW
    a = (-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502)
    b = (1.0, 1.5, 3.0, 3.5, 4.0, 7.5)
    series = sum(ai * tau ** bi for ai, bi in zip(a, b))
    return float(PcW * np.exp((TcW / T) * series))


def _validated_subcritical_saturation_temperature(T: float) -> float:
    T = float(T)
    if T <= 0.0:
        raise ValueError("Saturation properties require T > 0 K.")
    if T >= TcW:
        raise ValueError(f"Saturation properties are only defined below the critical temperature {TcW} K.")
    return T


@lru_cache(maxsize=4096)
def saturated_liquid_density_kg_m3(T: float) -> float:
    """Saturated-liquid density rho_l,sat [kg/m^3] from Wagner-Pruss Eq. (2.6)."""
    T = _validated_subcritical_saturation_temperature(T)
    tau = 1.0 - T / TcW
    b = (
        1.99274064,
        1.09965342,
        -0.510839303,
        -1.75493479,
        -45.5170352,
        -6.74694450e5,
    )
    exponents = (1.0 / 3.0, 2.0 / 3.0, 5.0 / 3.0, 16.0 / 3.0, 43.0 / 3.0, 110.0 / 3.0)
    series = 1.0 + sum(bi * tau ** ei for bi, ei in zip(b, exponents))
    return float(rho_cW * series)


@lru_cache(maxsize=4096)
def saturated_vapor_density_kg_m3(T: float) -> float:
    """Saturated-vapor density rho_v,sat [kg/m^3] from Wagner-Pruss Eq. (2.7)."""
    T = _validated_subcritical_saturation_temperature(T)
    tau = 1.0 - T / TcW
    c = (-2.03150240, -2.68302940, -5.38626492, -17.2991605, -44.7586581, -63.9201063)
    exponents = (2.0 / 6.0, 4.0 / 6.0, 8.0 / 6.0, 18.0 / 6.0, 37.0 / 6.0, 71.0 / 6.0)
    series = sum(ci * tau ** ei for ci, ei in zip(c, exponents))
    return float(rho_cW * np.exp(series))


@lru_cache(maxsize=4096)
def saturated_liquid_enthalpy_J_kg(T: float) -> float:
    """Saturated-liquid specific enthalpy h_l,sat [J/kg] from IAPWS-95."""
    T = _validated_subcritical_saturation_temperature(T)
    rho_l = saturated_liquid_density_kg_m3(T)
    return float(enthalpy_mass(rho_l, T))


@lru_cache(maxsize=4096)
def saturated_vapor_enthalpy_J_kg(T: float) -> float:
    """Saturated-vapor specific enthalpy h_v,sat [J/kg] from IAPWS-95."""
    T = _validated_subcritical_saturation_temperature(T)
    rho_v = saturated_vapor_density_kg_m3(T)
    return float(enthalpy_mass(rho_v, T))


@lru_cache(maxsize=4096)
def latent_heat_vaporization_saturated_J_kg(T: float) -> float:
    """Latent heat h_fg [J/kg] along saturation from IAPWS-95 saturated enthalpy difference."""
    T = _validated_subcritical_saturation_temperature(T)
    return float(saturated_vapor_enthalpy_J_kg(T) - saturated_liquid_enthalpy_J_kg(T))


def _self_test():
    # Example checks (same as MATLAB snippet)
    rho_test = 998.0
    T_test = 298.15
    k = thermal_conductivity(rho_test, T_test)
    mu = viscosity(rho_test, T_test)
    print("Example checks:")
    print(f"k(998,298.15) [W/mK] = {k:.6f}")
    print(f"mu(998,298.15) [Pa*s] = {mu:.9e}")

    rho_test = 0.0
    T_test = 873.15
    k = thermal_conductivity(rho_test, T_test)
    print(f"k(0,873.15) [W/mK] = {k:.6f}")

    # Validation: Table 4 (lambda2=0)
    T_ref_v4 = np.array([298.15, 298.15, 298.15, 873.15], dtype=float)
    rho_ref_v4 = np.array([0.0, 998.0, 1200.0, 0.0], dtype=float)
    lambda_ref_v4_mWmK = np.array([18.434, 607.71, 799.04, 79.103], dtype=float)

    lambda_calc_v4_mWmK = np.zeros_like(lambda_ref_v4_mWmK)
    for i in range(T_ref_v4.size):
        k_i = thermal_conductivity(float(rho_ref_v4[i]), float(T_ref_v4[i]))
        lambda_calc_v4_mWmK[i] = k_i * 1e3

    abs_err = lambda_calc_v4_mWmK - lambda_ref_v4_mWmK
    rel_err = 100.0 * abs_err / lambda_ref_v4_mWmK
    print("\nIAPWS ThCond verification (Table 4, lambda2 = 0):")
    for i in range(T_ref_v4.size):
        print(f"T={T_ref_v4[i]:8.2f} K  rho={rho_ref_v4[i]:8.1f}  ref={lambda_ref_v4_mWmK[i]:9.3f}  "
              f"calc={lambda_calc_v4_mWmK[i]:9.3f}  abs_err={abs_err[i]:9.3f}  rel_err%={rel_err[i]:9.3f}")

    # Validation: Table 5
    T_ref_v5 = np.full(8, 647.35, dtype=float)
    rho_ref_v5 = np.array([1, 122, 222, 272, 322, 372, 422, 750], dtype=float)

    lambda0_ref = np.full(8, 51.576, dtype=float)
    lambda1_ref = np.array([1.0068, 2.1445, 3.4841, 4.2234, 4.9682, 5.6961, 6.3973, 11.587], dtype=float)
    lambda2_ref = np.array([0.00013, 20.316, 188.09, 540.13, 1187.5, 356.53, 118.93, 3.3419], dtype=float)
    lambda_ref = np.array([51.93, 130.92, 367.79, 757.96, 1443.8, 650.32, 448.88, 600.96], dtype=float)

    l0_calc = np.zeros(8); l1_calc = np.zeros(8); l2_calc = np.zeros(8); l_calc = np.zeros(8)
    for i in range(8):
        k_i, l0_i, l1_i, l2_i = thermal_conductivity(float(rho_ref_v5[i]), float(T_ref_v5[i]), return_parts=True)
        l0_calc[i] = l0_i * 1e3
        l1_calc[i] = l1_i
        l2_calc[i] = l2_i * 1e3
        l_calc[i] = k_i * 1e3

    print("\nIAPWS ThCond verification (Table 5):")
    for i in range(8):
        print(f"rho={rho_ref_v5[i]:6.1f}  l0(ref/calc)={lambda0_ref[i]:8.3f}/{l0_calc[i]:8.3f}  "
              f"l1(ref/calc)={lambda1_ref[i]:7.4f}/{l1_calc[i]:7.4f}  "
              f"l2(ref/calc)={lambda2_ref[i]:8.3f}/{l2_calc[i]:8.3f}  "
              f"lam(ref/calc)={lambda_ref[i]:8.2f}/{l_calc[i]:8.2f}")

if __name__ == "__main__":
    _self_test()
