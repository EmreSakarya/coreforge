"""CoreForge thermal-hydraulics — steady-state closed-channel model.

Single-phase PWR subchannel analysis for the AVERAGE and the HOT
assembly, driven by the neutronics solution (radial peaking from the
block-power map, axial shape from the 3-D profile or a chopped cosine):

  * coolant temperature rise:  m_ch cp dT/dz = q'(z) n_pins
  * clad surface:  T_clad = T_cool + q'' / h,  h from Dittus-Boelter
        Nu = 0.023 Re^0.8 Pr^0.4
  * fuel:  ΔT_pellet = q'/(4π k_f)   (centerline, uniform-source rod)
           ΔT_gap    = q''_pellet / h_gap
           ΔT_clad   = q' ln(r_o/r_i)/(2π k_c)
  * pressure drop:  ΔP = f (H/D_h) ρ v²/2,  f = 0.184 Re^-0.2

Water properties are evaluated at a representative PWR condition
(~15.5 MPa, ~310 °C) and held constant (documented single-phase
approximation; boiling is flagged via the saturation-margin check,
not modelled).  Energy conservation is exact by construction and is
asserted in verify.py together with hand-calculable checks of the
conduction and convection formulas.
"""
import math

import numpy as np

# representative PWR water properties (15.5 MPa, ~310 C) — documented
RHO_W = 704.0        # kg/m3
CP_W = 5750.0        # J/kg-K
MU_W = 8.5e-5        # Pa-s
K_W = 0.54           # W/m-K
PR_W = MU_W * CP_W / K_W

# water saturation temperature vs pressure (IAPWS-IF97 points, degC)
_PSAT_MPa = [5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 15.5, 16.0, 17.0, 18.0]
_TSAT_C = [263.9, 285.8, 303.3, 318.1, 330.9, 342.1, 344.8, 347.4, 352.3, 357.0]

# a symmetric chopped-cosine axial shape has a hard peak/mean ceiling of
# pi/2; higher axial peaking needs an asymmetric shape (use the real 3-D
# axial profile instead).  Exposed so the UI can cap the input honestly.
FZ_MAX_COSINE = math.pi / 2.0


def t_sat(pressure_MPa):
    """Water saturation temperature [degC] at the given pressure."""
    import numpy as _np
    return float(_np.interp(pressure_MPa, _PSAT_MPa, _TSAT_C))

# fuel-rod thermal constants (UO2 / Zircaloy, hot, documented)
K_FUEL = 3.0         # W/m-K   (UO2, ~1000 K effective)
K_CLAD = 17.0        # W/m-K   (Zircaloy)
H_GAP = 6000.0       # W/m2-K  (He gap, closed-gap era value)
T_SAT_155 = 344.8    # C, saturation at 15.5 MPa


def _chopped_cosine(nz, fz):
    """Symmetric chopped-cosine axial shape (mean 1) with the requested
    peak/mean = fz.  fz is CLAMPED to the cosine ceiling pi/2 (a symmetric
    cosine cannot be more peaked); the caller is told via _cosine_clamped.
    """
    fz_eff = min(fz, FZ_MAX_COSINE - 1e-6)
    z = (np.arange(nz) + 0.5) / nz
    lo, hi = 1e-3, 0.999999
    for _ in range(60):
        c = 0.5 * (lo + hi)
        mean = np.sin(np.pi * c / 2) / (np.pi * c / 2)
        if 1.0 / mean < fz_eff:
            lo = c
        else:
            hi = c
    c = 0.5 * (lo + hi)
    s = np.cos(np.pi * c * (z - 0.5))
    return s / s.mean()


def channel(power_MW_assembly, H_m, n_pins=264, pin_d_mm=9.5,
            pin_pitch_mm=12.6, mflow_kg_s=90.0, t_in_C=292.0,
            shape=None, nz=40, fz=1.40, pressure_MPa=15.5):
    """Solve one closed channel (one assembly).  Returns dict of axial
    arrays (z, T_cool, T_clad, T_fuel_CL, q_lin) and scalar results."""
    fz_clamped = shape is None and fz > FZ_MAX_COSINE
    if shape is None:
        shape = _chopped_cosine(nz, fz)
    nz = len(shape)
    H = H_m
    dz = H / nz
    z = (np.arange(nz) + 0.5) * dz

    d = pin_d_mm * 1e-3
    p = pin_pitch_mm * 1e-3
    r_o = d / 2
    r_i = r_o - 0.6e-3                     # clad thickness 0.6 mm
    r_f = r_i - 0.08e-3                    # gap 80 um
    A_flow = n_pins * (p * p - np.pi * r_o**2)
    P_wet = n_pins * np.pi * d
    D_h = 4 * A_flow / P_wet
    v = mflow_kg_s / (RHO_W * A_flow)
    Re = RHO_W * v * D_h / MU_W
    h = 0.023 * Re**0.8 * PR_W**0.4 * K_W / D_h

    q_lin_avg = power_MW_assembly * 1e6 / (n_pins * H)   # W/m per pin
    q_lin = q_lin_avg * shape

    # coolant temperature (exact energy integration)
    dT = q_lin * n_pins * dz / (mflow_kg_s * CP_W)
    T_cool = t_in_C + np.cumsum(dT) - dT / 2

    q2_clad = q_lin / (np.pi * d)                        # W/m2 at clad OD
    T_clad = T_cool + q2_clad / h
    dT_clad = q_lin * np.log(r_o / r_i) / (2 * np.pi * K_CLAD)
    q2_pel = q_lin / (2 * np.pi * r_f)
    dT_gap = q2_pel / H_GAP
    dT_fuel = q_lin / (4 * np.pi * K_FUEL)
    T_fuel = T_clad + dT_clad + dT_gap + dT_fuel

    f = 0.184 * Re**-0.2
    dP = f * (H / D_h) * RHO_W * v**2 / 2                # Pa (friction)

    return dict(z=z, T_cool=T_cool, T_clad=T_clad, T_fuel=T_fuel,
                q_lin=q_lin, T_out=float(t_in_C + dT.sum()),
                fz_clamped=bool(fz_clamped), fz_max=float(FZ_MAX_COSINE),
                p_MPa=float(pressure_MPa),
                T_clad_max=float(T_clad.max()),
                T_fuel_max=float(T_fuel.max()),
                q_lin_peak=float(q_lin.max()),
                v_cool=float(v), Re=float(Re), htc=float(h),
                dP_kPa=float(dP / 1e3),
                sat_margin=float(t_sat(pressure_MPa) - T_clad.max()),
                t_sat=float(t_sat(pressure_MPa)))


def core_th(P_core_MW, n_fuel_assemblies, F_dh, H_m,
            axial_profile=None, mflow_total_kg_s=None, t_in_C=292.0,
            n_pins=264, pin_d_mm=9.5, pin_pitch_mm=12.6, fz=1.40,
            pressure_MPa=15.5):
    """Average- and hot-channel analysis of the whole core.
    F_dh: hot-assembly radial peaking (from the block-power map)."""
    if mflow_total_kg_s is None:
        mflow_total_kg_s = 90.0 * n_fuel_assemblies
    m_ch = mflow_total_kg_s / n_fuel_assemblies
    P_avg = P_core_MW / n_fuel_assemblies
    shape = (None if axial_profile is None
             else np.asarray(axial_profile, dtype=float))
    if shape is not None:
        shape = shape[shape > 1e-12]
        shape = shape / shape.mean()
    avg = channel(P_avg, H_m, n_pins, pin_d_mm, pin_pitch_mm,
                  m_ch, t_in_C, shape=shape, fz=fz, pressure_MPa=pressure_MPa)
    hot = channel(P_avg * F_dh, H_m, n_pins, pin_d_mm, pin_pitch_mm,
                  m_ch, t_in_C, shape=shape, fz=fz, pressure_MPa=pressure_MPa)
    # exact core energy balance (verification quantity)
    Q_balance_MW = mflow_total_kg_s * CP_W * (avg["T_out"] - t_in_C) / 1e6
    return dict(avg=avg, hot=hot, m_channel=m_ch,
                Q_balance_MW=float(Q_balance_MW),
                F_dh=F_dh, t_in=t_in_C, H=H_m,
                n_assemblies=n_fuel_assemblies)
